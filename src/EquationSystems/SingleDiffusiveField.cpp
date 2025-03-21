#include "SingleDiffusiveField.hpp"

namespace NESO::Solvers::tokamak
{

/// Name of class
static std::string class_name;
std::string SingleDiffusiveField::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "SingleDiffusiveField", SingleDiffusiveField::create,
        "Solves for a single diffusive field (n) with anisotropy");
/**
 * @brief Creates an instance of this class.
 */
static SU::EquationSystemSharedPtr create(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
{
    SU::EquationSystemSharedPtr p =
        MemoryManager<SingleDiffusiveField>::AllocateSharedPtr(session, graph);
    p->InitObject();
    return p;
}

SingleDiffusiveField::SingleDiffusiveField(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
    : TokamakSystem(session, graph)
{
    this->required_fld_names = {"n"};

    if (this->particles_enabled)
    {
        for (auto &[k, v] : particle_sys->get_species())
        {
            this->required_fld_names.push_back(v.name + "_src");
        }
    }
}

void SingleDiffusiveField::v_InitObject(bool DeclareFields)
{
    TokamakSystem::v_InitObject(DeclareFields);

    m_session->MatchSolverInfo("SpectralVanishingViscosity", "True",
                               m_useSpecVanVisc, false);
    m_session->LoadParameter("epsilon", m_epsilon, 1.0);
    if (m_useSpecVanVisc)
    {
        m_session->LoadParameter("SVVCutoffRatio", m_sVVCutoffRatio, 0.75);
        m_session->LoadParameter("SVVDiffCoeff", m_sVVDiffCoeff, 0.1);
        m_factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
        m_factors[StdRegions::eFactorSVVDiffCoeff] = m_sVVDiffCoeff / m_epsilon;
    }
    CalcDiffTensor();

    // Setup diffusion object
    m_ode.DefineOdeRhs(&SingleDiffusiveField::DoOdeRhs, this);

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            std::string diffName;
            m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
            m_diffusion = SolverUtils::GetDiffusionFactory().CreateInstance(
                diffName, diffName);
            m_diffusion->SetFluxVector(&SingleDiffusiveField::GetFluxVectorDiff,
                                       this);
            m_diffusion->InitObject(m_session, m_fields);
            break;
        }
        case MultiRegions::eGalerkin:
        {
            // Enable implicit solver for CG.
            m_ode.DefineImplicitSolve(&SingleDiffusiveField::ImplicitTimeIntCG,
                                      this);

            break;
        }
        default:
        {
            ASSERTL0(false, "Unknown projection scheme");
            break;
        }
    }
    if (this->particles_enabled)
    {
        for (auto &[k, v] : this->particle_sys->get_species())
        {
            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            this->src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_DENSITY"));
            this->components.push_back(0);
        }

        this->particle_sys->setup_project(this->src_fields);
    }
}
void SingleDiffusiveField::ImplicitTimeIntCG(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble time, const NekDouble lambda)
{

    int nvariables                       = inarray.size();
    int npoints                          = m_fields[0]->GetNpoints();
    m_factors[StdRegions::eFactorLambda] = 1.0 / lambda / m_epsilon;

    if (m_useSpecVanVisc)
    {
        m_factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
        m_factors[StdRegions::eFactorSVVDiffCoeff] = m_sVVDiffCoeff / m_epsilon;
    }

    // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
    // inarray = input: \hat{rhs} -> output: \hat{Y}
    // outarray = output: nabla^2 \hat{Y}
    // where \hat = modal coeffs
    if (m_intScheme->GetIntegrationSchemeType() == LibUtilities::eImplicit)
    {
        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Zero(npoints, outarray[i], 1);
        }

        for (auto &x : m_forcing)
        {
            x->Apply(m_fields, inarray, outarray, time);
        }
    }
    CalcDiffTensor();
    for (int i = 0; i < nvariables; ++i)
    {
        if (m_intScheme->GetIntegrationSchemeType() == LibUtilities::eImplicit)
        {
            // Multiply forcing term by -1 for definition of HelmSolve function
            Vmath::Smul(npoints, -1.0, outarray[i], 1, outarray[i], 1);

            // Multiply 1.0/timestep/lambda
            Vmath::Svtvp(npoints, -m_factors[StdRegions::eFactorLambda],
                         inarray[i], 1, outarray[i], 1, outarray[i], 1);
        }
        else
        {
            // Multiply 1.0/timestep/lambda
            Vmath::Smul(npoints, -m_factors[StdRegions::eFactorLambda],
                        inarray[i], 1, outarray[i], 1);
        }

        // Solve a system of equations with Helmholtz solver
        m_fields[i]->HelmSolve(outarray[i], m_fields[i]->UpdateCoeffs(),
                               m_factors, m_D);

        m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(), outarray[i]);

        m_fields[i]->SetPhysState(false);
    }
}

void SingleDiffusiveField::CalcKPar()
{
    // Change to fn of fields
    int npoints     = m_fields[0]->GetNpoints();
    NekDouble k_par = this->k_par;
    m_kpar          = Array<OneD, NekDouble>(npoints, k_par);
}

void SingleDiffusiveField::CalcKPerp()
{
    // Change to fn of fields
    int npoints      = m_fields[0]->GetNpoints();
    NekDouble k_perp = this->k_perp;
    m_session->LoadParameter("k_perp", k_perp, 1.0);
    m_kperp = Array<OneD, NekDouble>(npoints, k_perp);
}

void SingleDiffusiveField::CalcDiffTensor()
{
    int npoints = m_fields[0]->GetNpoints();
    CalcKPar();
    CalcKPerp();
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Array<OneD, NekDouble> d(npoints, 0.0);
            for (int k = 0; k < npoints; k++)
            {
                d[k] = (m_kpar[k] - m_kperp[k]) * b_unit[i][k] * b_unit[j][k];
                if (i == j)
                {
                    d[k] += m_kperp[k];
                }
            }
            m_D[vc[i][j]] = d;
        }
    }
}

/**
 * @brief Populate rhs array ( @p out_arr ) with diffused quantitites
 * @param in_arr physical values of all fields
 * @param[out] out_arr output array (RHSs of time integration equations)
 */
void SingleDiffusiveField::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time)
{
    if (m_explicitDiffusion)
    {
        CalcDiffTensor();
        size_t nvariables = in_arr.size();
        m_diffusion->Diffuse(nvariables, m_fields, in_arr, out_arr);
    }
    else
    {
        // RHS should be set to zero.
        for (int i = 0; i < out_arr.size(); ++i)
        {
            Vmath::Zero(out_arr[i].size(), &out_arr[i][0], 1);
        }
    }

    if (this->particles_enabled)
    {
        for (int i = 0; i < this->particle_sys->get_species().size(); ++i)
        {
            Vmath::Vadd(out_arr[0].size(), out_arr[0], 1,
                        this->src_fields[i]->GetPhys(), 1, out_arr[0], 1);
        }
    }

    // Add forcing terms
    for (auto &x : m_forcing)
    {
        x->Apply(m_fields, in_arr, out_arr, time);
    }
}

/**
 * @brief Construct the flux vector for the anisotropic diffusion problem.
 */
void SingleDiffusiveField::GetFluxVectorDiff(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    unsigned int nDim = qfield.size();
    unsigned int nPts = qfield[0][0].size();
    int n_idx         = this->field_to_index["n"];
    for (unsigned int j = 0; j < nDim; ++j)
    {
        // Calc diffusion of n with D tensor
        Vmath::Vmul(nPts, m_D[vc[j][0]].GetValue(), 1, qfield[0][n_idx], 1,
                    fluxes[j][n_idx], 1);
        for (unsigned int k = 1; k < nDim; ++k)
        {
            Vmath::Vvtvp(nPts, m_D[vc[j][k]].GetValue(), 1, qfield[k][n_idx], 1,
                         fluxes[j][n_idx], 1, fluxes[j][n_idx], 1);
        }
    }
}

void SingleDiffusiveField::load_params()
{
    TokamakSystem::load_params();
    // Adiabatic gamma
    m_session->LoadParameter("k_B", this->m_k_B);
    m_session->LoadParameter("k_par", this->k_par, 100.0);
    m_session->LoadParameter("k_perp", this->k_perp, 1.0);
}

void SingleDiffusiveField::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    TokamakSystem::v_ExtraFldOutput(fieldcoeffs, variables);
    const int nPhys   = m_fields[0]->GetNpoints();
    const int nCoeffs = m_fields[0]->GetNcoeffs();
    int i             = 0;
    if (this->particles_enabled)
    {
        for (auto &[k, v] : this->particle_sys->get_species())
        {
            variables.push_back(v.name + "_SOURCE_DENSITY");
            Array<OneD, NekDouble> SrcFwd(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[i++]->GetPhys(),
                                           SrcFwd);
            fieldcoeffs.push_back(SrcFwd);
        }
    }
}
} // namespace NESO::Solvers::tokamak