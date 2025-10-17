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
    this->n_indep_fields       = 0;
    this->n_fields_per_species = 1;
    this->required_fld_names   = {"n"};
}

void SingleDiffusiveField::v_InitObject(bool DeclareFields)
{
    TokamakSystem::v_InitObject(DeclareFields);
    this->ne = std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]);

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
            m_diffusion->InitObject(m_session, m_indfields);
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
        std::vector<Sym<REAL>> src_syms;
        std::vector<int> src_components;
        for (auto &[k, v] : this->particle_sys->get_species())
        {
            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_DENSITY"));
            src_components.push_back(0);
        }

        this->particle_sys->setup_evaluate_fields(this->E, this->B, this->ne,
                                                  this->Te, this->ve);

        this->particle_sys->finish_setup(this->src_fields, src_syms,
                                         src_components);
    }
}

void SingleDiffusiveField::ImplicitTimeIntCG(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble time, const NekDouble lambda)
{

    int nvariables                       = inarray.size();
    int npoints                          = m_indfields[0]->GetNpoints();
    m_factors[StdRegions::eFactorLambda] = 1.0 / lambda / m_epsilon;

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
            x->Apply(m_indfields, inarray, outarray, time);
        }
    }
    for (int i = 0; i < nvariables; ++i)
    {
        CalcDiffTensor(i);
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
        m_indfields[i]->HelmSolve(outarray[i], m_indfields[i]->UpdateCoeffs(),
                                  m_factors, m_D);

        m_indfields[i]->BwdTrans(m_indfields[i]->GetCoeffs(), outarray[i]);

        m_indfields[i]->SetPhysState(false);
    }
}

void SingleDiffusiveField::CalcKPar(int f)
{
    int npoints = m_fields[0]->GetNpoints();
    double Z;
    this->neso_config->load_species_parameter(f, "Charge", Z);
    m_kpar = Array<OneD, NekDouble>(npoints, this->k_par / (Z * Z));
    Vmath::Vdiv(npoints, m_kpar, 1, m_indfields[f]->GetPhys(), 1, m_kpar, 1);
}

void SingleDiffusiveField::CalcKPerp(int f)
{
    int npoints = m_fields[0]->GetNpoints();
    double Z, A;
    this->neso_config->load_species_parameter(f, "Charge", Z);
    this->neso_config->load_species_parameter(f, "Mass", A);
    m_kperp =
        Array<OneD, NekDouble>(npoints, this->k_perp * Z * Z * std::sqrt(A));
    Vmath::Vmul(npoints, m_kpar, 1, m_indfields[f]->GetPhys(), 1, m_kpar, 1);
    Vmath::Vdiv(npoints, m_kpar, 1, this->mag_B, 1, m_kpar, 1);
}

void SingleDiffusiveField::CalcDiffTensor(int f)
{
    int npoints = m_fields[0]->GetNpoints();
    CalcKPar(f);
    CalcKPerp(f);
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
    if (m_explicitDiffusion || m_projectionType == MR::eDiscontinuous)
    {
        size_t nvariables = in_arr.size();
        m_diffusion->Diffuse(nvariables, m_indfields, in_arr, out_arr);
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
            Vmath::Vadd(out_arr[i].size(), out_arr[i], 1,
                        this->src_fields[i]->GetPhys(), 1, out_arr[i], 1);
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
    unsigned int nFld = qfield[0].size();
    unsigned int nPts = qfield[0][0].size();

    for (unsigned int f = 0; f < nFld; ++f)
    {
        CalcDiffTensor(f);
        for (unsigned int j = 0; j < nDim; ++j)
        {
            // Calc diffusion of n with D tensor
            Vmath::Vmul(nPts, m_D[vc[j][0]].GetValue(), 1, qfield[0][f], 1,
                        fluxes[j][f], 1);
            for (unsigned int k = 1; k < nDim; ++k)
            {
                Vmath::Vvtvp(nPts, m_D[vc[j][k]].GetValue(), 1, qfield[k][f], 1,
                             fluxes[j][f], 1, fluxes[j][f], 1);
            }
        }
    }
}

void SingleDiffusiveField::load_params()
{
    TokamakSystem::load_params();
    NekDouble k_c, lambda, T_bg;

    m_session->LoadParameter("k_c", k_c);
    m_session->LoadParameter("lambda", lambda);
    m_session->LoadParameter("T_bg", T_bg);

    k_par = 6.0 * k_c * (sqrt(2.0 * pow(M_PI, 3)) / lambda) * constants::epsilon_0 *
            constants::epsilon_0 * constants::c * pow(this->Tnorm * T_bg, 2.5) / sqrt(constants::m_e);
    // Correct for microns in epsilon_0 and density scale
    k_par *= 1e12 / this->Nnorm;
    // Convert to solver length and time scale
    k_par /= (this->omega_c * this->mesh_length * this->mesh_length);
    // multiply k_par by Z^-2 n^-1 in solver

    k_perp = lambda / (6.0 * sqrt(this->Tnorm * T_bg * pow(M_PI, 3.0))) *
             (sqrt(constants::m_p) / constants::c) * pow((constants::e / constants::epsilon_0_si), 2);
    k_perp *= this->Nnorm;
    k_perp /= (this->omega_c * this->mesh_length * this->mesh_length);
    // multiply k_perp by A^0.5 Z^2 n B^-1 in solver
}

void SingleDiffusiveField::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    TokamakSystem::v_ExtraFldOutput(fieldcoeffs, variables);
    const int nPhys   = m_fields[0]->GetNpoints();
    const int nCoeffs = m_fields[0]->GetNcoeffs();

    if (this->particles_enabled)
    {
        int i = 0;
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