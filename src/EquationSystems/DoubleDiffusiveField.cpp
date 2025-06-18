#include "DoubleDiffusiveField.hpp"

namespace NESO::Solvers::tokamak
{

/// Name of class
static std::string class_name;
std::string DoubleDiffusiveField::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "DoubleDiffusiveField", DoubleDiffusiveField::create,
        "Solves for two diffusive fields (n, e) with anisotropy");
/**
 * @brief Creates an instance of this class.
 */
static SU::EquationSystemSharedPtr create(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
{
    SU::EquationSystemSharedPtr p =
        MemoryManager<DoubleDiffusiveField>::AllocateSharedPtr(session, graph);
    p->InitObject();
    return p;
}

DoubleDiffusiveField::DoubleDiffusiveField(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
    : TokamakSystem(session, graph)
{
    this->n_indep_fields       = 1;
    this->n_fields_per_species = 2;
    this->required_fld_names   = {"n"};
}

void DoubleDiffusiveField::v_InitObject(bool DeclareFields)
{
    TokamakSystem::v_InitObject(DeclareFields);

    this->ne = std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]);
    this->Te = MemoryManager<MR::DisContField>::AllocateSharedPtr(*ne);

    m_session->MatchSolverInfo("SpectralVanishingViscosity", "True",
                               m_useSpecVanVisc, false);
    m_session->LoadParameter("epsilon", m_epsilon, 1.0);
    m_D     = std::vector<StdRegions::VarCoeffMap>(m_indfields.size());
    m_kpar  = std::vector<Array<OneD, NekDouble>>(m_indfields.size());
    m_kperp = std::vector<Array<OneD, NekDouble>>(m_indfields.size());

    // Setup diffusion object

    if (m_useSpecVanVisc)
    {
        m_session->LoadParameter("SVVCutoffRatio", m_sVVCutoffRatio, 0.75);
        m_session->LoadParameter("SVVDiffCoeff", m_sVVDiffCoeff, 0.1);
        m_factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
        m_factors[StdRegions::eFactorSVVDiffCoeff] = m_sVVDiffCoeff / m_epsilon;
    }
    CalcDiffTensor();

    m_ode.DefineOdeRhs(&DoubleDiffusiveField::DoOdeRhs, this);

    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            std::string diffName;
            m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
            m_diffusion = SolverUtils::GetDiffusionFactory().CreateInstance(
                diffName, diffName);
            m_diffusion->SetFluxVector(&DoubleDiffusiveField::GetFluxVectorDiff,
                                       this);
            m_diffusion->InitObject(m_session, m_indfields);
            break;
        }
        case MultiRegions::eGalerkin:
        {
            // Enable implicit solver for CG.
            m_ode.DefineImplicitSolve(&DoubleDiffusiveField::ImplicitTimeIntCG,
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
            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_DENSITY"));
            src_components.push_back(0);
            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_ENERGY"));
            src_components.push_back(0);
        }

        this->particle_sys->setup_evaluate_fields(this->E, this->B, this->ne,
                                                  this->Te, this->ve);

        this->particle_sys->finish_setup(this->src_fields, src_syms,
                                         src_components);
    }
}

void DoubleDiffusiveField::ImplicitTimeIntCG(
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
    CalcDiffTensor();
    Vmath::Zero(npoints, m_fields[0]->UpdatePhys(), 1);

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

        if (i % 2 == 1)
        {
            Vmath::Vadd(npoints, inarray[i - 1], 1, m_fields[0]->UpdatePhys(),
                        1, m_fields[0]->UpdatePhys(), 1);
            Vmath::Vdiv(npoints, outarray[i], 1, inarray[i - 1], 1, outarray[i],
                        1);
        }
        else if (i == nvariables - 1)
        {
            Vmath::Vdiv(npoints, outarray[i], 1, m_fields[0]->GetPhys(), 1,
                        outarray[i], 1);
        }
        m_indfields[i]->HelmSolve(outarray[i], m_indfields[i]->UpdateCoeffs(),
                                  m_factors, m_D[i]);

        m_indfields[i]->BwdTrans(m_indfields[i]->GetCoeffs(), outarray[i]);
    }

    Vmath::Zero(npoints, m_fields[0]->UpdatePhys(), 1);
    for (int i = 0; i < nvariables; ++i)
    {
        if (i % 2 == 1)
        {
            Vmath::Vadd(npoints, outarray[i - 1], 1, m_fields[0]->UpdatePhys(),
                        1, m_fields[0]->UpdatePhys(), 1);
            Vmath::Vmul(npoints, outarray[i], 1, outarray[i - 1], 1,
                        outarray[i], 1);
        }
        else if (i == nvariables - 1)
        {
            Vmath::Vmul(npoints, outarray[i], 1, m_fields[0]->GetPhys(), 1,
                        outarray[i], 1);
        }
        m_indfields[i]->SetPhysState(false);
    }
}

void DoubleDiffusiveField::CalcKPar(int f)
{
    // Change to fn of fields
    int npoints     = m_fields[0]->GetNpoints();
    NekDouble k_par = this->k_par;
    m_kpar[f]       = Array<OneD, NekDouble>(npoints, k_par);
}

void DoubleDiffusiveField::CalcKPerp(int f)
{
    // Change to fn of fields
    int npoints      = m_fields[0]->GetNpoints();
    NekDouble k_perp = this->k_perp;
    m_session->LoadParameter("k_perp", k_perp, 1.0);
    m_kperp[f] = Array<OneD, NekDouble>(npoints, k_perp);
}

void DoubleDiffusiveField::CalcDiffTensor()
{
    int npoints    = m_fields[0]->GetNpoints();
    int nvariables = m_indfields.size();

    for (int f = 0; f < nvariables; ++f)
    {
        CalcKPar(f);
        CalcKPerp(f);
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                Array<OneD, NekDouble> d(npoints, 0.0);

                for (int k = 0; k < npoints; k++)
                {
                    d[k] = (m_kpar[f][k] - m_kperp[f][k]) * b_unit[i][k] *
                           b_unit[j][k];

                    if (i == j)
                    {
                        d[k] += m_kperp[f][k];
                    }
                }

                m_D[f][vc[i][j]] = d;
            }
        }
    }
}

/**
 * @brief Populate rhs array ( @p out_arr ) with diffused quantitites
 * @param in_arr physical values of all fields
 * @param[out] out_arr output array (RHSs of time integration equations)
 */
void DoubleDiffusiveField::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time)
{
    if (m_explicitDiffusion || m_projectionType == MR::eDiscontinuous)
    {
        int nPts          = GetNpoints();
        size_t nvariables = in_arr.size();

        Array<OneD, Array<OneD, NekDouble>> tmp(nvariables);
        tmp = in_arr;
        CalcDiffTensor();
        m_diffusion->Diffuse(nvariables, m_indfields, tmp, out_arr);
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
        for (int i = 0; i < this->src_fields.size(); ++i)
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
void DoubleDiffusiveField::GetFluxVectorDiff(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{

    unsigned int nDim = qfield.size();
    unsigned int nFld = qfield[0].size();
    unsigned int nPts = qfield[0][0].size();
    for (int f = 0; f < nFld; ++f)
    {
        for (unsigned int j = 0; j < nDim; ++j)
        {
            // Calc diffusion of n with D tensor and n field
            Vmath::Vmul(nPts, m_D[f][vc[j][0]].GetValue(), 1, qfield[0][f], 1,
                        fluxes[j][f], 1);
            for (unsigned int k = 1; k < nDim; ++k)
            {
                Vmath::Vvtvp(nPts, m_D[f][vc[j][k]].GetValue(), 1, qfield[k][f],
                             1, fluxes[j][f], 1, fluxes[j][f], 1);
            }
        }
    }
}

void DoubleDiffusiveField::load_params()
{
    TokamakSystem::load_params();
    // Adiabatic gamma

    m_session->LoadParameter("k_B", this->m_k_B);
    m_session->LoadParameter("k_par", this->k_par, 100.0);
    m_session->LoadParameter("k_perp", this->k_perp, 1.0);
}

void DoubleDiffusiveField::v_ExtraFldOutput(
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
            Array<OneD, NekDouble> SrcFwd1(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[i]->GetPhys(),
                                           SrcFwd1);
            fieldcoeffs.push_back(SrcFwd1);

            variables.push_back(v.name + "_SOURCE_ENERGY");
            Array<OneD, NekDouble> SrcFwd2(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[i + 1]->GetPhys(),
                                           SrcFwd2);
            fieldcoeffs.push_back(SrcFwd2);
            i += 2;
        }
    }
}
} // namespace NESO::Solvers::tokamak