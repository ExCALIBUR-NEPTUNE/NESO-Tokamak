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
    m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
        m_session, m_spacedim, m_graph);

    this->ne = std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]);
    this->Te = MemoryManager<MR::DisContField>::AllocateSharedPtr(*ne);

    pe_idx = 2 * this->n_species;
    int s  = 0;
    for (const auto &[k, v] : this->neso_config->get_species())
    {
        ni_idx.push_back(2 * s);
        pi_idx.push_back(1 + 2 * s);

        double charge, mass;
        this->neso_config->load_species_parameter(k, "Charge", charge);
        this->neso_config->load_species_parameter(k, "Mass", mass);
        m_varConv->charge[s] = charge;
        m_varConv->mass[s]   = mass;
        s++;
    }

    m_session->MatchSolverInfo("SpectralVanishingViscosity", "True",
                               m_useSpecVanVisc, false);
    m_session->LoadParameter("epsilon", m_epsilon, 1.0);

    // Setup diffusion object

    if (m_useSpecVanVisc)
    {
        m_session->LoadParameter("SVVCutoffRatio", m_sVVCutoffRatio, 0.75);
        m_session->LoadParameter("SVVDiffCoeff", m_sVVDiffCoeff, 0.1);
        m_factors[StdRegions::eFactorSVVCutoffRatio] = m_sVVCutoffRatio;
        m_factors[StdRegions::eFactorSVVDiffCoeff] = m_sVVDiffCoeff / m_epsilon;
    }

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
        this->src_fields.emplace_back(
            MemoryManager<MR::DisContField>::AllocateSharedPtr(
                *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
        src_syms.push_back(Sym<REAL>("ELECTRON_SOURCE_DENSITY"));
        src_components.push_back(0);

        this->src_fields.emplace_back(
            MemoryManager<MR::DisContField>::AllocateSharedPtr(
                *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
        src_syms.push_back(Sym<REAL>("ELECTRON_SOURCE_ENERGYY"));
        src_components.push_back(0);

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

        this->src_fields.emplace_back(
            MemoryManager<MR::DisContField>::AllocateSharedPtr(
                *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
        src_syms.push_back(Sym<REAL>("ELECTRON_SOURCE_ENERGY"));
        src_components.push_back(0);

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

        if (this->particles_enabled)
        {
            for (int i = 0; i < this->src_fields.size(); ++i)
            {
                Vmath::Vadd(outarray[i].size(), outarray[i], 1,
                            this->src_fields[i]->GetPhys(), 1, outarray[i], 1);
            }
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
                                  m_factors, m_D);

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

void DoubleDiffusiveField::CalcKPar(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr, int f)
{
    int npoints = m_fields[0]->GetNpoints();
    double Z;
    this->neso_config->load_species_parameter(f, "Charge", Z);

    for (int p = 0; p < npoints; ++p)
    {
        m_kpar[p] = this->k_ci * this->k_par * pow(in_arr[pi_idx[f]][p], 2.5) /
                    (Z * Z * in_arr[ni_idx[f]][p]);
    }
}

void DoubleDiffusiveField::CalcKPerp(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr, int f)
{
    int npoints = m_fields[0]->GetNpoints();
    double Z, A;
    this->neso_config->load_species_parameter(f, "Charge", Z);
    this->neso_config->load_species_parameter(f, "Mass", A);

    for (int p = 0; p < npoints; ++p)
    {
        m_kperp[p] = this->k_perp * Z * Z * std::sqrt(A) *
                     in_arr[ni_idx[f]][p] /
                     (sqrt(in_arr[pi_idx[f]][p]) * this->mag_B[p]);
    }
}

void DoubleDiffusiveField::CalcKappaPar(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr, int f)
{
    int npoints = m_fields[0]->GetNpoints();
    double Z, A;
    this->neso_config->load_species_parameter(f, "Charge", Z);
    this->neso_config->load_species_parameter(f, "Mass", A);
    for (int p = 0; p < npoints; ++p)
    {
        m_kpar[p] = this->k_ci * this->k_par * sqrt(this->me) *
                    pow(in_arr[pi_idx[f]][p], 2.5) / (sqrt(A) * Z * Z);
    }
}

void DoubleDiffusiveField::CalcKappaPerp(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr, int f)
{
    int npoints = m_fields[0]->GetNpoints();
    double Z, A;
    this->neso_config->load_species_parameter(f, "Charge", Z);
    this->neso_config->load_species_parameter(f, "Mass", A);
    for (int p = 0; p < npoints; ++p)
    {
        m_kperp[p] = this->k_perp * Z * Z * std::sqrt(A) *
                     in_arr[ni_idx[f]][p] /
                     (sqrt(in_arr[pi_idx[f]][p]) * this->mag_B[p]);
    }
}
void DoubleDiffusiveField::CalcKappaPar(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr)
{
    int npoints = m_fields[0]->GetNpoints();
    m_kpar      = Array<OneD, NekDouble>(npoints, k_par);
    m_kpar      = Array<OneD, NekDouble>(npoints, this->k_ce * this->k_par);
}

void DoubleDiffusiveField::CalcKappaPerp(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr)
{
    // Change to fn of fields
    int npoints      = m_fields[0]->GetNpoints();
    NekDouble k_perp = this->k_perp;
    m_kperp          = Array<OneD, NekDouble>(npoints, k_perp);
}

void DoubleDiffusiveField::CalcDiffTensor()
{
    int npoints    = m_fields[0]->GetNpoints();
    int nvariables = m_indfields.size();

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
void DoubleDiffusiveField::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time)
{
    m_varConv->GetElectronDensity(in_arr, m_fields[0]->UpdatePhys());

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

void DoubleDiffusiveField::DoDiffusion(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{

    size_t nvariables = inarray.size();
    size_t npointsIn  = GetNpoints();
    size_t npointsOut =
        npointsIn; // If ALE then outarray is in coefficient space
    size_t nTracePts = GetTraceTotPoints();

    // this should be preallocated
    Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvariables);
    for (size_t i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>(npointsOut, 0.0);
    }

    // Get primitive variables [n,T]
    Array<OneD, Array<OneD, NekDouble>> inarrayDiff(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inFwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inBwd(nvariables);

    for (size_t i = 0; i < nvariables; ++i)
    {
        inarrayDiff[i] = Array<OneD, NekDouble>{npointsIn};
        inFwd[i]       = Array<OneD, NekDouble>{nTracePts};
        inBwd[i]       = Array<OneD, NekDouble>{nTracePts};
    }

    // Extract temperature
    m_varConv->GetElectronTemperature(inarray, inarrayDiff[pe_idx]);

    int s = 0;
    double mass;

    for (const auto &[k, v] : this->neso_config->get_species())
    {
        Vmath::Vcopy(npointsIn, inarray[ni_idx[s]], 1, inarrayDiff[ni_idx[s]],
                     1);
        this->neso_config->load_species_parameter(k, "Mass", mass);
        m_varConv->GetIonTemperature(s, mass, inarray, inarrayDiff[pi_idx[s]]);

        s++;
    }

    // Repeat calculation for trace space
    if (pFwd == NullNekDoubleArrayOfArray || pBwd == NullNekDoubleArrayOfArray)
    {
        inFwd = NullNekDoubleArrayOfArray;
        inBwd = NullNekDoubleArrayOfArray;
    }
    else
    {
        m_varConv->GetElectronTemperature(pFwd, inFwd[pe_idx]);
        m_varConv->GetElectronTemperature(pBwd, inBwd[pe_idx]);
        s = 0;
        for (const auto &[k, v] : this->neso_config->get_species())
        {
            this->neso_config->load_species_parameter(k, "Mass", mass);
            Vmath::Vcopy(npointsIn, pFwd[ni_idx[s]], 1, inFwd[ni_idx[s]], 1);
            Vmath::Vcopy(npointsIn, pBwd[ni_idx[s]], 1, inBwd[ni_idx[s]], 1);

            m_varConv->GetIonTemperature(s, mass, pFwd, inFwd[pi_idx[s]]);
            m_varConv->GetIonTemperature(s, mass, pBwd, inBwd[pi_idx[s]]);

            s++;
        }
    }

    m_diffusion->Diffuse(nvariables, m_fields, inarrayDiff, outarrayDiff, inFwd,
                         inBwd);

    for (size_t i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(npointsOut, outarrayDiff[i], 1, outarray[i], 1, outarray[i],
                    1);
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
    int s             = 0;
    for (const auto &[k, v] : this->neso_config->get_species())
    {
        CalcKPar(in_arr, k);
        CalcKPerp(in_arr, k);
        CalcDiffTensor();

        for (unsigned int j = 0; j < nDim; ++j)
        {
            // Calc diffusion of n with D tensor and n field
            Vmath::Vmul(nPts, m_D[vc[j][0]].GetValue(), 1, qfield[0][ni_idx[s]],
                        1, fluxes[j][ni_idx[s]], 1);
            for (unsigned int k = 1; k < nDim; ++k)
            {
                Vmath::Vvtvp(nPts, m_D[vc[j][k]].GetValue(), 1,
                             qfield[k][ni_idx[s]], 1, fluxes[j][ni_idx[s]], 1,
                             fluxes[j][ni_idx[s]], 1);
            }
        }

        CalcKappaPar(in_arr, k);
        CalcKappaPerp(in_arr, k);
        CalcDiffTensor();
        for (unsigned int j = 0; j < nDim; ++j)
        {
            // Calc diffusion of n with D tensor and n field
            Vmath::Vmul(nPts, m_D[vc[j][0]].GetValue(), 1, qfield[0][pi_idx[s]],
                        1, fluxes[j][pi_idx[s]], 1);
            for (unsigned int k = 1; k < nDim; ++k)
            {
                Vmath::Vvtvp(nPts, m_D[vc[j][k]].GetValue(), 1,
                             qfield[k][pi_idx[s]], 1, fluxes[j][pi_idx[s]], 1,
                             fluxes[j][pi_idx[s]], 1);
            }
        }
    }

    CalcKappaPar(in_arr);
    CalcKappaPerp(in_arr);
    CalcDiffTensor();
    for (unsigned int j = 0; j < nDim; ++j)
    {
        Vmath::Vmul(nPts, m_D[vc[j][0]].GetValue(), 1, qfield[0][pe_idx], 1,
                    fluxes[j][pe_idx], 1);
        for (unsigned int k = 1; k < nDim; ++k)
        {
            Vmath::Vvtvp(nPts, m_D[vc[j][k]].GetValue(), 1, qfield[k][pe_idx],
                         1, fluxes[j][pe_idx], 1, fluxes[j][pe_idx], 1);
        }
    }
}

void DoubleDiffusiveField::load_params()
{
    TokamakSystem::load_params();
    NekDouble lambda, T_bg;

    m_session->LoadParameter("k_ci", this->k_ci, 3.9);
    m_session->LoadParameter("k_ce", this->k_ce, 3.16);
    m_session->LoadParameter("lambda", lambda);
    m_session->LoadParameter("T_bg", T_bg);

    k_par = 6.0 * (sqrt(2.0 * pow(M_PI, 3)) / lambda) * constants::epsilon_0 * constants::epsilon_0 *
            constants::c * pow(this->Tnorm, 2.5) / sqrt(constants::m_e);
    // Correct for microns in epsilon_0 and density scale
    k_par *= 1e12 / this->Nnorm;
    // Convert to solver length and time scale
    k_par /= (this->omega_c * this->mesh_length * this->mesh_length);
    // multiply k_par by Z^-2 n^-1 T^2.5 in solver

    k_perp = lambda / (6.0 * sqrt(this->Tnorm * pow(M_PI, 3.0))) *
             (sqrt(constants::m_p) / constants::c) * pow((constants::e / constants::epsilon_0_si), 2);
    k_perp *= this->Nnorm;
    k_perp /= (this->omega_c * this->mesh_length * this->mesh_length);
    // multiply k_perp by A^0.5 Z^2 n B^-1 T^-0.5 in solver
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