#include "DoubleDiffusiveField.hpp"

namespace NESO::Solvers::tokamak
{

/// Name of class
static std::string class_name;
std::string DoubleDiffusiveField::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "DoubleDiffusiveField", DoubleDiffusiveField::create,
        "Solves for two diffusive fields (n, p) with anisotropy");
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
    this->required_fld_names = {"n", "p"};

    if (this->particles_enabled)
    {
        this->required_fld_names.push_back("E_src");
        for (auto &[k, v] : particle_sys->get_species())
        {
            this->required_fld_names.push_back(v.name + "_src");
        }
    }
}

void DoubleDiffusiveField::v_InitObject(bool DeclareFields)
{
    TokamakSystem::v_InitObject(DeclareFields);

    m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
        m_session, m_spacedim, m_graph);

    // Setup diffusion object
    std::string diffName;
    m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
    m_diffusion =
        SolverUtils::GetDiffusionFactory().CreateInstance(diffName, diffName);
    m_diffusion->SetFluxVector(&DoubleDiffusiveField::GetFluxVectorDiff, this);
    m_diffusion->InitObject(m_session, m_fields);
    m_ode.DefineOdeRhs(&DoubleDiffusiveField::DoOdeRhs, this);

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

            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_DENSITY"));
            src_components.push_back(0);

            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_DENSITY"));
            src_components.push_back(0);
            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_ENERGY"));
            src_components.push_back(0);
        }
        std::vector<MR::DisContFieldSharedPtr> src_fields =
            this->density_src_fields;
        src_fields.insert(src_fields.end(), this->energy_src_fields.begin(),
                          this->energy_src_fields.end());
        this->particle_sys->finish_setup(src_fields, src_syms, src_components);
    }
}

void DoubleDiffusiveField::CalcKPar()
{
    // Change to fn of fields
    int npoints     = m_fields[0]->GetNpoints();
    NekDouble k_par = this->k_par;
    m_kpar          = Array<OneD, NekDouble>(npoints, k_par);
}

void DoubleDiffusiveField::CalcKPerp()
{
    // Change to fn of fields
    int npoints      = m_fields[0]->GetNpoints();
    NekDouble k_perp = this->k_perp;
    m_kperp          = Array<OneD, NekDouble>(npoints, k_perp);
}

void DoubleDiffusiveField::CalcKappaPar()
{
    // Change to fn of T
    int npoints         = m_fields[0]->GetNpoints();
    NekDouble kappa_par = this->kappa_par;
    m_kappapar          = Array<OneD, NekDouble>(npoints, 2 * kappa_par / 3);
}

void DoubleDiffusiveField::CalcKappaPerp()
{
    // Change to fn of T
    int npoints          = m_fields[0]->GetNpoints();
    NekDouble kappa_perp = this->kappa_perp;

    m_kappaperp = Array<OneD, NekDouble>(npoints, 2 * kappa_perp / 3);
}

void DoubleDiffusiveField::CalcDiffTensor()
{
    int npoints = m_fields[0]->GetNpoints();
    CalcKPar();
    CalcKPerp();
    CalcKappaPar();
    CalcKappaPerp();
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Array<OneD, NekDouble> d(npoints, 0.0);
            Array<OneD, NekDouble> kappa(npoints, 0.0);

            for (int k = 0; k < npoints; k++)
            {
                d[k] = (m_kpar[k] - m_kperp[k]) * b_unit[i][k] * b_unit[j][k];
                kappa[k] = (m_kappapar[k] - m_kappaperp[k]) * b_unit[i][k] *
                           b_unit[j][k];
                if (i == j)
                {
                    d[k] += m_kperp[k];
                    kappa[k] += m_kappaperp[k];
                }
            }
            m_D[vc[i][j]]     = d;
            m_kappa[vc[i][j]] = kappa;
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
    int nPts          = GetNpoints();
    size_t nvariables = in_arr.size();
    size_t nTracePts  = GetTraceTotPoints();
    int n_idx         = this->field_to_index["n"];
    int p_idx         = this->field_to_index["p"];
    int T_idx         = this->field_to_index["T"];

    Array<OneD, Array<OneD, NekDouble>> tmp(nvariables);
    tmp = in_arr;

    // Store forwards/backwards space along trace space
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvariables);

    for (size_t i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        m_fields[i]->GetFwdBwdTracePhys(tmp[i], Fwd[i], Bwd[i]);
    }

    DoDiffusion(tmp, out_arr, Fwd, Bwd);

    if (this->particles_enabled)
    {
        for (int i = 0; i < this->particle_sys->get_species().size(); ++i)
        {
            Vmath::Vadd(out_arr[0].size(), out_arr[0], 1,
                        this->src_fields[2 * i]->GetPhys(), 1, out_arr[0], 1);
            Vmath::Vadd(out_arr[1].size(), out_arr[1], 1,
                        this->src_fields[2 * i + 1]->GetPhys(), 1, out_arr[1],
                        1);
        }
    }

    // Add forcing terms
    for (auto &x : m_forcing)
    {
        x->Apply(m_fields, tmp, out_arr, time);
    }
}

void DoubleDiffusiveField::DoDiffusion(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    CalcDiffTensor();

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
    Vmath::Vcopy(npointsIn, inarray[0], 1, inarrayDiff[0], 1);

    // Extract temperature
    m_varConv->GetElectronTemperature(inarray, inarrayDiff[1]);

    // Repeat calculation for trace space
    if (pFwd == NullNekDoubleArrayOfArray || pBwd == NullNekDoubleArrayOfArray)
    {
        inFwd = NullNekDoubleArrayOfArray;
        inBwd = NullNekDoubleArrayOfArray;
    }
    else
    {
        Vmath::Vcopy(npointsIn, pFwd[0], 1, inFwd[0], 1);
        Vmath::Vcopy(npointsIn, pBwd[0], 1, inBwd[0], 1);

        m_varConv->GetElectronTemperature(pFwd, inFwd[1]);
        m_varConv->GetElectronTemperature(pBwd, inBwd[1]);
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
    unsigned int nPts = qfield[0][0].size();

    for (unsigned int j = 0; j < nDim; ++j)
    {
        // Calc n flux
        Vmath::Vmul(nPts, m_D[vc[j][0]].GetValue(), 1, qfield[0][0], 1,
                    fluxes[j][0], 1);
        // Calc energy flux
        Vmath::Vmul(nPts, m_kappa[vc[j][0]].GetValue(), 1, qfield[0][1], 1,
                    fluxes[j][1], 1);
        for (unsigned int k = 1; k < nDim; ++k)
        {
            Vmath::Vvtvp(nPts, m_D[vc[j][k]].GetValue(), 1, qfield[k][0], 1,
                         fluxes[j][0], 1, fluxes[j][0], 1);
            Vmath::Vvtvp(nPts, m_kappa[vc[j][k]].GetValue(), 1, qfield[k][1], 1,
                         fluxes[j][1], 1, fluxes[j][1], 1);
        }
    }
}

void DoubleDiffusiveField::load_params()
{
    TokamakSystem::load_params();
    // Adiabatic gamma
    m_session->LoadParameter("gamma", this->m_gamma);
    m_session->LoadParameter("k_B", this->m_k_B);
    m_session->LoadParameter("k_par", this->k_par, 100.0);
    m_session->LoadParameter("k_perp", this->k_perp, 1.0);
    m_session->LoadParameter("kappa_par", this->kappa_par, 100.0);
    m_session->LoadParameter("kappa_perp", this->kappa_perp, 1.0);
}

void DoubleDiffusiveField::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    TokamakSystem::v_ExtraFldOutput(fieldcoeffs, variables);

    bool extraFields;
    m_session->MatchSolverInfo("OutputExtraFields", "True", extraFields, true);
    if (extraFields)
    {
        const int nPhys   = m_fields[0]->GetNpoints();
        const int nCoeffs = m_fields[0]->GetNcoeffs();

        Array<OneD, Array<OneD, NekDouble>> tmp(m_fields.size());

        for (int i = 0; i < m_fields.size(); ++i)
        {
            tmp[i] = m_fields[i]->GetPhys();
        }

        Array<OneD, NekDouble> pressure(nPhys), temperature(nPhys);
        m_varConv->GetElectronPressure(tmp, pressure);
        m_varConv->GetElectronTemperature(tmp, temperature);

        Array<OneD, NekDouble> pFwd(nCoeffs), TFwd(nCoeffs);

        m_fields[0]->FwdTransLocalElmt(pressure, pFwd);
        m_fields[0]->FwdTransLocalElmt(temperature, TFwd);

        variables.push_back("pe");
        variables.push_back("Te");

        fieldcoeffs.push_back(pFwd);
        fieldcoeffs.push_back(TFwd);

        if (this->particles_enabled)
        {
            int i = 0;
            for (auto s : this->particle_sys->get_species())
            {
                variables.push_back(s.second.name + "_SOURCE_DENSITY");
                Array<OneD, NekDouble> SrcFwd1(nCoeffs);
                m_fields[0]->FwdTransLocalElmt(
                    this->src_fields[2 * i + 2]->GetPhys(), SrcFwd1);
                fieldcoeffs.push_back(SrcFwd1);

                variables.push_back(s.second.name + "_SOURCE_ENERGY");
                Array<OneD, NekDouble> SrcFwd2(nCoeffs);
                m_fields[0]->FwdTransLocalElmt(
                    this->src_fields[2 * i + 3]->GetPhys(), SrcFwd2);
                fieldcoeffs.push_back(SrcFwd2);
                ++i;
            }
        }
    }
}
} // namespace NESO::Solvers::tokamak