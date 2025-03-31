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
    this->required_fld_names = {"n", "p", "T"};

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
        for (auto &[k, v] : this->particle_sys->get_species())
        {
            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            this->energy_src_fields.emplace_back(
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
    Vmath::Vdiv(nPts, tmp[p_idx], 1, tmp[n_idx], 1, tmp[T_idx], 1);

    CalcDiffTensor();
    m_diffusion->Diffuse(nvariables, m_fields, tmp, out_arr);

    // Add forcing terms
    for (auto &x : m_forcing)
    {
        x->Apply(m_fields, tmp, out_arr, time);
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

    int n_idx = this->field_to_index["n"];
    int p_idx = this->field_to_index["p"];
    int T_idx = this->field_to_index["T"];

    for (unsigned int j = 0; j < nDim; ++j)
    {
        // Calc diffusion of n with D tensor and n field
        Vmath::Vmul(nPts, m_D[vc[j][0]].GetValue(), 1, qfield[0][n_idx], 1,
                    fluxes[j][n_idx], 1);
        // Calc diffusion of p with kappa tensor and T field
        Vmath::Vmul(nPts, m_kappa[vc[j][0]].GetValue(), 1, qfield[0][T_idx], 1,
                    fluxes[j][p_idx], 1);
        for (unsigned int k = 1; k < nDim; ++k)
        {
            Vmath::Vvtvp(nPts, m_D[vc[j][k]].GetValue(), 1, qfield[k][n_idx], 1,
                         fluxes[j][n_idx], 1, fluxes[j][n_idx], 1);
            Vmath::Vvtvp(nPts, m_kappa[vc[j][k]].GetValue(), 1,
                         qfield[k][T_idx], 1, fluxes[j][p_idx], 1,
                         fluxes[j][p_idx], 1);
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

bool DoubleDiffusiveField::v_PostIntegrate(int step)
{
    unsigned int nPts = GetNpoints();
    int n_idx         = this->field_to_index["n"];
    int p_idx         = this->field_to_index["p"];
    int T_idx         = this->field_to_index["T"];
    Vmath::Vdiv(nPts, m_fields[p_idx]->GetPhys(), 1, m_fields[n_idx]->GetPhys(),
                1, m_fields[T_idx]->UpdatePhys(), 1);
    m_fields[T_idx]->FwdTrans(m_fields[T_idx]->GetPhys(),
                              m_fields[T_idx]->UpdateCoeffs());
    return TokamakSystem::v_PostIntegrate(step);
}

void DoubleDiffusiveField::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    TokamakSystem::v_ExtraFldOutput(fieldcoeffs, variables);
    const int nPhys   = m_fields[0]->GetNpoints();
    const int nCoeffs = m_fields[0]->GetNcoeffs();
    for (auto s : this->particle_sys->get_species())
    {
        variables.push_back(s.second.name + "_SOURCE_DENSITY");
        Array<OneD, NekDouble> SrcFwd1(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(
            this->density_src_fields[s.first]->GetPhys(), SrcFwd1);
        fieldcoeffs.push_back(SrcFwd1);

        variables.push_back(s.second.name + "_SOURCE_ENERGY");
        Array<OneD, NekDouble> SrcFwd2(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(
            this->energy_src_fields[s.first]->GetPhys(), SrcFwd2);
        fieldcoeffs.push_back(SrcFwd2);
    }
}
} // namespace NESO::Solvers::tokamak