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
        this->required_fld_names = {"n_src", "E_src"};
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
}

void DoubleDiffusiveField::CalcKPar()
{
    // Change to fn of fields
    int npoints = m_fields[0]->GetNpoints();
    NekDouble k_par;
    m_session->LoadParameter("k_par", k_par, 100.0);
    m_kpar = Array<OneD, NekDouble>(npoints, k_par);
}

void DoubleDiffusiveField::CalcKPerp()
{
    // Change to fn of fields
    int npoints = m_fields[0]->GetNpoints();
    NekDouble k_perp;
    m_session->LoadParameter("k_perp", k_perp, 1.0);
    m_kperp = Array<OneD, NekDouble>(npoints, k_perp);
}

void DoubleDiffusiveField::CalcDiffTensor()
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

void DoubleDiffusiveField::CalcKappaPar()
{
    // Change to fn of T
    int npoints = m_fields[0]->GetNpoints();
    NekDouble kappa_par;
    m_session->LoadParameter("kappa_par", kappa_par, 100.0);
    m_kappapar = Array<OneD, NekDouble>(npoints, 2 * kappa_par / 3);
}

void DoubleDiffusiveField::CalcKappaPerp()
{
    // Change to fn of T
    int npoints = m_fields[0]->GetNpoints();
    NekDouble kappa_perp;
    m_session->LoadParameter("kappa_perp", kappa_perp, 1.0);
    m_kappaperp = Array<OneD, NekDouble>(npoints, 2 * kappa_perp / 3);
}

void DoubleDiffusiveField::CalcKappaTensor()
{
    int npoints = m_fields[0]->GetNpoints();
    CalcKappaPar();
    CalcKappaPerp();
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Array<OneD, NekDouble> kappa(npoints, 0.0);
            for (int k = 0; k < npoints; k++)
            {
                kappa[k] = (m_kappapar[k] - m_kappaperp[k]) * b_unit[i][k] *
                           b_unit[j][k];
                if (i == j)
                {
                    kappa[k] += m_kappaperp[k];
                }
            }
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
    int npts          = GetNpoints();
    size_t nvariables = in_arr.size();

    Array<OneD, MultiRegions::ExpListSharedPtr> diff_fields(nvariables);
    Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvariables);

    CalcDiffTensor();
    CalcKappaTensor();
    m_diffusion->Diffuse(nvariables, m_fields, in_arr, out_arr);

    // Add forcing terms
    for (auto &x : m_forcing)
    {
        x->Apply(m_fields, in_arr, out_arr, time);
    }
    // Recalculate T from p and n
    int T_idx = this->field_to_index["T"];
    int n_idx = this->field_to_index["n"];
    int p_idx = this->field_to_index["p"];
    Vmath::Vdiv(npts, out_arr[p_idx], 1, out_arr[n_idx], 1, out_arr[T_idx], 1);
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

    int T_idx = this->field_to_index["T"];
    int n_idx = this->field_to_index["n"];
    int p_idx = this->field_to_index["p"];
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
}

} // namespace NESO::Solvers::tokamak