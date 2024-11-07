#include "SingleDiffusiveField.hpp"

namespace NESO::Solvers::tokamak
{

/// Name of class
static std::string class_name;
std::string SingleDiffusiveField::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "Single Diffusive Field", SingleDiffusiveField::create,
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
        this->required_fld_names.push_back("n_src");
    }
}

void SingleDiffusiveField::v_InitObject(bool DeclareFields)
{
    TokamakSystem::v_InitObject(DeclareFields);
    // Setup diffusion object
    std::string diffName;
    m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
    m_diffusion =
        SolverUtils::GetDiffusionFactory().CreateInstance(diffName, diffName);
    m_diffusion->SetFluxVector(&SingleDiffusiveField::GetFluxVectorDiff, this);
    m_diffusion->InitObject(m_session, m_fields);
    m_ode.DefineOdeRhs(&SingleDiffusiveField::DoOdeRhs, this);
}

void SingleDiffusiveField::CalcKPar()
{
    // Change to fn of fields
    int npoints = m_fields[0]->GetNpoints();
    NekDouble k_par;
    m_session->LoadParameter("k_par", k_par, 100.0);
    m_kpar = Array<OneD, NekDouble>(npoints, k_par);
}

void SingleDiffusiveField::CalcKPerp()
{
    // Change to fn of fields
    int npoints = m_fields[0]->GetNpoints();
    NekDouble k_perp;
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
    int npts          = GetNpoints();
    size_t nvariables = in_arr.size();

    Array<OneD, MultiRegions::ExpListSharedPtr> diff_fields(nvariables);
    Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>(out_arr[i].size(), 0.0);
        diff_fields[i]  = m_fields[i];
    }
    CalcDiffTensor();
    m_diffusion->Diffuse(nvariables, diff_fields, in_arr, outarrayDiff);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(out_arr[i].size(), outarrayDiff[i], 1, out_arr[i], 1,
                    out_arr[i], 1);
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
} // namespace NESO::Solvers::tokamak