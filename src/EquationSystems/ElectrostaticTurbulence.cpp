#include "ElectrostaticTurbulence.hpp"

namespace NESO::Solvers::tokamak
{

/// Name of class
static std::string class_name;
std::string ElectrostaticTurbulence::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "ElectrostaticTurbulence", ElectrostaticTurbulence::create,
        "Solves electrostatic turbulence with anisotropic diffusion");
/**
 * @brief Creates an instance of this class.
 */
static SU::EquationSystemSharedPtr create(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
{
    SU::EquationSystemSharedPtr p =
        MemoryManager<ElectrostaticTurbulence>::AllocateSharedPtr(session,
                                                                  graph);
    p->InitObject();
    return p;
}

ElectrostaticTurbulence::ElectrostaticTurbulence(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
    : TokamakSystem(session, graph), ExB_vel(graph->GetSpaceDimension())
{
    this->required_fld_names = {"phi", "w", "Te", "ne", "Ti", "ni"};
    this->int_fld_names      = {"w", "Te", "ne", "Ti", "ni"};

    if (this->particles_enabled)
    {
        this->required_fld_names.push_back("n_src");
        this->required_fld_names.push_back("E_src");
    }
}

void ElectrostaticTurbulence::v_InitObject(bool DeclareFields)
{
    TokamakSystem::v_InitObject(DeclareFields);

    // Since we are starting from a setup where each field is defined to be a
    // discontinuous field (and thus support DG), the first thing we do is to
    // recreate the phi field so that it is continuous, in order to support the
    // Poisson solve. Note that you can still perform a Poisson solve using a
    // discontinuous field, which is done via the hybridisable discontinuous
    // Galerkin (HDG) approach.
    int phi_idx       = this->field_to_index["phi"];
    m_fields[phi_idx] = MemoryManager<MR::ContField>::AllocateSharedPtr(
        m_session, m_graph, m_session->GetVariable(phi_idx), true, true);

    std::string diffName;
    m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
    m_diffusion =
        SolverUtils::GetDiffusionFactory().CreateInstance(diffName, diffName);
    m_diffusion->SetFluxVector(&ElectrostaticTurbulence::GetFluxVectorDiff,
                               this);
    m_diffusion->InitObject(m_session, m_fields);

    // Create storage for ExB vel
    int npts = GetNpoints();
    for (int i = 0; i < m_graph->GetSpaceDimension(); ++i)
    {
        this->ExB_vel[i] = Array<OneD, NekDouble>(npts, 0.0);
    }

    // Define the normal velocity fields.
    // These are populated at each step (by reference) in calls to GetVnAdv()
    if (m_fields[0]->GetTrace())
    {
        auto nTrace         = GetTraceNpoints();
        this->norm_vel_elec = Array<OneD, NekDouble>(nTrace);
    }

    // Create Riemann solvers (one per advection object) and set normal velocity
    // callback functions
    this->riemann_solver = SU::GetRiemannSolverFactory().CreateInstance(
        this->riemann_solver_type, m_session);
    this->riemann_solver->SetScalar(
        "Vn", &ElectrostaticTurbulence::GetAdvVelNormElec, this);

    // Setup advection object
    m_advection = SU::GetAdvectionFactory().CreateInstance(this->adv_type,
                                                           this->adv_type);
    m_advection->SetFluxVector(&ElectrostaticTurbulence::GetFluxVectorElec,
                               this);
    m_advection->SetRiemannSolver(this->riemann_solver);
    m_advection->InitObject(m_session, m_fields);
    m_ode.DefineOdeRhs(&ElectrostaticTurbulence::DoOdeRhs, this);
}

void ElectrostaticTurbulence::CalcKPar()
{
    // Change to fn of fields
    int npoints = m_fields[0]->GetNpoints();
    NekDouble k_par;
    m_session->LoadParameter("k_par", k_par, 100.0);
    m_kpar = Array<OneD, NekDouble>(npoints, k_par);
}

void ElectrostaticTurbulence::CalcKPerp()
{
    // Change to fn of fields
    int npoints = m_fields[0]->GetNpoints();
    NekDouble k_perp;
    m_session->LoadParameter("k_perp", k_perp, 1.0);
    m_kperp = Array<OneD, NekDouble>(npoints, k_perp);
}

void ElectrostaticTurbulence::CalcDiffTensor()
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

void ElectrostaticTurbulence::CalcKappaPar()
{
    // Change to fn of T
    int npoints = m_fields[0]->GetNpoints();
    NekDouble kappa_par;
    m_session->LoadParameter("kappa_par", kappa_par, 100.0);
    m_kappapar = Array<OneD, NekDouble>(npoints, kappa_par);
}

void ElectrostaticTurbulence::CalcKappaPerp()
{
    // Change to fn of T
    int npoints = m_fields[0]->GetNpoints();
    NekDouble kappa_perp;
    m_session->LoadParameter("kappa_perp", kappa_perp, 1.0);
    m_kappaperp = Array<OneD, NekDouble>(npoints, kappa_perp);
}

void ElectrostaticTurbulence::CalcKappaTensor()
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
 * @brief Populate rhs array ( @p out_arr ) where electrostatic turbulence
 * equations (Hasegana-Wakatani) are used for the advection velocities.
 *
 * @param in_arr physical values of all fields
 * @param[out] out_arr output array (RHSs of time integration equations)
 */
void ElectrostaticTurbulence::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time)
{

    // Get field indices
    int npts    = GetNpoints();
    int phi_idx = this->field_to_index["phi"];
    int w_idx   = this->field_to_index["w"];
    int ne_idx  = this->field_to_index["ne"];

    if (this->m_explicitAdvection)
    {
        zero_array_of_arrays(out_arr);
    }

    // Helmholtz Solve for electrostatic potential
    SolvePhi(in_arr);
    // Calculate E
    ComputeE();

    // Calculate ExB velocity
    Vmath::Vvtvvtm(npts, E[1]->GetPhys(), 1, B[2]->GetPhys(), 1,
                   E[2]->GetPhys(), 1, B[1]->GetPhys(), 1, this->ExB_vel[0], 1);

    Vmath::Vvtvvtm(npts, E[2]->GetPhys(), 1, B[0]->GetPhys(), 1,
                   E[0]->GetPhys(), 1, B[2]->GetPhys(), 1, this->ExB_vel[1], 1);

    if (this->m_graph->GetMeshDimension() == 3)
    {
        Vmath::Vvtvvtm(npts, E[0]->GetPhys(), 1, B[1]->GetPhys(), 1,
                       E[1]->GetPhys(), 1, B[0]->GetPhys(), 1, this->ExB_vel[2],
                       1);
    }

    // Advect T, ne and w

    size_t nvariables = in_arr.size();

    m_advection->Advect(nvariables, m_fields, this->ExB_vel, in_arr, out_arr,
                        time);
    for (auto i = 0; i < nvariables; ++i)
    {
        Vmath::Neg(npts, out_arr[i], 1);
    }

    // HW equation
    //  Add \alpha*(\phi-ne) to RHS
    Array<OneD, NekDouble> HWterm_2D_alpha(npts);
    Vmath::Vsub(npts, m_fields[phi_idx]->GetPhys(), 1,
                m_fields[ne_idx]->GetPhys(), 1, HWterm_2D_alpha, 1);
    Vmath::Smul(npts, this->alpha, HWterm_2D_alpha, 1, HWterm_2D_alpha, 1);
    Vmath::Vadd(npts, out_arr[w_idx], 1, HWterm_2D_alpha, 1, out_arr[w_idx], 1);
    Vmath::Vadd(npts, out_arr[ne_idx], 1, HWterm_2D_alpha, 1, out_arr[ne_idx],
                1);

    // Get poloidal field (y in HW equations)
    Array<OneD, NekDouble> E_poloidal(npts);
    Vmath::Vvtvvtp(npts, B_pol[0], 1, E[0]->GetPhys(), 1, B_pol[1], 1,
                   E[1]->GetPhys(), 1, E_poloidal, 1);
    // TODO normalise

    // Add -\kappa*\dpartial\phi/\dpartial y to RHS
    Vmath::Svtvp(npts, -this->kappa, E_poloidal, 1, out_arr[ne_idx], 1,
                 out_arr[ne_idx], 1);

    Array<OneD, MR::ExpListSharedPtr> diff_fields(nvariables);
    Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>(out_arr[i].size(), 0.0);
        diff_fields[i]  = m_fields[i];
    }
    CalcDiffTensor();
    CalcKappaTensor();
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
 * @brief Calls HelmSolve to solve for the electric potential
 *
 * @param in_arr Array of physical field values
 */
void ElectrostaticTurbulence::SolvePhi(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr)
{

    // Field indices
    int npts    = GetNpoints();
    int phi_idx = this->field_to_index["phi"];
    int w_idx   = this->field_to_index["w"];

    StdRegions::ConstFactorMap factors;
    // Helmholtz => Poisson (lambda = 0)
    factors[StdRegions::eFactorLambda] = 0.0;
    // Solve for phi. Output of this routine is in coefficient (spectral)
    // space, so backwards transform to physical space since we'll need that
    // for the advection step & computing drift velocity.
    m_fields[phi_idx]->HelmSolve(in_arr[w_idx],
                                 m_fields[phi_idx]->UpdateCoeffs(), factors,
                                 m_phi_varcoeff);
    m_fields[phi_idx]->BwdTrans(m_fields[phi_idx]->GetCoeffs(),
                                m_fields[phi_idx]->UpdatePhys());
}

/**
 * @brief Calculates initial potential and gradient
 */
void ElectrostaticTurbulence::CalcInitPhi()
{
    Array<OneD, Array<OneD, NekDouble>> in_arr(m_fields.size());
    for (auto ii = 0; ii < m_fields.size(); ii++)
    {
        in_arr[ii] = m_fields[ii]->GetPhys();
    }
    SolvePhi(in_arr);
    if (this->particles_enabled)
    {
        ComputeE();
    }
}

/**
 * @brief Compute the gradient of phi for evaluation at the particle positions.
 * (Stop-gap until NP's gradient-evaluate can be extended to 3D).
 */
void ElectrostaticTurbulence::ComputeE()
{
    int phi_idx = this->field_to_index["phi"];
    m_fields[phi_idx]->PhysDeriv(m_fields[phi_idx]->GetPhys(),
                                 E[0]->UpdatePhys(), E[1]->UpdatePhys(),
                                 E[2]->UpdatePhys());
    int npoints = m_fields[phi_idx]->GetTotPoints();
    Vmath::Smul(npoints, 1.0, E[0]->GetPhys(), 1, E[0]->UpdatePhys(), 1);
    Vmath::Smul(npoints, 1.0, E[1]->GetPhys(), 1, E[1]->UpdatePhys(), 1);
    Vmath::Smul(npoints, 1.0, E[2]->GetPhys(), 1, E[2]->UpdatePhys(), 1);

    E[0]->FwdTrans(E[0]->GetPhys(), E[0]->UpdateCoeffs());
    E[1]->FwdTrans(E[1]->GetPhys(), E[1]->UpdateCoeffs());
    E[2]->FwdTrans(E[2]->GetPhys(), E[2]->UpdateCoeffs());
}

/**
 *  @brief Compute components of advection velocities normal to trace elements
 * (faces, in 3D).
 *
 * @param[in,out] trace_vel_norm Trace normal velocities for each field
 * @param         adv_vel        Advection velocities for each field
 */
Array<OneD, NekDouble> &ElectrostaticTurbulence::GetAdvVelNorm(
    Array<OneD, NekDouble> &trace_vel_norm,
    const Array<OneD, Array<OneD, NekDouble>> &adv_vel)
{
    // Number of trace (interface) points
    int num_trace_pts = GetTraceNpoints();
    // Auxiliary variable to compute normal velocities
    Array<OneD, NekDouble> tmp(num_trace_pts);

    // Ensure output array is zeroed
    Vmath::Zero(num_trace_pts, trace_vel_norm, 1);

    // Compute advection vel dot trace normals and store
    for (int i = 0; i < adv_vel.size(); ++i)
    {
        m_fields[0]->ExtractTracePhys(adv_vel[i], tmp);
        Vmath::Vvtvp(num_trace_pts, m_traceNormals[i], 1, tmp, 1,
                     trace_vel_norm, 1, trace_vel_norm, 1);
    }
    return trace_vel_norm;
}

/**
 *  @brief Compute trace-normal advection velocities for the plasma density.
 */
Array<OneD, NekDouble> &ElectrostaticTurbulence::GetAdvVelNormElec()
{
    return GetAdvVelNorm(this->norm_vel_elec, this->ExB_vel);
}

/**
 *  @brief Construct flux array.
 *
 * @param  field_vals Physical values for each advection field
 * @param  adv_vel    Advection velocities for each advection field
 * @param[out] flux       Flux array
 */
void ElectrostaticTurbulence::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    const Array<OneD, Array<OneD, NekDouble>> &adv_vel,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    ASSERTL1(
        fluxes[0].size() == adv_vel.size(),
        "Dimension of flux array and advection velocity array do not match");
    int npts = field_vals[0].size();

    for (auto i = 0; i < fluxes.size(); ++i)
    {
        for (auto j = 0; j < fluxes[0].size(); ++j)
        {
            Vmath::Vmul(npts, field_vals[i], 1, adv_vel[j], 1, fluxes[i][j], 1);
        }
    }
}

/**
 * @brief Compute the flux vector for advection in the plasma density and
 * momentum equations.
 *
 * @param field_vals Array of Fields ptrs
 * @param[out] flux Resulting flux array
 */
void ElectrostaticTurbulence::GetFluxVectorElec(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    GetFluxVector(field_vals, this->ExB_vel, flux);
}

/**
 * @brief Construct the flux vector for the anisotropic diffusion problem.
 */
void ElectrostaticTurbulence::GetFluxVectorDiff(
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

/**
 * @brief After reading ICs, calculate phi and grad(phi)
 */
void ElectrostaticTurbulence::v_SetInitialConditions(NekDouble init_time,
                                                     bool dump_ICs,
                                                     const int domain)
{
    TokamakSystem::v_SetInitialConditions(init_time, dump_ICs, domain);
    CalcInitPhi();
}

void ElectrostaticTurbulence::load_params()
{
    TokamakSystem::load_params();
    // HW alpha
    m_session->LoadParameter("HW_alpha", this->alpha);
    // HW kappa
    m_session->LoadParameter("HW_kappa", this->kappa);
    int npoints = m_fields[0]->GetNpoints();

    /// Perp Laplace Tensor
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Array<OneD, NekDouble> phi(npoints, 0.0);
            for (int k = 0; k < npoints; k++)
            {
                phi[k] = -b_unit[i][k] * b_unit[j][k];
                if (i == j)
                {
                    phi[k] += 1;
                }
            }
            m_phi_varcoeff[vc[i][j]] = phi;
        }
    }
}
} // namespace NESO::Solvers::tokamak