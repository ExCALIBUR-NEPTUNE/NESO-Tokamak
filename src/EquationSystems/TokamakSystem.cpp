#include "TokamakSystem.hpp"
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <boost/core/ignore_unused.hpp>

namespace NESO::Solvers::tokamak
{

/// Name of class
static std::string class_name;
std::string TokamakSystem::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "Tokamak", TokamakSystem::create,
        "(2D) Tokamak equation system. Runs in either a 2D or 3D "
        "domain.");
/**
 * @brief Creates an instance of this class.
 */
static SU::EquationSystemSharedPtr create(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
{
    SU::EquationSystemSharedPtr p =
        MemoryManager<TokamakSystem>::AllocateSharedPtr(session, graph);
    p->InitObject();
    return p;
}

TokamakSystem::TokamakSystem(const LU::SessionReaderSharedPtr &session,
                             const SD::MeshGraphSharedPtr &graph)
    : TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>(session, graph),
      ExB_vel(graph->GetSpaceDimension())
{
    this->required_fld_names = {"T",        "ne",       "w",       "phi",
                                "gradphi0", "gradphi1", "gradphi2"};
    this->int_fld_names      = {"T", "ne", "w"};

    // Determine whether energy,enstrophy recording is enabled (output freq is
    // set)
    this->energy_enstrophy_recording_enabled =
        session->DefinesParameter("growth_rates_recording_step");
}

void TokamakSystem::calc_init_phi_and_gradphi()
{
    Array<OneD, Array<OneD, NekDouble>> in_arr(m_fields.size());
    for (auto ii = 0; ii < m_fields.size(); ii++)
    {
        in_arr[ii] = m_fields[ii]->GetPhys();
    }
    solve_phi(in_arr);
    if (this->particles_enabled)
    {
        compute_grad_phi();
    }
}

/**
 * @brief Compute the gradient of phi for evaluation at the particle positions.
 * (Stop-gap until NP's gradient-evaluate can be extended to 3D).
 */
void TokamakSystem::compute_grad_phi()
{
    int phi_idx      = this->field_to_index["phi"];
    int gradphi0_idx = this->field_to_index["gradphi0"];
    int gradphi1_idx = this->field_to_index["gradphi1"];
    int gradphi2_idx = this->field_to_index["gradphi2"];
    m_fields[phi_idx]->PhysDeriv(m_fields[phi_idx]->GetPhys(),
                                 m_fields[gradphi0_idx]->UpdatePhys(),
                                 m_fields[gradphi1_idx]->UpdatePhys(),
                                 m_fields[gradphi2_idx]->UpdatePhys());

    m_fields[gradphi0_idx]->FwdTrans(m_fields[gradphi0_idx]->GetPhys(),
                                     m_fields[gradphi0_idx]->UpdateCoeffs());
    m_fields[gradphi1_idx]->FwdTrans(m_fields[gradphi1_idx]->GetPhys(),
                                     m_fields[gradphi1_idx]->UpdateCoeffs());
    m_fields[gradphi2_idx]->FwdTrans(m_fields[gradphi2_idx]->GetPhys(),
                                     m_fields[gradphi2_idx]->UpdateCoeffs());
}
/**
 * @brief Perform projection into correct polynomial space.
 *
 * @details This routine projects the @p in_arr input and ensures the @p
 * out_arr output lives in the correct space. Since we are hard-coding DG,
 * this corresponds to a simple copy from in to out, since no elemental
 * connectivity is required and the output of the RHS function is
 * polynomial.
 *
 * @param in_arr Unprojected values
 * @param[out] out_arr Projected values
 * @param time Current simulation time
 *
 */
void TokamakSystem::DoOdeProjection(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time)
{
    int num_vars = in_arr.size();
    int npoints  = GetNpoints();
    for (int i = 0; i < num_vars; ++i)
    {
        Vmath::Vcopy(npoints, in_arr[i], 1, out_arr[i], 1);
    }
    SetBoundaryConditions(out_arr, time);
}

/**
 * @brief Populate rhs array ( @p out_arr ) for explicit time integration of
 * the 2D Hasegawa Wakatani equations.
 *
 * @param in_arr physical values of all fields
 * @param[out] out_arr output array (RHSs of time integration equations)
 */
void TokamakSystem::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time)
{

    // Get field indices
    int npts         = GetNpoints();
    int gradphi0_idx = this->field_to_index["gradphi0"];
    int gradphi1_idx = this->field_to_index["gradphi1"];
    int gradphi2_idx = this->field_to_index["gradphi2"];
    int ne_idx       = this->field_to_index["ne"];
    int phi_idx      = this->field_to_index["phi"];
    int w_idx        = this->field_to_index["w"];

    if (this->m_explicitAdvection)
    {
        zero_array_of_arrays(out_arr);
    }

    // Solve for electrostatic potential
    solve_phi(in_arr);
    // Calculate grad_phi
    compute_grad_phi();

    // Calculate ExB velocity
    // Array<OneD, NekDouble> inv_B_sq(npts);
    // Vmath::Vdiv(npts, 1.0, &mag_B[0], 1, &inv_B_sq[0], 1);
    // Vmath::Vdiv(npts, &inv_B_sq[0], 1, &mag_B[0], 1, &inv_B_sq[0], 1);

    // Vmath::Vmul(npts, this->B[1], 1, m_fields[gradphi2_idx]->GetPhys(), 1,
    //             this->ExB_vel[0], 1);
    Vmath::Vvtvvtm(npts, this->B[1], 1, m_fields[gradphi2_idx]->GetPhys(), 1,
                   this->B[2], 1, m_fields[gradphi1_idx]->GetPhys(), 1,
                   this->ExB_vel[0], 1);
    // Vmath::Vmul(npts, ExB_vel[0], 1, inv_B_sq, 1, ExB_vel[0], 1);

    // Vmath::Vmul(npts, this->B[0], 1, m_fields[gradphi2_idx]->GetPhys(), 1,
    //             this->ExB_vel[1], 1);
    Vmath::Vvtvvtm(npts, this->B[2], 1, m_fields[gradphi0_idx]->GetPhys(), 1,
                   this->B[0], 1, m_fields[gradphi2_idx]->GetPhys(), 1,
                   this->ExB_vel[1], 1);

    // Vmath::Vmul(npts, ExB_vel[1], 1, inv_B_sq, 1, ExB_vel[1], 1);

    if (this->m_graph->GetMeshDimension() == 3)
    {
        Vmath::Vvtvvtm(npts, this->B[0], 1, m_fields[gradphi1_idx]->GetPhys(),
                       1, this->B[1], 1, m_fields[gradphi0_idx]->GetPhys(), 1,
                       this->ExB_vel[2], 1);

        // Vmath::Vmul(npts, ExB_vel[2], 1, inv_B_sq, 1, ExB_vel[2], 1);
    }

    // Advect T, ne and w

    size_t nvariables = in_arr.size();
    size_t nTracePts  = GetTraceTotPoints();
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvariables);

    for (size_t i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        m_fields[i]->GetFwdBwdTracePhys(in_arr[i], Fwd[i], Bwd[i]);
    }

    this->adv_obj->Advect(nvariables, m_fields, this->ExB_vel, in_arr, out_arr,
                          time, Fwd, Bwd);
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

    // Get poloidal gradient (y in HW equations)
    Array<OneD, NekDouble> gradphi_poloidal(npts);
    Vmath::Vvtvvtp(npts, B_pol[0], 1, m_fields[gradphi0_idx]->GetPhys(), 1,
                   B_pol[1], 1, m_fields[gradphi1_idx]->GetPhys(), 1,
                   gradphi_poloidal, 1);
    // TODO normalise

    // Add -\kappa*\dpartial\phi/\dpartial y to RHS
    Vmath::Svtvp(npts, -this->kappa, gradphi_poloidal, 1, out_arr[ne_idx], 1,
                 out_arr[ne_idx], 1);

    if (m_explicitDiffusion)
    {
        Array<OneD, MultiRegions::ExpListSharedPtr> diff_fields(nvariables);
        Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvariables);
        for (int i = 0; i < nvariables; ++i)
        {
            outarrayDiff[i] = Array<OneD, NekDouble>(out_arr[i].size(), 0.0);
            diff_fields[i]  = m_fields[i];
        }

        m_diffusion->Diffuse(nvariables, diff_fields, in_arr, outarrayDiff, Fwd,
                             Bwd);

        for (int i = 0; i < nvariables; ++i)
        {
            Vmath::Vadd(out_arr[i].size(), outarrayDiff[i], 1, out_arr[i], 1,
                        out_arr[i], 1);
        }
    }
}

/**
 *  @brief Compute components of advection velocities normal to trace elements
 * (faces, in 3D).
 *
 * @param[in,out] trace_vel_norm Trace normal velocities for each field
 * @param         adv_vel        Advection velocities for each field
 */
Array<OneD, NekDouble> &TokamakSystem::get_adv_vel_norm(
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
Array<OneD, NekDouble> &TokamakSystem::get_adv_vel_norm_elec()
{
    return get_adv_vel_norm(this->norm_vel_elec, this->ExB_vel);
}

/**
 *  @brief Construct flux array.
 *
 * @param  field_vals Physical values for each advection field
 * @param  adv_vel    Advection velocities for each advection field
 * @param[out] flux       Flux array
 */
void TokamakSystem::get_flux_vector(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    const Array<OneD, Array<OneD, NekDouble>> &adv_vel,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    ASSERTL1(
        flux[0].size() == adv_vel.size(),
        "Dimension of flux array and advection velocity array do not match");
    int npts = field_vals[0].size();

    for (auto i = 0; i < flux.size(); ++i)
    {
        for (auto j = 0; j < flux[0].size(); ++j)
        {
            Vmath::Vmul(npts, field_vals[i], 1, adv_vel[j], 1, flux[i][j], 1);
        }
    }
}

/**
 * @brief Construct the flux vector for the anisotropic diffusion problem.
 */
void TokamakSystem::GetFluxVectorDiff(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscous_tensor)
{
    unsigned int nDim              = qfield.size();
    unsigned int nConvectiveFields = qfield[0].size();
    unsigned int nPts              = qfield[0][0].size();

    StdRegions::VarCoeffType vc[3][3] = {
        {StdRegions::eVarCoeffD00, StdRegions::eVarCoeffD01,
         StdRegions::eVarCoeffD02},
        {StdRegions::eVarCoeffD01, StdRegions::eVarCoeffD11,
         StdRegions::eVarCoeffD12},
        {StdRegions::eVarCoeffD02, StdRegions::eVarCoeffD12,
         StdRegions::eVarCoeffD22}};

    for (unsigned int i = 0; i < nConvectiveFields; ++i)
    {
        for (unsigned int j = 0; j < nDim; ++j)
        {
            Vmath::Vmul(nPts, m_varcoeff[vc[j][0]].GetValue(), 1, qfield[0][i],
                        1, viscous_tensor[j][i], 1);
            for (unsigned int k = 1; k < nDim; ++k)
            {
                Vmath::Vvtvp(nPts, m_varcoeff[vc[j][k]].GetValue(), 1,
                             qfield[k][i], 1, viscous_tensor[j][i], 1,
                             viscous_tensor[j][i], 1);
            }
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
void TokamakSystem::get_flux_vector_elec(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    get_flux_vector(field_vals, this->ExB_vel, flux);
}

/**
 * @brief Load all required session parameters into member variables.
 */
void TokamakSystem::load_params()
{
    TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>::load_params();

    // Type of advection to use -- in theory we also support flux reconstruction
    // for quad-based meshes, or you can use a standard convective term if you
    // were fully continuous in space. Default is DG.
    m_session->LoadSolverInfo("AdvectionType", this->adv_type, "WeakDG");

    int npoints  = m_fields[0]->GetNpoints();
    this->B      = Array<OneD, Array<OneD, NekDouble>>(3);
    this->B_pol  = Array<OneD, Array<OneD, NekDouble>>(3);
    this->b_unit = Array<OneD, Array<OneD, NekDouble>>(3);
    std::vector<std::string> Bstring;
    Bstring.push_back("Bx");
    Bstring.push_back("By");
    Bstring.push_back("Bz");
    Bstring.resize(3);

    if (m_session->DefinesFunction("MagneticMeanField"))
    {
        Array<OneD, Array<OneD, NekDouble>> B2D(3);
        GetFunction("MagneticMeanField")->Evaluate(Bstring, B2D);
        for (int d = 0; d < 3; ++d)
        {
            B[d]     = Array<OneD, NekDouble>(npoints, 0.0);
            B_pol[d] = Array<OneD, NekDouble>(npoints, 0.0);
        }

        if (this->m_graph->GetMeshDimension() == 3)
        {
            int npoints_2d = B2D[0].size();
            int increments = npoints / npoints_2d;
            Array<OneD, NekDouble> Bx(npoints_2d);
            Array<OneD, NekDouble> Bz(npoints_2d);
            for (int i = 0; i < increments; ++i)
            {
                Vmath::Vcopy(npoints_2d, &B2D[1][0], 1, &B[1][i * npoints_2d],
                             1);

                NekDouble theta = i * 2 * M_PI / increments;
                for (int j = 0; j < npoints_2d; ++j)
                {
                    Bx[j] = cos(theta) * B2D[0][j] - sin(theta) * B2D[2][j];
                    Bz[j] = cos(theta) * B2D[2][j] + sin(theta) * B2D[0][j];
                }
                Vmath::Vcopy(npoints_2d, &Bx[0], 1, &B[0][i * npoints_2d], 1);
                Vmath::Vcopy(npoints_2d, &Bz[0], 1, &B[2][i * npoints_2d], 1);

                Vmath::Vcopy(npoints_2d, &B2D[1][0], 1,
                             &B_pol[1][i * npoints_2d], 1);
                Vmath::Vcopy(npoints_2d, &B2D[0][0], 1,
                             &B_pol[0][i * npoints_2d], 1);
                Vmath::Vcopy(npoints_2d, &B2D[2][0], 1,
                             &B_pol[2][i * npoints_2d], 1);
            }
        }
        else
        {
            this->B     = B2D;
            this->B_pol = B2D;
        }
    }
    else if (m_session->DefinesFunction("MagneticField"))
    {
        GetFunction("MagneticField")->Evaluate(Bstring, B);
    }

    this->mag_B = Array<OneD, NekDouble>(npoints, 0.0);
    for (int i = 0; i < 3; ++i)
    {
        Vmath::Vvtvp(npoints, &this->B[i][0], 1, &this->B[i][0], 1, &mag_B[0],
                     1, &mag_B[0], 1);
    }

    for (auto idim = 0; idim < 3; idim++)
    {
        b_unit[idim] = Array<OneD, NekDouble>(npoints, 0.0);
        for (int k = 0; k < npoints; ++k)
        {
            this->b_unit[idim][k] =
                (this->mag_B[k] > 0) ? this->B[idim][k] / this->mag_B[k] : 0.0;
        }
    }

    /// Anisotropic Diffusivity Tensor
    m_session->LoadParameter("k_par", m_kpar, 100.0);
    m_session->LoadParameter("k_perp", m_kperp, 1.0);

    Array<OneD, NekDouble> d00(npoints, 1.0);
    Array<OneD, NekDouble> d01(npoints, 0.0);
    Array<OneD, NekDouble> d02(npoints, 0.0);
    Array<OneD, NekDouble> d10(npoints, 0.0);
    Array<OneD, NekDouble> d11(npoints, 1.0);
    Array<OneD, NekDouble> d12(npoints, 0.0);
    Array<OneD, NekDouble> d20(npoints, 0.0);
    Array<OneD, NekDouble> d21(npoints, 0.0);
    Array<OneD, NekDouble> d22(npoints, 1.0);

    for (int k = 0; k < npoints; k++)
    {
        d00[k] = (m_kpar - m_kperp) * B[0][k] * B[0][k] + m_kperp;
        d01[k] = (m_kpar - m_kperp) * B[0][k] * B[1][k];
        d02[k] = (m_kpar - m_kperp) * B[0][k] * B[2][k];
        d10[k] = (m_kpar - m_kperp) * B[1][k] * B[0][k];
        d11[k] = (m_kpar - m_kperp) * B[1][k] * B[1][k] + m_kperp;
        d12[k] = (m_kpar - m_kperp) * B[1][k] * B[2][k];
        d20[k] = (m_kpar - m_kperp) * B[2][k] * B[0][k];
        d21[k] = (m_kpar - m_kperp) * B[2][k] * B[1][k];
        d22[k] = (m_kpar - m_kperp) * B[2][k] * B[2][k] + m_kperp;
    }
    m_varcoeff[StdRegions::eVarCoeffD00] = d00;
    m_varcoeff[StdRegions::eVarCoeffD01] = d01;
    m_varcoeff[StdRegions::eVarCoeffD02] = d02;
    m_varcoeff[StdRegions::eVarCoeffD10] = d10;
    m_varcoeff[StdRegions::eVarCoeffD11] = d11;
    m_varcoeff[StdRegions::eVarCoeffD12] = d12;
    m_varcoeff[StdRegions::eVarCoeffD20] = d20;
    m_varcoeff[StdRegions::eVarCoeffD21] = d21;
    m_varcoeff[StdRegions::eVarCoeffD22] = d22;

    /// Anisotropic Diffusivity Tensor
    Array<OneD, NekDouble> phi_d00(npoints, 1.0);
    Array<OneD, NekDouble> phi_d01(npoints, 0.0);
    Array<OneD, NekDouble> phi_d02(npoints, 0.0);
    Array<OneD, NekDouble> phi_d10(npoints, 0.0);
    Array<OneD, NekDouble> phi_d11(npoints, 1.0);
    Array<OneD, NekDouble> phi_d12(npoints, 0.0);
    Array<OneD, NekDouble> phi_d20(npoints, 0.0);
    Array<OneD, NekDouble> phi_d21(npoints, 0.0);
    Array<OneD, NekDouble> phi_d22(npoints, 1.0);

    for (int k = 0; k < npoints; k++)
    {
        d00[k] = -B[0][k] * B[0][k] + 1;
        d01[k] = -B[0][k] * B[1][k];
        d02[k] = -B[0][k] * B[2][k];
        d10[k] = -B[1][k] * B[0][k];
        d11[k] = -B[1][k] * B[1][k] + 1;
        d12[k] = -B[1][k] * B[2][k];
        d20[k] = -B[2][k] * B[0][k];
        d21[k] = -B[2][k] * B[1][k];
        d22[k] = -B[2][k] * B[2][k] + 1;
    }

    m_phi_varcoeff[StdRegions::eVarCoeffD00] = phi_d00;
    m_phi_varcoeff[StdRegions::eVarCoeffD01] = phi_d01;
    m_phi_varcoeff[StdRegions::eVarCoeffD02] = phi_d02;
    m_phi_varcoeff[StdRegions::eVarCoeffD10] = phi_d10;
    m_phi_varcoeff[StdRegions::eVarCoeffD11] = phi_d11;
    m_phi_varcoeff[StdRegions::eVarCoeffD12] = phi_d12;
    m_phi_varcoeff[StdRegions::eVarCoeffD20] = phi_d20;
    m_phi_varcoeff[StdRegions::eVarCoeffD21] = phi_d21;
    m_phi_varcoeff[StdRegions::eVarCoeffD22] = phi_d22;

    // HW alpha (required)
    m_session->LoadParameter("HW_alpha", this->alpha);

    // HW kappa (required)
    m_session->LoadParameter("HW_kappa", this->kappa);

    // Type of Riemann solver to use. Default = "Upwind"
    m_session->LoadSolverInfo("UpwindType", this->riemann_solver_type,
                              "Upwind");

    // Particle-related parameters
    m_session->LoadParameter("particle_output_freq", particle_output_freq, 0);
    m_session->LoadParameter("num_particle_steps_per_fluid_step",
                             this->num_part_substeps, 1);
    this->part_timestep = m_timestep / this->num_part_substeps;
}

/**
 * @brief Calls HelmSolve to solve for the electric potential
 *
 * @param in_arr Array of physical field values
 */
void TokamakSystem::solve_phi(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr)
{

    // Field indices
    int npts    = GetNpoints();
    int w_idx   = this->field_to_index["w"];
    int phi_idx = this->field_to_index["phi"];
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

void TokamakSystem::v_GenerateSummary(SU::SummaryList &s)
{
    TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>::v_GenerateSummary(s);

    SU::AddSummaryItem(s, "Riemann solver", this->riemann_solver_type);

    SU::AddSummaryItem(s, "HW alpha", this->alpha);
    SU::AddSummaryItem(s, "HW kappa", this->kappa);
}

/**
 * @brief Post-construction class initialisation.
 *
 * @param create_field if true, create a new field object and add it to
 * m_fields. Optional, defaults to true.
 */
void TokamakSystem::v_InitObject(bool create_field)
{
    TimeEvoEqnSysBase::v_InitObject(create_field);

    // Since we are starting from a setup where each field is defined to be a
    // discontinuous field (and thus support DG), the first thing we do is to
    // recreate the phi field so that it is continuous, in order to support the
    // Poisson solve. Note that you can still perform a Poisson solve using a
    // discontinuous field, which is done via the hybridisable discontinuous
    // Galerkin (HDG) approach.
    int phi_idx       = this->field_to_index["phi"];
    m_fields[phi_idx] = MemoryManager<MR::ContField>::AllocateSharedPtr(
        m_session, m_graph, m_session->GetVariable(phi_idx), true, true);

    // Create storage for ExB vel
    int npts = GetNpoints();
    for (int i = 0; i < m_graph->GetSpaceDimension(); ++i)
    {
        this->ExB_vel[i] = Array<OneD, NekDouble>(npts, 0.0);
    }

    // Type of advection class to be used. By default, we only support the
    // discontinuous projection, since this is the only approach we're
    // considering for this solver.
    ASSERTL0(m_projectionType == MR::eDiscontinuous,
             "Unsupported projection type: only discontinuous"
             " projection supported.");

    if (m_explicitDiffusion)
    {
        std::string diffName;
        m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
        m_diffusion = SolverUtils::GetDiffusionFactory().CreateInstance(
            diffName, diffName);
        m_diffusion->SetFluxVector(&TokamakSystem::GetFluxVectorDiff, this);
        m_diffusion->InitObject(m_session, m_fields);
    }

    // Turn off forward-transform of initial conditions.
    m_homoInitialFwd = false;

    // Define the normal velocity fields.
    // These are populated at each step (by reference) in calls to GetVnAdv()
    if (m_fields[0]->GetTrace())
    {
        auto nTrace         = GetTraceNpoints();
        this->norm_vel_elec = Array<OneD, NekDouble>(nTrace);
    }

    // Advection objects
    // Need one per advection velocity
    this->adv_obj = SU::GetAdvectionFactory().CreateInstance(this->adv_type,
                                                             this->adv_type);

    // Set callback functions to compute flux vectors
    this->adv_obj->SetFluxVector(&TokamakSystem::get_flux_vector_elec, this);

    // Create Riemann solvers (one per advection object) and set normal velocity
    // callback functions
    this->riemann_solver = SU::GetRiemannSolverFactory().CreateInstance(
        this->riemann_solver_type, m_session);
    this->riemann_solver->SetScalar("Vn", &TokamakSystem::get_adv_vel_norm_elec,
                                    this);

    // Tell advection objects about the Riemann solvers and finish init
    this->adv_obj->SetRiemannSolver(this->riemann_solver);
    this->adv_obj->InitObject(m_session, m_fields);

    // Bind projection function for time integration object
    m_ode.DefineProjection(&TokamakSystem::DoOdeProjection, this);

    // Store DisContFieldSharedPtr casts of fields in a map, indexed by name
    int idx = 0;
    for (auto &field_name : m_session->GetVariables())
    {
        this->discont_fields[field_name] =
            std::dynamic_pointer_cast<MR::DisContField>(m_fields[idx]);
        idx++;
    }

    if (this->particles_enabled)
    {
        // Set up object to evaluate density field
        this->particle_sys->setup_evaluate_grad_phi(
            this->discont_fields.at("gradphi0"),
            this->discont_fields.at("gradphi1"),
            this->discont_fields.at("gradphi2"));
    }

    // Bind RHS function for explicit time integration
    m_ode.DefineOdeRhs(&TokamakSystem::DoOdeRhs, this);
    // Bind RHS function for implicit time integration
    m_ode.DefineImplicitSolve(&TokamakSystem::DoImplicitSolve, this);

    if (!m_explicitAdvection)
    {
        init_nonlin_sys_solver();
    }
    for (size_t i = 0; i < m_fields.size(); ++i)
    {

        bool Set = false;
        Array<OneD, const SD::BoundaryConditionShPtr> BndConds;
        Array<OneD, MR::ExpListSharedPtr> BndExp;
        size_t cnt = 0;
        BndConds   = m_fields[i]->GetBndConditions();
        BndExp     = m_fields[i]->GetBndCondExpansions();

        for (size_t n = 0; n < BndConds.size(); ++n)
        {
            std::string type = BndConds[n]->GetUserDefined();
            if (type.rfind("Oblique", 0) == 0)
            {
                ASSERTL0(
                    BndConds[n]->GetBoundaryConditionType() == SD::eRobin,
                    "Oblique boundary condition must be of type Robin <R>");
                // std::string::size_type indxBeg = type.find_first_of(':') + 1;
                // std::string fieldcomps =
                //     type.substr(indxBeg, std::string::npos);
            }
            if (BndConds[n]->GetBoundaryConditionType() == SD::ePeriodic)
            {
                continue;
            }

            if (!type.empty())
            {
                if (!Set)
                {
                    Array<OneD, Array<OneD, NekDouble>> magneticFieldBndElmt(3);

                    for (int d = 0; d < 3; d++)
                    {
                        m_fields[0]->ExtractPhysToBndElmt(
                            n, B[d], magneticFieldBndElmt[d]);
                    }
                    m_bndConds.push_back(
                        GetTokamakBndCondFactory().CreateInstance(
                            type, m_session, m_fields, magneticFieldBndElmt,
                            m_spacedim, n, cnt));
                    Set = true;
                }
            }
            cnt += BndExp[n]->GetExpSize();
        }
    }
    SetBoundaryConditionsBwdWeight();

    // Create diagnostic for recording growth rates
    if (this->energy_enstrophy_recording_enabled)
    {
        this->energy_enstrophy_recorder =
            std::make_shared<GrowthRatesRecorder<MultiRegions::DisContField>>(
                m_session, 2, this->discont_fields["ne"],
                this->discont_fields["w"], this->discont_fields["phi"],
                GetNpoints(), this->alpha, this->kappa);
    }
}

void TokamakSystem::init_nonlin_sys_solver()
{
    int npts_tot = 3 * m_fields[0]->GetNpoints();

    // Create the key to hold settings for nonlin solver
    LibUtilities::NekSysKey key = LibUtilities::NekSysKey();

    // Load required LinSys parameters:
    m_session->LoadParameter("NekLinSysMaxIterations",
                             key.m_NekLinSysMaxIterations, 30);
    m_session->LoadParameter("LinSysMaxStorage", key.m_LinSysMaxStorage, 30);
    m_session->LoadParameter("LinSysRelativeTolInNonlin",
                             key.m_NekLinSysTolerance, 5.0E-2);
    m_session->LoadParameter("GMRESMaxHessMatBand", key.m_KrylovMaxHessMatBand,
                             31);

    // Load required NonLinSys parameters:
    m_session->LoadParameter("jacobi_free_eps", this->jacobi_free_eps, 5.0E-8);
    m_session->LoadParameter("NekNonlinSysMaxIterations",
                             key.m_NekNonlinSysMaxIterations, 10);
    m_session->LoadParameter("NewtonRelativeIteTol",
                             key.m_NekNonLinSysTolerance, 1.0E-12);
    WARNINGL0(!m_session->DefinesParameter("NewtonAbsoluteIteTol"),
              "Please specify NewtonRelativeIteTol instead of "
              "NewtonAbsoluteIteTol in XML session file");
    m_session->LoadParameter("NonlinIterTolRelativeL2",
                             key.m_NonlinIterTolRelativeL2, 1.0E-3);
    m_session->LoadSolverInfo("LinSysIterSolverTypeInNonlin",
                              key.m_LinSysIterSolverTypeInNonlin, "GMRES");

    LibUtilities::NekSysOperators nekSysOp;
    nekSysOp.DefineNekSysResEval(&TokamakSystem::nonlin_sys_evaluator_1D, this);
    nekSysOp.DefineNekSysLhsEval(&TokamakSystem::matrix_multiply_matrix_free,
                                 this);
    nekSysOp.DefineNekSysPrecon(&TokamakSystem::do_null_precon, this);

    // Initialize non-linear system
    this->nonlin_sys =
        LibUtilities::GetNekNonlinSysIterFactory().CreateInstance(
            "Newton", m_session, m_comm->GetRowComm(), npts_tot, key);
    this->nonlin_sys->SetSysOperators(nekSysOp);
}

void TokamakSystem::DoImplicitSolve(
    const Array<OneD, const Array<OneD, NekDouble>> &inpnts,
    Array<OneD, Array<OneD, NekDouble>> &outpnt, const NekDouble time,
    const NekDouble lambda)
{
    this->time_int_lambda   = lambda;
    this->bnd_evaluate_time = time;
    unsigned int npoints    = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> in_arr(3 * npoints);
    Array<OneD, NekDouble> outarray(3 * npoints, 0.0);
    Array<OneD, NekDouble> tmp;

    for (int i = 0; i < 3; ++i)
    {
        int noffset = i * npoints;
        Vmath::Vcopy(npoints, inpnts[i], 1, tmp = in_arr + noffset, 1);
    }

    implicit_time_int_1D(in_arr, outarray);

    for (int i = 0; i < 3; ++i)
    {
        int noffset = i * npoints;
        Vmath::Vcopy(npoints, outarray + noffset, 1, outpnt[i], 1);
    }
}

void TokamakSystem::implicit_time_int_1D(
    const Array<OneD, const NekDouble> &in_arr, Array<OneD, NekDouble> &out)
{
    calc_ref_vals(in_arr);

    this->nonlin_sys->SetRhsMagnitude(this->in_arr_norm);

    this->newton_its_counter +=
        this->nonlin_sys->SolveSystem(in_arr.size(), in_arr, out, 0);

    this->lin_its_counter += this->nonlin_sys->GetNtotLinSysIts();

    this->imp_stages_counter++;
}

void TokamakSystem::calc_ref_vals(const Array<OneD, const NekDouble> &in_arr)
{
    unsigned int npoints = m_fields[0]->GetNpoints();

    Array<OneD, NekDouble> magnitdEstimat(3, 0.0);

    for (int i = 0; i < 3; ++i)
    {
        int offset = i * npoints;
        magnitdEstimat[i] =
            Vmath::Dot(npoints, in_arr + offset, in_arr + offset);
    }
    m_comm->GetSpaceComm()->AllReduce(magnitdEstimat,
                                      Nektar::LibUtilities::ReduceSum);

    this->in_arr_norm = 0.0;
    for (int i = 0; i < 3; ++i)
    {
        this->in_arr_norm += magnitdEstimat[i];
    }
}

void TokamakSystem::nonlin_sys_evaluator_1D(
    const Array<OneD, const NekDouble> &in_arr, Array<OneD, NekDouble> &out,
    [[maybe_unused]] const bool &flag)
{
    unsigned int npoints = m_fields[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble>> in2D(3);
    Array<OneD, Array<OneD, NekDouble>> out2D(3);
    for (int i = 0; i < 3; ++i)
    {
        int offset = i * npoints;
        in2D[i]    = in_arr + offset;
        out2D[i]   = out + offset;
    }
    nonlin_sys_evaluator(in2D, out2D);
}

void TokamakSystem::nonlin_sys_evaluator(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out)
{
    unsigned int npoints = m_fields[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble>> inpnts(3);
    for (int i = 0; i < 3; ++i)
    {
        inpnts[i] = Array<OneD, NekDouble>(npoints, 0.0);
    }

    DoOdeProjection(in_arr, inpnts, this->bnd_evaluate_time);
    DoOdeRhs(inpnts, out, this->bnd_evaluate_time);

    for (int i = 0; i < 3; ++i)
    {
        Vmath::Svtvp(npoints, -this->time_int_lambda, out[i], 1, in_arr[i], 1,
                     out[i], 1);
        Vmath::Vsub(npoints, out[i], 1,
                    this->nonlin_sys->GetRefSourceVec() + i * npoints, 1,
                    out[i], 1);
    }
}

void TokamakSystem::matrix_multiply_matrix_free(
    const Array<OneD, const NekDouble> &in_arr, Array<OneD, NekDouble> &out,
    [[maybe_unused]] const bool &flag)
{
    const Array<OneD, const NekDouble> solref =
        this->nonlin_sys->GetRefSolution();
    const Array<OneD, const NekDouble> resref =
        this->nonlin_sys->GetRefResidual();

    unsigned int ntotal   = in_arr.size();
    NekDouble magninarray = Vmath::Dot(ntotal, in_arr, in_arr);
    m_comm->GetSpaceComm()->AllReduce(magninarray,
                                      Nektar::LibUtilities::ReduceSum);
    NekDouble eps = this->jacobi_free_eps *
                    sqrt((sqrt(this->in_arr_norm) + 1.0) / magninarray);

    Array<OneD, NekDouble> solplus{ntotal};
    Array<OneD, NekDouble> resplus{ntotal};

    Vmath::Svtvp(ntotal, eps, in_arr, 1, solref, 1, solplus, 1);
    nonlin_sys_evaluator_1D(solplus, resplus, flag);
    Vmath::Vsub(ntotal, resplus, 1, resref, 1, out, 1);
    Vmath::Smul(ntotal, 1.0 / eps, out, 1, out, 1);
}

void TokamakSystem::do_null_precon(const Array<OneD, const NekDouble> &in_arr,
                                   Array<OneD, NekDouble> &outarray,
                                   [[maybe_unused]] const bool &flag)
{
    Vmath::Vcopy(in_arr.size(), in_arr, 1, outarray, 1);
}

/**
 * @brief Override v_PostIntegrate to do particle output and compute
 * diagnostics, if enabled
 * @param step Time step number
 */
bool TokamakSystem::v_PostIntegrate(int step)
{
    if (energy_enstrophy_recording_enabled)
    {
        this->energy_enstrophy_recorder->compute(step);
    }

    this->solver_callback_handler.call_post_integrate(this);

    // Writes a step of the particle trajectory.
    if (this->particles_enabled && particle_output_freq > 0 &&
        (step % particle_output_freq) == 0)
    {
        this->particle_sys->write(step);
    }
    return TimeEvoEqnSysBase<SU::UnsteadySystem,
                             ParticleSystem>::v_PostIntegrate(step);
}

/**
 * @brief Override v_PreIntegrate to do particle system integration, and initial
 * set up for mass recording diagnostic (first call only)
 *
 * @param step Time step number
 */
bool TokamakSystem::v_PreIntegrate(int step)
{
    this->solver_callback_handler.call_pre_integrate(this);

    if (this->particles_enabled)
    {
        // Integrate the particle system to the requested time.
        this->particle_sys->evaluate_fields();
        this->particle_sys->integrate(m_time + m_timestep, this->part_timestep);
    }

    return UnsteadySystem::v_PreIntegrate(step);
}

/**
 * @brief After reading ICs, calculate phi and grad(phi)
 */
void TokamakSystem::v_SetInitialConditions(NekDouble init_time, bool dump_ICs,
                                           const int domain)
{
    TimeEvoEqnSysBase<SU::UnsteadySystem,
                      ParticleSystem>::v_SetInitialConditions(init_time,
                                                              dump_ICs, domain);
    calc_init_phi_and_gradphi();
    if (this->particle_sys)
    {
        this->particle_sys->initialise_particles_from_fields();
    }
}

void TokamakSystem::SetBoundaryConditions(
    Array<OneD, Array<OneD, NekDouble>> &physarray, NekDouble time)
{
    if (!m_bndConds.empty())
    {
        // Loop over user-defined boundary conditions
        for (auto &x : m_bndConds)
        {
            x->Apply(physarray, time);
        }
    }
}

void TokamakSystem::SetBoundaryConditionsBwdWeight()
{
    if (m_bndConds.size())
    {
        // Loop over user-defined boundary conditions
        for (auto &x : m_bndConds)
        {
            x->ApplyBwdWeight();
        }
    }
}

} // namespace NESO::Solvers::tokamak
