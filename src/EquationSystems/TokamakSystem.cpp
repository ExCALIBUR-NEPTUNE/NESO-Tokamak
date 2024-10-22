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
    this->required_fld_names = {"phi", "w", "Te", "ne"};
    this->int_fld_names      = {"w", "Te", "ne"};

    // Determine whether energy,enstrophy recording is enabled (output freq is
    // set)
    this->energy_enstrophy_recording_enabled =
        session->DefinesParameter("growth_rates_recording_step");
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

    int npoints = m_fields[0]->GetNpoints();
    Array<OneD, MR::ExpListSharedPtr> B_in(3);
    this->B_pol  = Array<OneD, Array<OneD, NekDouble>>(3);
    this->b_unit = Array<OneD, Array<OneD, NekDouble>>(3);
    std::vector<std::string> Bstring;
    Bstring.push_back("Bx");
    Bstring.push_back("By");
    Bstring.push_back("Bz");
    Bstring.resize(3);

    if (m_session->DefinesFunction("MagneticMeanField"))
    {
        if (this->m_graph->GetMeshDimension() == 2)
        {
            GetFunction("MagneticMeanField")->Evaluate(Bstring, B_in);
            for (int d = 0; d < 3; ++d)
            {
                this->B_pol[d] = B[d]->GetPhys();
            }
        }
        else if (this->m_graph->GetMeshDimension() == 3)
        {
            Array<OneD, Array<OneD, NekDouble>> B2D(3);
            GetFunction("MagneticMeanField")->Evaluate(Bstring, B2D);

            Array<OneD, Array<OneD, NekDouble>> B3D(3);
            for (int d = 0; d < 3; ++d)
            {
                B3D[d]   = Array<OneD, NekDouble>(npoints, 0.0);
                B_pol[d] = Array<OneD, NekDouble>(npoints, 0.0);
            }

            int npoints_2d = B2D[0].size();
            int increments = npoints / npoints_2d;
            Array<OneD, NekDouble> Bx(npoints_2d);
            Array<OneD, NekDouble> Bz(npoints_2d);
            for (int i = 0; i < increments; ++i)
            {
                Vmath::Vcopy(npoints_2d, &B2D[1][0], 1, &B3D[1][i * npoints_2d],
                             1);

                NekDouble theta = i * 2 * M_PI / increments;
                for (int j = 0; j < npoints_2d; ++j)
                {
                    Bx[j] = cos(theta) * B2D[0][j] - sin(theta) * B2D[2][j];
                    Bz[j] = cos(theta) * B2D[2][j] + sin(theta) * B2D[0][j];
                }
                Vmath::Vcopy(npoints_2d, &Bx[0], 1, &B3D[0][i * npoints_2d], 1);
                Vmath::Vcopy(npoints_2d, &Bz[0], 1, &B3D[2][i * npoints_2d], 1);

                Vmath::Vcopy(npoints_2d, &B2D[1][0], 1,
                             &B_pol[1][i * npoints_2d], 1);
                Vmath::Vcopy(npoints_2d, &B2D[0][0], 1,
                             &B_pol[0][i * npoints_2d], 1);
                Vmath::Vcopy(npoints_2d, &B2D[2][0], 1,
                             &B_pol[2][i * npoints_2d], 1);
            }
        }
    }
    else if (m_session->DefinesFunction("MagneticField"))
    {
        GetFunction("MagneticField")->Evaluate(Bstring, B_in);
    }
    this->B = Array<OneD, MR::DisContFieldSharedPtr>(3);
    for (int d = 0; d < 3; ++d)
    {
        B[d] = std::static_pointer_cast<MR::DisContField>(B_in[d]);
    }

    m_session->LoadSolverInfo("Mode", m_mode, "MeanFluid");

    this->mag_B = Array<OneD, NekDouble>(npoints, 0.0);
    for (int i = 0; i < 3; ++i)
    {
        Vmath::Vvtvp(npoints, B[i]->GetPhys(), 1, B[i]->GetPhys(), 1, mag_B, 1,
                     mag_B, 1);
    }

    for (auto idim = 0; idim < 3; idim++)
    {
        b_unit[idim] = Array<OneD, NekDouble>(npoints, 0.0);
        for (int k = 0; k < npoints; ++k)
        {
            this->b_unit[idim][k] = (this->mag_B[k] > 0)
                                        ? B[idim]->GetPhys()[k] / this->mag_B[k]
                                        : 0.0;
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
        d00[k] = (m_kpar - m_kperp) * b_unit[0][k] * b_unit[0][k] + m_kperp;
        d01[k] = (m_kpar - m_kperp) * b_unit[0][k] * b_unit[1][k];
        d02[k] = (m_kpar - m_kperp) * b_unit[0][k] * b_unit[2][k];
        d10[k] = (m_kpar - m_kperp) * b_unit[1][k] * b_unit[0][k];
        d11[k] = (m_kpar - m_kperp) * b_unit[1][k] * b_unit[1][k] + m_kperp;
        d12[k] = (m_kpar - m_kperp) * b_unit[1][k] * b_unit[2][k];
        d20[k] = (m_kpar - m_kperp) * b_unit[2][k] * b_unit[0][k];
        d21[k] = (m_kpar - m_kperp) * b_unit[2][k] * b_unit[1][k];
        d22[k] = (m_kpar - m_kperp) * b_unit[2][k] * b_unit[2][k] + m_kperp;
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
        d00[k] = -b_unit[0][k] * b_unit[0][k] + 1;
        d01[k] = -b_unit[0][k] * b_unit[1][k];
        d02[k] = -b_unit[0][k] * b_unit[2][k];
        d10[k] = -b_unit[1][k] * b_unit[0][k];
        d11[k] = -b_unit[1][k] * b_unit[1][k] + 1;
        d12[k] = -b_unit[1][k] * b_unit[2][k];
        d20[k] = -b_unit[2][k] * b_unit[0][k];
        d21[k] = -b_unit[2][k] * b_unit[1][k];
        d22[k] = -b_unit[2][k] * b_unit[2][k] + 1;
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

    m_session->LoadParameter("lambda", this->lambda);

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
 * @brief Calculates initial potential and gradient
 */
void TokamakSystem::CalcInitPhi()
{
    Array<OneD, Array<OneD, NekDouble>> in_arr(m_fields.size());
    for (auto ii = 0; ii < m_fields.size(); ii++)
    {
        in_arr[ii] = m_fields[ii]->GetPhys();
    }
    SolvePhi(in_arr);
    if (this->particles_enabled)
    {
        ComputeGradPhi();
    }
}

/**
 * @brief Calls HelmSolve to solve for the electric potential
 *
 * @param in_arr Array of physical field values
 */
void TokamakSystem::SolvePhi(
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
 * @brief Compute the gradient of phi for evaluation at the particle positions.
 * (Stop-gap until NP's gradient-evaluate can be extended to 3D).
 */
void TokamakSystem::ComputeGradPhi()
{
    int phi_idx = this->field_to_index["phi"];
    m_fields[phi_idx]->PhysDeriv(
        m_fields[phi_idx]->GetPhys(), m_grad_phi[0]->UpdatePhys(),
        m_grad_phi[1]->UpdatePhys(), m_grad_phi[2]->UpdatePhys());

    m_grad_phi[0]->FwdTrans(m_grad_phi[0]->GetPhys(),
                            m_grad_phi[0]->UpdateCoeffs());
    m_grad_phi[1]->FwdTrans(m_grad_phi[1]->GetPhys(),
                            m_grad_phi[1]->UpdateCoeffs());
    m_grad_phi[2]->FwdTrans(m_grad_phi[2]->GetPhys(),
                            m_grad_phi[2]->UpdateCoeffs());
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
 * @brief Populate rhs array ( @p out_arr ) where mean fluid
 * equations are used for the advection velocities.
 *
 * @param in_arr physical values of all fields
 * @param[out] out_arr output array (RHSs of time integration equations)
 */
void TokamakSystem::DoOdeRhsMF(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time)
{
    // Get field indices
    int npts    = GetNpoints();
    int phi_idx = this->field_to_index["phi"];
    int Te_idx  = this->field_to_index["Te"];
    int ne_idx  = this->field_to_index["ne"];

    for (int i = 0; i < npts; ++i)
    {
        m_fields[phi_idx]->UpdatePhys()[i] =
            this->lambda + m_fields[Te_idx]->GetPhys()[i] *
                               std::log(m_fields[ne_idx]->GetPhys()[i]);
    }

    if (this->m_explicitAdvection)
    {
        zero_array_of_arrays(out_arr);
    }

    Array<OneD, Array<OneD, NekDouble>> gradn(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        gradn[i] = Array<OneD, NekDouble>(npts, 0.0);
    }
    if (m_spacedim == 2)
    {
        m_fields[phi_idx]->PhysDeriv(m_fields[ne_idx]->GetPhys(), gradn[0],
                                     gradn[1]);
    }
    else if (m_spacedim == 3)
    {
        m_fields[ne_idx]->PhysDeriv(m_fields[ne_idx]->GetPhys(), gradn[0],
                                    gradn[1], gradn[2]);
    }
    for (int i = 0; i < m_spacedim; ++i)
    {
        Vmath::Smul(npts, m_kperp, gradn[i], 1, this->ExB_vel[i], 1);
    }

    size_t nvariables = in_arr.size();

    m_advection->Advect(nvariables, m_fields, this->ExB_vel, in_arr, out_arr,
                        time);
    for (auto i = 0; i < nvariables; ++i)
    {
        Vmath::Neg(npts, out_arr[i], 1);
    }

    Array<OneD, MultiRegions::ExpListSharedPtr> diff_fields(nvariables);
    Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>(out_arr[i].size(), 0.0);
        diff_fields[i]  = m_fields[i];
    }

    m_diffusion->Diffuse(nvariables, diff_fields, in_arr, outarrayDiff);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(out_arr[i].size(), outarrayDiff[i], 1, out_arr[i], 1,
                    out_arr[i], 1);
    }
}

/**
 * @brief Populate rhs array ( @p out_arr ) where electrostatic turbulence
 * equations (Hasegana-Wakatani) are used for the advection velocities.
 *
 * @param in_arr physical values of all fields
 * @param[out] out_arr output array (RHSs of time integration equations)
 */
void TokamakSystem::DoOdeRhsET(
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
    // Calculate grad_phi
    ComputeGradPhi();

    // Calculate ExB velocity
    Vmath::Vvtvvtm(npts, B[1]->GetPhys(), 1, m_grad_phi[2]->GetPhys(), 1,
                   B[2]->GetPhys(), 1, m_grad_phi[1]->GetPhys(), 1,
                   this->ExB_vel[0], 1);

    Vmath::Vvtvvtm(npts, B[2]->GetPhys(), 1, m_grad_phi[0]->GetPhys(), 1,
                   B[0]->GetPhys(), 1, m_grad_phi[2]->GetPhys(), 1,
                   this->ExB_vel[1], 1);

    if (this->m_graph->GetMeshDimension() == 3)
    {
        Vmath::Vvtvvtm(npts, B[0]->GetPhys(), 1, m_grad_phi[1]->GetPhys(), 1,
                       B[1]->GetPhys(), 1, m_grad_phi[0]->GetPhys(), 1,
                       this->ExB_vel[2], 1);
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

    // Get poloidal gradient (y in HW equations)
    Array<OneD, NekDouble> gradphi_poloidal(npts);
    Vmath::Vvtvvtp(npts, B_pol[0], 1, m_grad_phi[0]->GetPhys(), 1, B_pol[1], 1,
                   m_grad_phi[1]->GetPhys(), 1, gradphi_poloidal, 1);
    // TODO normalise

    // Add -\kappa*\dpartial\phi/\dpartial y to RHS
    Vmath::Svtvp(npts, -this->kappa, gradphi_poloidal, 1, out_arr[ne_idx], 1,
                 out_arr[ne_idx], 1);

    Array<OneD, MultiRegions::ExpListSharedPtr> diff_fields(nvariables);
    Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>(out_arr[i].size(), 0.0);
        diff_fields[i]  = m_fields[i];
    }

    m_diffusion->Diffuse(nvariables, diff_fields, in_arr, outarrayDiff);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(out_arr[i].size(), outarrayDiff[i], 1, out_arr[i], 1,
                    out_arr[i], 1);
    }
}

/**
 *  @brief Compute components of advection velocities normal to trace elements
 * (faces, in 3D).
 *
 * @param[in,out] trace_vel_norm Trace normal velocities for each field
 * @param         adv_vel        Advection velocities for each field
 */
Array<OneD, NekDouble> &TokamakSystem::GetAdvVelNorm(
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
Array<OneD, NekDouble> &TokamakSystem::GetAdvVelNormElec()
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
void TokamakSystem::GetFluxVector(
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
 * @brief Compute the flux vector for advection in the plasma density and
 * momentum equations.
 *
 * @param field_vals Array of Fields ptrs
 * @param[out] flux Resulting flux array
 */
void TokamakSystem::GetFluxVectorElec(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux)
{
    GetFluxVector(field_vals, this->ExB_vel, flux);
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

void TokamakSystem::v_GenerateSummary(SU::SummaryList &s)
{
    TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>::v_GenerateSummary(s);

    SU::AddSummaryItem(s, "Riemann solver", this->riemann_solver_type);
    SU::AddSummaryItem(s, "Turbulence mode", this->m_mode);
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

    m_grad_phi = Array<OneD, MR::DisContFieldSharedPtr>(m_spacedim);
    for (int i = 0; i < m_spacedim; ++i)
    {
        m_grad_phi[i] = MemoryManager<MR::DisContField>::AllocateSharedPtr(
            *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]));
    }

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

    // Setup diffusion object
    std::string diffName;
    m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
    m_diffusion =
        SolverUtils::GetDiffusionFactory().CreateInstance(diffName, diffName);
    m_diffusion->SetFluxVector(&TokamakSystem::GetFluxVectorDiff, this);
    m_diffusion->InitObject(m_session, m_fields);

    // Turn off forward-transform of initial conditions.
    m_homoInitialFwd = false;

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
    this->riemann_solver->SetScalar("Vn", &TokamakSystem::GetAdvVelNormElec,
                                    this);

    // Setup advection object
    m_advection = SU::GetAdvectionFactory().CreateInstance(this->adv_type,
                                                           this->adv_type);
    m_advection->SetFluxVector(&TokamakSystem::GetFluxVectorElec, this);
    m_advection->SetRiemannSolver(this->riemann_solver);
    m_advection->InitObject(m_session, m_fields);

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
            m_grad_phi[0], m_grad_phi[1], m_grad_phi[2]);
        this->particle_sys->setup_evaluate_B(B[0], B[1], B[2]);
    }

    // Bind RHS function for explicit time integration
    if (m_mode == "MeanFluid")
    {
        m_ode.DefineOdeRhs(&TokamakSystem::DoOdeRhsMF, this);
    }
    else if (m_mode == "ElectrostaticTurbulent")
    {
        m_ode.DefineOdeRhs(&TokamakSystem::DoOdeRhsET, this);
    }

    // Bind projection function for time integration object
    m_ode.DefineProjection(&TokamakSystem::DoOdeProjection, this);

    if (!m_explicitAdvection)
    {
        m_implHelper =
            std::make_shared<ImplicitHelper>(m_session, m_fields, m_ode, 3);
        m_implHelper->InitialiseNonlinSysSolver();
        m_ode.DefineImplicitSolve(&ImplicitHelper::ImplicitTimeInt,
                                  m_implHelper);
    }

    // Setup boundary conditions
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
            }
            if (type.rfind("Sheath", 0) == 0)
            {
                ASSERTL0(BndConds[n]->GetBoundaryConditionType() == SD::eRobin,
                         "Sheath boundary condition must be of type Robin <R>");
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
                            n, b_unit[d], magneticFieldBndElmt[d]);
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
    CalcInitPhi();
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
