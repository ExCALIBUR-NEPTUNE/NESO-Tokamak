#include "HWITSystem.hpp"
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <LibUtilities/TimeIntegration/TimeIntegrationScheme.h>
#include <boost/core/ignore_unused.hpp>

namespace NESO::Solvers::hw_impurity_transport {

/// Name of class
static std::string class_name;
std::string HWITSystem::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "HWIT", HWITSystem::create,
        "(2D) Hasegawa-Wakatani equation system. Runs in either a 2D or 3D "
        "domain.");
/**
 * @brief Creates an instance of this class.
 */
static SU::EquationSystemSharedPtr
create(const LU::SessionReaderSharedPtr &session,
       const SD::MeshGraphSharedPtr &graph) {
  SU::EquationSystemSharedPtr p =
      MemoryManager<HWITSystem>::AllocateSharedPtr(session, graph);
  p->InitObject();
  return p;
}

HWITSystem::HWITSystem(const LU::SessionReaderSharedPtr &session,
                       const SD::MeshGraphSharedPtr &graph)
    : TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>(session, graph),
      m_ExB_vel(graph->GetSpaceDimension()) {
  this->required_fld_names = {"ne",       "w",        "phi",
                              "gradphi0", "gradphi1", "gradphi2"};
  this->int_fld_names = {"ne", "w"};

  // Frequency of growth rate recording. Set zero to disable.
  m_diag_growth_rates_recording_enabled =
      session->DefinesParameter("growth_rates_recording_step");

  // Frequency of mass recording. Set zero to disable.
  m_diag_mass_recording_enabled =
      session->DefinesParameter("mass_recording_step");
}

void HWITSystem::calc_init_phi_and_gradphi() {
  Array<OneD, Array<OneD, NekDouble>> in_arr(m_fields.size());
  for (auto ii = 0; ii < m_fields.size(); ii++) {
    in_arr[ii] = m_fields[ii]->GetPhys();
  }
  solve_phi(in_arr);
  if (this->particles_enabled) {
    compute_grad_phi();
  }
}

/**
 * @brief Compute the gradient of phi for evaluation at the particle positions.
 * (Stop-gap until NP's gradient-evaluate can be extended to 3D).
 */
void HWITSystem::compute_grad_phi() {
  int phi_idx = m_field_to_index.get_idx("phi");
  int gradphi0_idx = m_field_to_index.get_idx("gradphi0");
  int gradphi1_idx = m_field_to_index.get_idx("gradphi1");
  int gradphi2_idx = m_field_to_index.get_idx("gradphi2");
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
void HWITSystem::do_ode_projection(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time) {
  int num_vars = in_arr.size();
  int npoints = GetNpoints();
  SetBoundaryConditions(time);
  for (int i = 0; i < num_vars; ++i) {
    Vmath::Vcopy(npoints, in_arr[i], 1, out_arr[i], 1);
  }
}

/**
 * @brief Populate rhs array ( @p out_arr ) for explicit time integration of
 * the 2D Hasegawa Wakatani equations.
 *
 * @param in_arr physical values of all fields
 * @param[out] out_arr output array (RHSs of time integration equations)
 */
void HWITSystem::explicit_time_int(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time) {

  // Get field indices
  int npts = GetNpoints();
  int gradphi0_idx = m_field_to_index.get_idx("gradphi0");
  int gradphi1_idx = m_field_to_index.get_idx("gradphi1");
  int gradphi2_idx = m_field_to_index.get_idx("gradphi2");
  int ne_idx = m_field_to_index.get_idx("ne");
  int phi_idx = m_field_to_index.get_idx("phi");
  int w_idx = m_field_to_index.get_idx("w");
  // // Check in_arr for NaNs - SLLLLLLOOOOOW
  // for (auto &var : {"ne", "w"}) {
  //   auto fidx = m_field_to_index.get_idx(var);
  //   for (auto ii = 0; ii < in_arr[fidx].size(); ii++) {
  //     std::stringstream err_msg;
  //     err_msg << "Found NaN in field " << var;
  //     NESOASSERT(std::isfinite(in_arr[fidx][ii]), err_msg.str().c_str());
  //   }
  // }

  if (this->m_explicitAdvection) {
    zero_out_array(out_arr);
  }

  // Solve for electrostatic potential
  solve_phi(in_arr);
  // Calculate grad_phi
  compute_grad_phi();

  // v_ExB = grad(phi) X B / |B|^2
  Vmath::Svtsvtp(npts, -m_B[2] / m_Bmag / m_Bmag,
                 m_fields[gradphi1_idx]->GetPhys(), 1, m_B[1] / m_Bmag / m_Bmag,
                 m_fields[gradphi2_idx]->GetPhys(), 1, m_ExB_vel[0], 1);
  Vmath::Svtsvtp(npts, -m_B[0] / m_Bmag / m_Bmag,
                 m_fields[gradphi2_idx]->GetPhys(), 1, m_B[2] / m_Bmag / m_Bmag,
                 m_fields[gradphi0_idx]->GetPhys(), 1, m_ExB_vel[1], 1);
  if (this->m_graph->GetMeshDimension() == 3) {

    Vmath::Svtsvtp(npts, -m_B[1] / m_Bmag / m_Bmag,
                   m_fields[gradphi0_idx]->GetPhys(), 1,
                   m_B[0] / m_Bmag / m_Bmag, m_fields[gradphi1_idx]->GetPhys(),
                   1, m_ExB_vel[2], 1);
  }
  // Advect ne and w
  this->adv_obj->Advect(2, m_fields, m_ExB_vel, in_arr, out_arr, time);
  for (auto i = 0; i < 2; ++i) {
    Vmath::Neg(npts, out_arr[i], 1);
  }

  // Add \alpha*(\phi-ne) to RHS
  Array<OneD, NekDouble> HWterm_2D_alpha(npts);
  Vmath::Vsub(npts, m_fields[phi_idx]->GetPhys(), 1,
              m_fields[ne_idx]->GetPhys(), 1, HWterm_2D_alpha, 1);
  Vmath::Smul(npts, m_alpha, HWterm_2D_alpha, 1, HWterm_2D_alpha, 1);
  Vmath::Vadd(npts, out_arr[w_idx], 1, HWterm_2D_alpha, 1, out_arr[w_idx], 1);
  Vmath::Vadd(npts, out_arr[ne_idx], 1, HWterm_2D_alpha, 1, out_arr[ne_idx], 1);

  // Add -\kappa*\dpartial\phi/\dpartial y to RHS
  Vmath::Svtvp(npts, -m_kappa, m_fields[gradphi1_idx]->GetPhys(), 1,
               out_arr[ne_idx], 1, out_arr[ne_idx], 1);
}

/**
 *  @brief Compute components of advection velocities normal to trace elements
 * (faces, in 3D).
 *
 * @param[in,out] trace_vel_norm Trace normal velocities for each field
 * @param         adv_vel        Advection velocities for each field
 */
Array<OneD, NekDouble> &HWITSystem::get_adv_vel_norm(
    Array<OneD, NekDouble> &trace_vel_norm,
    const Array<OneD, Array<OneD, NekDouble>> &adv_vel) {
  // Number of trace (interface) points
  int num_trace_pts = GetTraceNpoints();
  // Auxiliary variable to compute normal velocities
  Array<OneD, NekDouble> tmp(num_trace_pts);

  // Ensure output array is zeroed
  Vmath::Zero(num_trace_pts, trace_vel_norm, 1);

  // Compute advection vel dot trace normals and store
  for (int i = 0; i < adv_vel.size(); ++i) {
    m_fields[0]->ExtractTracePhys(adv_vel[i], tmp);
    Vmath::Vvtvp(num_trace_pts, m_traceNormals[i], 1, tmp, 1, trace_vel_norm, 1,
                 trace_vel_norm, 1);
  }
  return trace_vel_norm;
}

/**
 *  @brief Compute trace-normal advection velocities for the plasma density.
 */
Array<OneD, NekDouble> &HWITSystem::get_adv_vel_norm_elec() {
  return get_adv_vel_norm(m_norm_vel_elec, m_ExB_vel);
}

/**
 *  @brief Construct flux array.
 *
 * @param  field_vals Physical values for each advection field
 * @param  adv_vel    Advection velocities for each advection field
 * @param[out] flux       Flux array
 */
void HWITSystem::get_flux_vector(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    const Array<OneD, Array<OneD, NekDouble>> &adv_vel,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux) {
  ASSERTL1(flux[0].size() == adv_vel.size(),
           "Dimension of flux array and advection velocity array do not match");
  int npts = field_vals[0].size();

  for (auto i = 0; i < flux.size(); ++i) {
    for (auto j = 0; j < flux[0].size(); ++j) {
      Vmath::Vmul(npts, field_vals[i], 1, adv_vel[j], 1, flux[i][j], 1);
    }
  }
}

/**
 * @brief Construct the flux vector for the diffusion problem.
 * @todo not implemented
 */
void HWITSystem::get_flux_vector_diff(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscous_tensor) {
  std::cout << "*** GetFluxVectorDiff not defined! ***" << std::endl;
}

/**
 * @brief Compute the flux vector for advection in the plasma density and
 * momentum equations.
 *
 * @param field_vals Array of Fields ptrs
 * @param[out] flux Resulting flux array
 */
void HWITSystem::get_flux_vector_elec(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux) {
  get_flux_vector(field_vals, m_ExB_vel, flux);
}

/**
 * @brief Load all required session parameters into member variables.
 */
void HWITSystem::load_params() {
  TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>::load_params();

  // Type of advection to use -- in theory we also support flux reconstruction
  // for quad-based meshes, or you can use a standard convective term if you
  // were fully continuous in space. Default is DG.
  m_session->LoadSolverInfo("AdvectionType", m_adv_type, "WeakDG");

  // Magnetic field strength
  m_B = std::vector<NekDouble>(m_graph->GetSpaceDimension(), 0);
  m_session->LoadParameter("Bx", m_B[0], 0.0);
  m_session->LoadParameter("By", m_B[1], 0.0);
  m_session->LoadParameter("Bz", m_B[2], 1.0);
  // *** Bx=By=0 is assumed elsewhere in the code ***
  NESOASSERT(m_B[0] == 0, "Fluid solver doesn't yet correctly handle Bx != 0");
  NESOASSERT(m_B[1] == 0, "Fluid solver doesn't yet correctly handle By != 0");

  // Coefficient factors for potential solve
  m_session->LoadParameter("d00", m_d00, 1);
  m_session->LoadParameter("d11", m_d11, 1);
  m_session->LoadParameter("d22", m_d22, 1);

  // Type of Riemann solver to use. Default = "Upwind"
  m_session->LoadSolverInfo("UpwindType", m_riemann_solver_type, "Upwind");

  // Particle-related parameters
  m_session->LoadParameter("num_particle_steps_per_fluid_step",
                           m_num_part_substeps, 1);
  m_session->LoadParameter("particle_num_write_particle_steps",
                           m_num_write_particle_steps, 0);
  m_part_timestep = m_timestep / m_num_part_substeps;

  // Compute some properties derived from params
  m_Bmag = std::sqrt(m_B[0] * m_B[0] + m_B[1] * m_B[1] + m_B[2] * m_B[2]);
  m_b_unit = std::vector<NekDouble>(m_graph->GetSpaceDimension());
  for (auto idim = 0; idim < m_b_unit.size(); idim++) {
    m_b_unit[idim] = (m_Bmag > 0) ? m_B[idim] / m_Bmag : 0.0;
  }

  // alpha (required)
  m_session->LoadParameter("HW_alpha", m_alpha);

  // kappa (required)
  m_session->LoadParameter("HW_kappa", m_kappa);
}

/**
 * @brief Utility function to print the size of a 1D Nektar array.
 * @param arr Array to print the size of
 * @param label Label to include in the output message
 * @param all_tasks If true, print message on all tasks, else print only on
 * task 0 (default=false)
 */
void HWITSystem::print_arr_size(const Array<OneD, NekDouble> &arr,
                                std::string label, bool all_tasks) {
  if (m_session->GetComm()->TreatAsRankZero() || all_tasks) {
    if (!label.empty()) {
      std::cout << label << " ";
    }
    std::cout << "size = " << arr.size() << std::endl;
  }
}
/**
 * @brief Calls HelmSolve to solve for the electric potential
 *
 * @param in_arr Array of physical field values
 */
void HWITSystem::solve_phi(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr) {

  // Field indices
  int npts = GetNpoints();
  int w_idx = m_field_to_index.get_idx("w");
  int phi_idx = m_field_to_index.get_idx("phi");

  // Set up factors for electrostatic potential solve
  StdRegions::ConstFactorMap factors;
  // Helmholtz => Poisson (lambda = 0)
  factors[StdRegions::eFactorLambda] = 0.0;
  // Set coefficient factors
  factors[StdRegions::eFactorCoeffD00] = m_d00;
  factors[StdRegions::eFactorCoeffD11] = m_d11;
  factors[StdRegions::eFactorCoeffD22] = m_d22;

  // Solve for phi. Output of this routine is in coefficient (spectral)
  // space, so backwards transform to physical space since we'll need that
  // for the advection step & computing drift velocity.
  m_fields[phi_idx]->HelmSolve(in_arr[w_idx], m_fields[phi_idx]->UpdateCoeffs(),
                               factors);
  m_fields[phi_idx]->BwdTrans(m_fields[phi_idx]->GetCoeffs(),
                              m_fields[phi_idx]->UpdatePhys());
}

void HWITSystem::v_GenerateSummary(SU::SummaryList &s) {
  TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>::v_GenerateSummary(s);

  std::stringstream tmpss;
  tmpss << "[" << m_d00 << "," << m_d11 << "," << m_d22 << "]";
  SU::AddSummaryItem(s, "Helmsolve coeffs.", tmpss.str());

  SU::AddSummaryItem(s, "Riemann solver", m_riemann_solver_type);
  // Particle stuff
  // SU::AddSummaryItem(s, "Num. part. substeps", m_num_part_substeps);
  SU::AddSummaryItem(s, "Part. output freq", m_num_write_particle_steps);
  tmpss = std::stringstream();
  tmpss << "[" << m_B[0] << "," << m_B[1] << "," << m_B[2] << "]";
  SU::AddSummaryItem(s, "B", tmpss.str());
  SU::AddSummaryItem(s, "|B|", m_Bmag);

  SU::AddSummaryItem(s, "HW alpha", m_alpha);
  SU::AddSummaryItem(s, "HW kappa", m_kappa);
}

/**
 * @brief Post-construction class initialisation.
 *
 * @param create_field if true, create a new field object and add it to
 * m_fields. Optional, defaults to true.
 */
void HWITSystem::v_InitObject(bool create_field) {
  TimeEvoEqnSysBase::v_InitObject(create_field);

  // Since we are starting from a setup where each field is defined to be a
  // discontinuous field (and thus support DG), the first thing we do is to
  // recreate the phi field so that it is continuous, in order to support the
  // Poisson solve. Note that you can still perform a Poisson solve using a
  // discontinuous field, which is done via the hybridisable discontinuous
  // Galerkin (HDG) approach.
  int phi_idx = m_field_to_index.get_idx("phi");
  m_fields[phi_idx] = MemoryManager<MR::ContField>::AllocateSharedPtr(
      m_session, m_graph, m_session->GetVariable(phi_idx), true, true);

  // Create storage for ExB vel
  int npts = GetNpoints();
  for (int i = 0; i < m_graph->GetSpaceDimension(); ++i) {
    m_ExB_vel[i] = Array<OneD, NekDouble>(npts, 0.0);
  }

  // Type of advection class to be used. By default, we only support the
  // discontinuous projection, since this is the only approach we're
  // considering for this solver.
  ASSERTL0(m_projectionType == MR::eDiscontinuous,
           "Unsupported projection type: only discontinuous"
           " projection supported.");

  // Do not forwards transform initial condition.
  m_homoInitialFwd = false; ////

  // Define the normal velocity fields.
  // These are populated at each step (by reference) in calls to GetVnAdv()
  if (m_fields[0]->GetTrace()) {
    auto nTrace = GetTraceNpoints();
    m_norm_vel_elec = Array<OneD, NekDouble>(nTrace);
  }

  // Advection objects
  // Need one per advection velocity
  this->adv_obj =
      SU::GetAdvectionFactory().CreateInstance(m_adv_type, m_adv_type);

  // Set callback functions to compute flux vectors
  this->adv_obj->SetFluxVector(&HWITSystem::get_flux_vector_elec, this);

  // Create Riemann solvers (one per advection object) and set normal velocity
  // callback functions
  m_riemann_elec = SU::GetRiemannSolverFactory().CreateInstance(
      m_riemann_solver_type, m_session);
  m_riemann_elec->SetScalar("Vn", &HWITSystem::get_adv_vel_norm_elec, this);

  // Tell advection objects about the Riemann solvers and finish init
  this->adv_obj->SetRiemannSolver(m_riemann_elec);
  this->adv_obj->InitObject(m_session, m_fields);

  // Bind projection function for time integration object
  m_ode.DefineProjection(&HWITSystem::do_ode_projection, this);

  // Store DisContFieldSharedPtr casts of fields in a map, indexed by name
  int idx = 0;
  for (auto &field_name : m_session->GetVariables()) {
    m_discont_fields[field_name] =
        std::dynamic_pointer_cast<MR::DisContField>(m_fields[idx]);
    idx++;
  }

  if (this->particles_enabled) {
    // Set up object to evaluate density field
    this->particle_sys->setup_evaluate_grad_phi(
        m_discont_fields.at("gradphi0"), m_discont_fields.at("gradphi1"),
        m_discont_fields.at("gradphi2"));
  }

  // Bind RHS function for explicit time integration
  m_ode.DefineOdeRhs(&HWITSystem::explicit_time_int, this);
  // Bind RHS function for implicit time integration
  m_ode.DefineImplicitSolve(&HWITSystem::implicit_time_int, this);

  if (!m_explicitAdvection) {
    init_nonlin_sys_solver();
  }

  // Create diagnostic for recording growth rates
  if (m_diag_growth_rates_recording_enabled) {
    m_diag_growth_rates_recorder =
        std::make_shared<GrowthRatesRecorder<MultiRegions::DisContField>>(
            m_session, 2, m_discont_fields["ne"], m_discont_fields["w"],
            m_discont_fields["phi"], GetNpoints(), m_alpha, m_kappa);
  }
}

void HWITSystem::init_nonlin_sys_solver() {
  int npts_tot = 2 * m_fields[0]->GetNpoints();

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
  m_session->LoadParameter("NewtonRelativeIteTol", key.m_NekNonLinSysTolerance,
                           1.0E-12);
  WARNINGL0(!m_session->DefinesParameter("NewtonAbsoluteIteTol"),
            "Please specify NewtonRelativeIteTol instead of "
            "NewtonAbsoluteIteTol in XML session file");
  m_session->LoadParameter("NonlinIterTolRelativeL2",
                           key.m_NonlinIterTolRelativeL2, 1.0E-3);
  m_session->LoadSolverInfo("LinSysIterSolverTypeInNonlin",
                            key.m_LinSysIterSolverTypeInNonlin, "GMRES");

  LibUtilities::NekSysOperators nekSysOp;
  nekSysOp.DefineNekSysResEval(&HWITSystem::nonlin_sys_evaluator_1D, this);
  nekSysOp.DefineNekSysLhsEval(&HWITSystem::matrix_multiply_matrix_free, this);
  nekSysOp.DefineNekSysPrecon(&HWITSystem::do_null_precon, this);

  // Initialize non-linear system
  this->nonlin_sys = LibUtilities::GetNekNonlinSysIterFactory().CreateInstance(
      "Newton", m_session, m_comm->GetRowComm(), npts_tot, key);
  this->nonlin_sys->SetSysOperators(nekSysOp);
}

void HWITSystem::implicit_time_int(
    const Array<OneD, const Array<OneD, NekDouble>> &inpnts,
    Array<OneD, Array<OneD, NekDouble>> &outpnt, const NekDouble time,
    const NekDouble lambda) {
  this->time_int_lambda = lambda;
  this->bnd_evaluate_time = time;
  unsigned int npoints = m_fields[0]->GetNpoints();

  Array<OneD, NekDouble> in_arr(2 * npoints);
  Array<OneD, NekDouble> outarray(2 * npoints, 0.0);
  Array<OneD, NekDouble> tmp;

  for (int i = 0; i < 2; ++i) {
    int noffset = i * npoints;
    Vmath::Vcopy(npoints, inpnts[i], 1, tmp = in_arr + noffset, 1);
  }

  implicit_time_int_1D(in_arr, outarray);

  for (int i = 0; i < 2; ++i) {
    int noffset = i * npoints;
    Vmath::Vcopy(npoints, outarray + noffset, 1, outpnt[i], 1);
  }
}

void HWITSystem::implicit_time_int_1D(
    const Array<OneD, const NekDouble> &in_arr, Array<OneD, NekDouble> &out) {
  calc_ref_vals(in_arr);

  this->nonlin_sys->SetRhsMagnitude(this->in_arr_norm);

  this->newton_its_counter +=
      this->nonlin_sys->SolveSystem(in_arr.size(), in_arr, out, 0);

  this->lin_its_counter += this->nonlin_sys->GetNtotLinSysIts();

  this->imp_stages_counter++;
}

void HWITSystem::calc_ref_vals(const Array<OneD, const NekDouble> &in_arr) {
  unsigned int npoints = m_fields[0]->GetNpoints();

  Array<OneD, NekDouble> magnitdEstimat(2, 0.0);

  for (int i = 0; i < 2; ++i) {
    int offset = i * npoints;
    magnitdEstimat[i] = Vmath::Dot(npoints, in_arr + offset, in_arr + offset);
  }
  m_comm->GetSpaceComm()->AllReduce(magnitdEstimat,
                                    Nektar::LibUtilities::ReduceSum);

  this->in_arr_norm = 0.0;
  for (int i = 0; i < 2; ++i) {
    this->in_arr_norm += magnitdEstimat[i];
  }
}

void HWITSystem::nonlin_sys_evaluator_1D(
    const Array<OneD, const NekDouble> &in_arr, Array<OneD, NekDouble> &out,
    [[maybe_unused]] const bool &flag) {
  unsigned int npoints = m_fields[0]->GetNpoints();
  Array<OneD, Array<OneD, NekDouble>> in2D(2);
  Array<OneD, Array<OneD, NekDouble>> out2D(2);
  for (int i = 0; i < 2; ++i) {
    int offset = i * npoints;
    in2D[i] = in_arr + offset;
    out2D[i] = out + offset;
  }
  nonlin_sys_evaluator(in2D, out2D);
}

void HWITSystem::nonlin_sys_evaluator(
    const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &out) {
  unsigned int npoints = m_fields[0]->GetNpoints();
  Array<OneD, Array<OneD, NekDouble>> inpnts(2);
  for (int i = 0; i < 2; ++i) {
    inpnts[i] = Array<OneD, NekDouble>(npoints, 0.0);
  }

  do_ode_projection(in_arr, inpnts, this->bnd_evaluate_time);
  explicit_time_int(inpnts, out, this->bnd_evaluate_time);

  for (int i = 0; i < 2; ++i) {
    Vmath::Svtvp(npoints, -this->time_int_lambda, out[i], 1, in_arr[i], 1,
                 out[i], 1);
    Vmath::Vsub(npoints, out[i], 1,
                this->nonlin_sys->GetRefSourceVec() + i * npoints, 1, out[i],
                1);
  }
}

void HWITSystem::matrix_multiply_matrix_free(
    const Array<OneD, const NekDouble> &in_arr, Array<OneD, NekDouble> &out,
    [[maybe_unused]] const bool &flag) {
  const Array<OneD, const NekDouble> solref =
      this->nonlin_sys->GetRefSolution();
  const Array<OneD, const NekDouble> resref =
      this->nonlin_sys->GetRefResidual();

  unsigned int ntotal = in_arr.size();
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

void HWITSystem::do_null_precon(const Array<OneD, const NekDouble> &in_arr,
                                Array<OneD, NekDouble> &outarray,
                                [[maybe_unused]] const bool &flag) {
  Vmath::Vcopy(in_arr.size(), in_arr, 1, outarray, 1);
}

/**
 * @brief Override v_PostIntegrate to do particle output and compute
 * diagnostics, if enabled
 * @param step Time step number
 */
bool HWITSystem::v_PostIntegrate(int step) {
  if (m_diag_growth_rates_recording_enabled) {
    m_diag_growth_rates_recorder->compute(step);
  }

  m_solver_callback_handler.call_post_integrate(this);

  // Writes a step of the particle trajectory.
  if (this->particles_enabled && m_num_write_particle_steps > 0 &&
      (step % m_num_write_particle_steps) == 0) {
    this->particle_sys->write(step);
  }
  return TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>::v_PostIntegrate(
      step);
}

/**
 * @brief Override v_PreIntegrate to do particle system integration, and initial
 * set up for mass recording diagnostic (first call only)
 *
 * @param step Time step number
 */
bool HWITSystem::v_PreIntegrate(int step) {
  m_solver_callback_handler.call_pre_integrate(this);

  if (this->particles_enabled) {
    // Integrate the particle system to the requested time.
    this->particle_sys->evaluate_fields();
    this->particle_sys->integrate(m_time + m_timestep, m_part_timestep);
  }

  return UnsteadySystem::v_PreIntegrate(step);
}

/**
 * @brief After reading ICs, calculate phi and grad(phi)
 */
void HWITSystem::v_SetInitialConditions(NekDouble init_time, bool dump_ICs,
                                        const int domain) {
  TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>::v_SetInitialConditions(
      init_time, dump_ICs, domain);
  calc_init_phi_and_gradphi();
  if(this->particle_sys){
    this->particle_sys->initialise_particles_from_fields();
  }
}

/**
 * @brief Convenience function to zero a Nektar Array of 1D Arrays.
 *
 * @param out_arr Array of 1D arrays to be zeroed
 *
 */
void HWITSystem::zero_out_array(Array<OneD, Array<OneD, NekDouble>> &out_arr) {
  for (auto ifld = 0; ifld < out_arr.size(); ifld++) {
    Vmath::Zero(out_arr[ifld].size(), out_arr[ifld], 1);
  }
}
} // namespace NESO::Solvers::hw_impurity_transport
