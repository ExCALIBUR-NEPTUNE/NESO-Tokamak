#include "ReducedBraginskii.hpp"
#include "../RiemannSolvers/TokamakSolver.hpp"
#include <SolverUtils/Advection/AdvectionNonConservative.h>
#include <SolverUtils/Advection/AdvectionWeakDG.h>

namespace NESO::Solvers::tokamak
{

/// Name of class
static std::string class_name;
std::string ReducedBraginskii::class_name =
    SU::GetEquationSystemFactory().RegisterCreatorFunction(
        "ReducedBraginskii", ReducedBraginskii::create,
        "Solves reduced Braginskii equations");
/**
 * @brief Creates an instance of this class.
 */
static SU::EquationSystemSharedPtr create(
    const LU::SessionReaderSharedPtr &session,
    const SD::MeshGraphSharedPtr &graph)
{
    SU::EquationSystemSharedPtr p =
        MemoryManager<ReducedBraginskii>::AllocateSharedPtr(session, graph);
    p->InitObject();
    return p;
}

ReducedBraginskii::ReducedBraginskii(const LU::SessionReaderSharedPtr &session,
                                     const SD::MeshGraphSharedPtr &graph)
    : TokamakSystem(session, graph)
{
    this->n_indep_fields       = 1; // p_e
    this->n_fields_per_species = 3; // n_i, v_i, p_i
}

void ReducedBraginskii::v_InitObject(bool DeclareFields)
{
    TokamakSystem::v_InitObject(DeclareFields);
    m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
        m_session, m_spacedim, m_graph);

    std::string diffName;
    m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
    m_diffusion =
        SolverUtils::GetDiffusionFactory().CreateInstance(diffName, diffName);
    m_diffusion->SetFluxVector(&ReducedBraginskii::GetFluxVectorDiff, this);

    // workaround for bug in DiffusionLDG
    m_difffields = Array<OneD, MR::ExpListSharedPtr>(m_indfields.size());
    for (int f = 0; f < m_difffields.size(); ++f)
    {
        m_difffields[f] = m_indfields[f];
    }

    m_diffusion->InitObject(m_session, m_difffields);

    // Create storage for velocities
    int npts = GetNpoints();
    pe_idx   = 3 * this->n_species;

    // Parallel velocities
    this->v_e_par = Array<OneD, NekDouble>(npts, 0.0);
    this->v_i_par = std::vector<Array<OneD, NekDouble>>(n_species);

    // Per-field advection velocities
    this->adv_vel = Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(
        m_indfields.size());

    for (int i = 0; i < this->adv_vel.size(); ++i)
    {
        this->adv_vel[i] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    }

    for (int d = 0; d < m_spacedim; ++d)
    {
        this->adv_vel[pe_idx][d] = Array<OneD, NekDouble>(npts, 0.0);
    }

    int s = 0;
    for (const auto &[k, v] : this->neso_config->get_species())
    {
        ni_idx.push_back(3 * s);
        vi_idx.push_back(1 + 3 * s);
        pi_idx.push_back(2 + 3 * s);

        this->v_i_par[s] = Array<OneD, NekDouble>(npts, 0.0);
        for (int d = 0; d < m_graph->GetSpaceDimension(); ++d)
        {
            this->adv_vel[ni_idx[s]][d] = Array<OneD, NekDouble>(npts, 0.0);
            this->adv_vel[vi_idx[s]][d] = Array<OneD, NekDouble>(npts, 0.0);
            this->adv_vel[pi_idx[s]][d] = Array<OneD, NekDouble>(npts, 0.0);
        }

        double charge, mass;
        this->neso_config->load_species_parameter(k, "Charge", charge);
        this->neso_config->load_species_parameter(k, "Mass", mass);
        m_varConv->charge[s] = charge;
        m_varConv->mass[s]   = mass;
        s++;
    }

    m_varConv->ni_idx = ni_idx;
    m_varConv->vi_idx = vi_idx;
    m_varConv->pi_idx = pi_idx;
    m_varConv->pe_idx = pe_idx;

    for (auto &[k, bnd] : m_bndConds->GetBounds())
    {
        bnd->pe_idx      = this->pe_idx;
        bnd->ni_idx      = this->ni_idx;
        bnd->vi_idx      = this->vi_idx;
        bnd->pi_idx      = this->pi_idx;
        bnd->neso_config = this->neso_config;
    }

    if (m_indfields[0]->GetTrace())
    {
        auto nTrace = GetTraceNpoints();

        this->trace_vel_norm =
            Array<OneD, Array<OneD, NekDouble>>(adv_vel.size());
        this->trace_b_norm = Array<OneD, NekDouble>(nTrace, 0.0);
        this->adv_vel_trace =
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(adv_vel.size());
        for (int i = 0; i < this->trace_vel_norm.size(); ++i)
        {
            this->adv_vel_trace[i] =
                Array<OneD, Array<OneD, NekDouble>>(m_spacedim);

            this->trace_vel_norm[i] = Array<OneD, NekDouble>(nTrace, 0.0);
            for (int d = 0; d < m_spacedim; ++d)
            {
                this->adv_vel_trace[i][d] = Array<OneD, NekDouble>(nTrace, 0.0);
            }
        }
    }

    // Create Riemann solver and set normal velocity
    // callback functions
    this->riemann_solver = SU::GetRiemannSolverFactory().CreateInstance(
        this->riemann_solver_type, m_session);
    auto t = std::dynamic_pointer_cast<TokamakSolver>(this->riemann_solver);
    // t->ni_idx    = ni_idx;
    // t->vi_idx    = vi_idx;
    // t->pi_idx    = pi_idx;
    // t->pe_idx    = pe_idx;

    this->riemann_solver->SetVector("Vn", &ReducedBraginskii::GetAdvVelNorm,
                                    this);

    // Setup advection object
    m_session->LoadSolverInfo("AdvectionType", this->adv_type, "WeakDG");
    m_advection = SU::GetAdvectionFactory().CreateInstance(this->adv_type,
                                                           this->adv_type);
    m_advection->SetFluxVector(&ReducedBraginskii::GetFluxVector, this);
    m_advection->SetRiemannSolver(this->riemann_solver);
    m_advection->InitObject(m_session, m_indfields);

    m_ode.DefineOdeRhs(&ReducedBraginskii::DoOdeRhs, this);

    if (this->particles_enabled)
    {
        std::vector<Sym<REAL>> src_syms;
        std::vector<int> src_components;

        s = 0;
        for (const auto &[k, v] : this->particle_sys->get_species())
        {
            ni_src_idx.push_back(s * (2 + m_spacedim));
            pi_src_idx.push_back(1 + s * (2 + m_spacedim));
            vi_src_idx.push_back(2 + s * (2 + m_spacedim));

            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_DENSITY"));
            src_components.push_back(0);

            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));

            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_ENERGY"));
            src_components.push_back(0);

            for (int d = 0; d < this->m_spacedim; ++d)
            {
                this->src_fields.emplace_back(
                    MemoryManager<MR::DisContField>::AllocateSharedPtr(
                        *std::dynamic_pointer_cast<MR::DisContField>(
                            m_fields[0])));
                src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_MOMENTUM"));
                src_components.push_back(d);
            }
            s++;
        }
        this->src_fields.emplace_back(
            MemoryManager<MR::DisContField>::AllocateSharedPtr(
                *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
        src_syms.push_back(Sym<REAL>("ELECTRON_SOURCE_ENERGY"));
        src_components.push_back(0);

        this->particle_sys->finish_setup(this->src_fields, src_syms,
                                         src_components);
    }
}

bool ReducedBraginskii::v_PostIntegrate(int step)
{
    m_fields[0]->FwdTrans(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs());
    m_fields[1]->FwdTrans(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs());

    // Writes a step of the particle trajectory.

    return TokamakSystem::v_PostIntegrate(step);
}

/**
 * @brief Populate rhs array ( @p outarray )
 *
 * @param inarray physical values of all fields
 * @param[out] outarray output array (RHSs of time integration equations)
 */
void ReducedBraginskii::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Get field indices
    int npts      = GetNpoints();
    int nTracePts = GetTraceTotPoints();
    for (int f = 0; f < outarray.size(); ++f)
    {
        Vmath::Zero(npts, outarray[f], 1);
    }

    int nvariables = inarray.size();

    m_varConv->GetElectronDensity(inarray, m_fields[0]->UpdatePhys());

    // Store forwards/backwards space along trace space
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        m_indfields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
    }

    // Calculate E
    ComputeE();
    // Calculate ExB, parallel and diamagnetic velocities
    CalcVelocities(inarray, m_fields[0]->GetPhys(), this->adv_vel);

    // Perform advection
    DoAdvection(inarray, outarray, time, Fwd, Bwd);

    m_bndConds->Update(inarray, time);

    // DoExtra(inarray, outarray);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Neg(npts, outarray[i], 1);
    }

    CalcKappaTensor();

    // Perform Diffusion
    DoDiffusion(inarray, outarray, Fwd, Bwd);

    if (this->particles_enabled)
    {
        DoParticles(inarray, outarray);
    }

    // Add forcing terms
    for (auto &x : m_forcing)
    {
        x->Apply(m_fields, inarray, outarray, time);
    }
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void ReducedBraginskii::DoAdvection(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = inarray.size();
    Array<OneD, Array<OneD, NekDouble>> advVel(m_spacedim);

    m_advection->Advect(nvariables, m_indfields, advVel, inarray, outarray,
                        time, pFwd, pBwd);
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void ReducedBraginskii::DoParticles(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int npts = GetNpoints();

    // Add contribution to electron energy
    Vmath::Vadd(npts, outarray[pe_idx], 1, this->src_fields[0]->GetPhys(), 1,
                outarray[pe_idx], 1);

    int s = 0;
    for (const auto &[k, v] : this->neso_config->get_species())
    {
        double mass;
        this->neso_config->load_species_parameter(k, "Mass", mass);
        //  Add contribution to ion density
        Vmath::Vadd(npts, outarray[ni_idx[s]], 1,
                    this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                    outarray[ni_idx[s]], 1);

        // Add contribution to ion energy
        Vmath::Vadd(npts, outarray[pi_idx[s]], 1,
                    this->src_fields[pi_src_idx[s]]->GetPhys(), 1,
                    outarray[pi_idx[s]], 1);
        // Add number density source contribution to ion energy
        Array<OneD, NekDouble> dynamic_energy(npts);
        m_varConv->GetIonDynamicEnergy(s, mass, inarray, dynamic_energy);
        Vmath::Vvtvp(npts, dynamic_energy, 1,
                     this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                     outarray[pi_idx[s]], 1, outarray[pi_idx[s]], 1);

        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(npts, this->b_unit[d], 1,
                         this->src_fields[vi_src_idx[s] + d]->GetPhys(), 1,
                         outarray[vi_idx[s]], 1, outarray[vi_idx[s]], 1);
        }
        s++;
    }
}

/**
 * @brief Compute the gradient of phi for evaluation at the particle positions.
 */
void ReducedBraginskii::ComputeE()
{
    int npts = GetNpoints();

    Vmath::Neg(npts, this->E[0]->UpdatePhys(), 1);
    Vmath::Neg(npts, this->E[1]->UpdatePhys(), 1);
    Vmath::Neg(npts, this->E[2]->UpdatePhys(), 1);

    this->E[0]->FwdTrans(this->E[0]->GetPhys(), this->E[0]->UpdateCoeffs());
    this->E[1]->FwdTrans(this->E[1]->GetPhys(), this->E[1]->UpdateCoeffs());
    this->E[2]->FwdTrans(this->E[2]->GetPhys(), this->E[2]->UpdateCoeffs());
}

void ReducedBraginskii::CalcVelocities(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &adv_vel)
{
    int npts = inarray[0].size();
    for (int f = 0; f < adv_vel.size(); ++f)
    {
        for (int d = 0; d < adv_vel[f].size(); ++d)
        {
            Vmath::Zero(npts, adv_vel[f][d], 1);
        }
    }

    // Zero Electron velocity
    Vmath::Zero(npts, m_fields[1]->UpdatePhys(), 1);
    int s = 0;
    for (const auto &[k, v] : this->neso_config->get_species())
    {
        double charge = 1;
        double mass   = 1;
        this->neso_config->load_species_parameter(k, "Charge", charge);
        this->neso_config->load_species_parameter(k, "Mass", mass);

        // Calculate Ion parallel velocities
        Vmath::Smul(npts, 1.0 / mass, inarray[vi_idx[s]], 1, this->v_i_par[s],
                    1);
        Vmath::Svtvp(npts, charge, this->v_i_par[s], 1, m_fields[1]->GetPhys(),
                     1, m_fields[1]->UpdatePhys(), 1);

        Vmath::Vdiv(npts, this->v_i_par[s], 1, inarray[ni_idx[s]], 1,
                    this->v_i_par[s], 1);

        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vmul(npts, this->b_unit[d], 1, this->v_i_par[s], 1,
                        adv_vel[ni_idx[s]][d], 1);
            Vmath::Vcopy(npts, adv_vel[ni_idx[s]][d], 1, adv_vel[vi_idx[s]][d],
                         1);
            Vmath::Vcopy(npts, adv_vel[ni_idx[s]][d], 1, adv_vel[pi_idx[s]][d],
                         1);
        }
        s++;
    }

    Vmath::Vdiv(npts, m_fields[1]->GetPhys(), 1, ne, 1, this->v_e_par, 1);

    for (int d = 0; d < m_spacedim; ++d)
    {
        Vmath::Vmul(npts, this->b_unit[d], 1, this->v_e_par, 1,
                    adv_vel[pe_idx][d], 1);
    }
}

/**
 *  @brief Compute components of advection velocities normal to trace elements
 * (faces, in 3D).
 *
 * @param[in,out] trace_vel_norm Trace normal velocities for each field
 * @param         adv_vel_trace        Advection velocities for each field
 */
Array<OneD, Array<OneD, NekDouble>> &ReducedBraginskii::GetAdvVelNorm()
{
    // Number of trace (interface) points
    int num_trace_pts = GetTraceNpoints();
    // Auxiliary variable to compute normal velocities

    // Compute advection vel dot trace normals and store
    for (int j = 0; j < this->adv_vel_trace.size(); ++j)
    {
        Array<OneD, Array<OneD, NekDouble>> normals(m_spacedim);

        for (int d = 0; d < m_spacedim; ++d)
        {
            normals[d] = Array<OneD, NekDouble>(num_trace_pts);
        }
        m_indfields[j]->GetTrace()->GetNormals(normals);
        // Ensure output array is zeroed
        Vmath::Zero(num_trace_pts, this->trace_vel_norm[j], 1);
        for (int d = 0; d < this->adv_vel_trace[j].size(); ++d)
        {
            m_indfields[j]->ExtractTracePhys(this->adv_vel[j][d],
                                             this->adv_vel_trace[j][d]);
            Vmath::Vvtvp(num_trace_pts, normals[d], 1,
                         this->adv_vel_trace[j][d], 1, this->trace_vel_norm[j],
                         1, this->trace_vel_norm[j], 1);
        }
    }
    return this->trace_vel_norm;
}

/**
 *  @brief Construct flux array.
 *
 * @param  field_vals Physical values for each advection field
 * @param[out] flux       Flux array
 */
void ReducedBraginskii::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    int npts       = field_vals[0].size();
    int nVariables = field_vals.size();

    for (int d = 0; d < m_spacedim; ++d)
    {
        // Electron Energy Flux
        Vmath::Vmul(npts, field_vals[pe_idx], 1, this->adv_vel[pe_idx][d], 1,
                    fluxes[pe_idx][d], 1);
    }

    int s = 0;
    for (const auto &[k, v] : this->neso_config->get_species())
    {
        double charge, mass;
        this->neso_config->load_species_parameter(k, "Charge", charge);
        this->neso_config->load_species_parameter(k, "Mass", mass);

        for (int d = 0; d < m_spacedim; ++d)
        {

            // Ion Density Flux
            Vmath::Vmul(npts, field_vals[ni_idx[s]], 1,
                        this->adv_vel[ni_idx[s]][d], 1, fluxes[ni_idx[s]][d],
                        1);
            // Ion Momentum Flux
            Vmath::Vmul(npts, field_vals[vi_idx[s]], 1,
                        this->adv_vel[vi_idx[s]][d], 1, fluxes[vi_idx[s]][d],
                        1);
            // Ion Energy Flux
            Vmath::Vmul(npts, field_vals[pi_idx[s]], 1,
                        this->adv_vel[pi_idx[s]][d], 1, fluxes[pi_idx[s]][d],
                        1);
        }

        s++;
    }
}

void ReducedBraginskii::DoDiffusion(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = inarray.size();
    int npointsIn  = GetNpoints();
    int npointsOut = npointsIn;
    int nTracePts  = GetTraceTotPoints();

    // this should be preallocated
    Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>(npointsOut, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> inarrayDiff(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inFwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inBwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        inarrayDiff[i] = Array<OneD, NekDouble>(npointsIn, 0.0);
        inFwd[i]       = Array<OneD, NekDouble>(nTracePts, 0.0);
        inBwd[i]       = Array<OneD, NekDouble>(nTracePts, 0.0);
    }

    // Extract temperature
    m_varConv->GetElectronTemperature(inarray, inarrayDiff[pe_idx]);

    int s = 0;
    double mass;

    for (const auto &[k, v] : this->neso_config->get_species())
    {
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

            m_varConv->GetIonTemperature(s, mass, pFwd, inFwd[pi_idx[s]]);
            m_varConv->GetIonTemperature(s, mass, pBwd, inBwd[pi_idx[s]]);

            s++;
        }
    }

    m_diffusion->Diffuse(nvariables, m_difffields, inarrayDiff, outarrayDiff,
                         inFwd, inBwd);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(npointsOut, outarrayDiff[i], 1, outarray[i], 1, outarray[i],
                    1);
    }
}

void ReducedBraginskii::CalcKPar()
{
    // Change to fn of fields
    int npoints = m_fields[0]->GetNpoints();
    NekDouble k_par;
    m_session->LoadParameter("k_par", k_par, 100.0);
    m_kpar = Array<OneD, NekDouble>(npoints, k_par);
}

void ReducedBraginskii::CalcKPerp()
{
    // Change to fn of fields
    int npoints = m_fields[0]->GetNpoints();
    NekDouble k_perp;
    m_session->LoadParameter("k_perp", k_perp, 1.0);
    m_kperp = Array<OneD, NekDouble>(npoints, k_perp);
}

void ReducedBraginskii::CalcDiffTensor()
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

void ReducedBraginskii::CalcKappaPar()
{
    // Change to fn of T
    int npoints = m_fields[0]->GetNpoints();
    NekDouble kappa_par;
    m_session->LoadParameter("kappa_par", kappa_par, 100.0);
    m_kappapar = Array<OneD, NekDouble>(npoints, kappa_par);
}

void ReducedBraginskii::CalcKappaPerp()
{
    // Change to fn of T
    int npoints = m_fields[0]->GetNpoints();
    NekDouble kappa_perp;
    m_session->LoadParameter("kappa_perp", kappa_perp, 1.0);
    m_kappaperp = Array<OneD, NekDouble>(npoints, kappa_perp);
}

void ReducedBraginskii::CalcKappaTensor()
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
 * @brief Construct the flux vector for the anisotropic diffusion problem.
 */
void ReducedBraginskii::GetFluxVectorDiff(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    unsigned int nPts = inarray[0].size();

    for (int j = 0; j < m_spacedim; ++j)
    {
        // Calc diffusion of n with D tensor
        Vmath::Vmul(nPts, m_kappa[vc[j][0]].GetValue(), 1, qfield[0][pe_idx], 1,
                    fluxes[j][pe_idx], 1);
        for (int k = 1; k < m_spacedim; ++k)
        {
            Vmath::Vvtvp(nPts, m_kappa[vc[j][k]].GetValue(), 1,
                         qfield[k][pe_idx], 1, fluxes[j][pe_idx], 1,
                         fluxes[j][pe_idx], 1);
        }
    }
    int s = 0;
    for (const auto &[k, v] : this->neso_config->get_species())
    {
        for (int j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vmul(nPts, m_kappa[vc[j][0]].GetValue(), 1,
                        qfield[0][pi_idx[s]], 1, fluxes[j][pi_idx[s]], 1);
            for (int k = 1; k < m_spacedim; ++k)
            {
                Vmath::Vvtvp(nPts, m_kappa[vc[j][k]].GetValue(), 1,
                             qfield[k][pi_idx[s]], 1, fluxes[j][pi_idx[s]], 1,
                             fluxes[j][pi_idx[s]], 1);
            }
        }
        s++;
    }
}

/**
 * @brief Populate rhs array ( @p outarray )
 *
 * @param inarray physical values of all fields
 * @param[out] outarray output array (RHSs of time integration equations)
 */
void ReducedBraginskii::DoOdeImplicitRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    int nvariables = inarray.size();
    int ncoeffs    = m_fields[0]->GetNcoeffs();

    Array<OneD, Array<OneD, NekDouble>> tmpOut(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        tmpOut[i] = Array<OneD, NekDouble>(ncoeffs);
    }

    DoOdeRhsCoeff(inarray, tmpOut, time);

    for (int i = 0; i < nvariables; ++i)
    {
        m_fields[i]->BwdTrans(tmpOut[i], outarray[i]);
    }
}

/**
 * @brief Compute the right-hand side.
 */
void ReducedBraginskii::DoOdeRhsCoeff(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{

    int nvariables = inarray.size();
    int nTracePts  = GetTraceTotPoints();
    int ncoeffs    = GetNcoeffs();

    // Store forwards/backwards space along trace space
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        m_indfields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
    }

    // Calculate advection
    DoAdvectionCoeff(inarray, outarray, time, Fwd, Bwd);

    // Negate results
    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Neg(ncoeffs, outarray[i], 1);
    }

    // Add diffusion terms

    DoDiffusionCoeff(inarray, outarray, Fwd, Bwd);

    if (this->particles_enabled)
    {
        DoParticlesCoeff(inarray, outarray);
    }

    // Add forcing terms
    for (auto &x : m_forcing)
    {
        x->ApplyCoeff(m_indfields, inarray, outarray, time);
    }
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void ReducedBraginskii::DoAdvectionCoeff(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = inarray.size();
    Array<OneD, Array<OneD, NekDouble>> advVel(m_spacedim);
    std::dynamic_pointer_cast<SU::AdvectionWeakDG>(m_advection)
        ->AdvectCoeffs(nvariables, m_indfields, advVel, inarray, outarray, time,
                       pFwd, pBwd);
}

void ReducedBraginskii::DoDiffusionCoeff(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, const Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, const Array<OneD, NekDouble>> &pBwd)
{
    size_t nvariables = inarray.size();
    size_t npoints    = GetNpoints();
    size_t ncoeffs    = GetNcoeffs();
    size_t nTracePts  = GetTraceTotPoints();

    Array<OneD, Array<OneD, NekDouble>> outarrayDiff{nvariables};
    for (int i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>{ncoeffs, 0.0};
    }

    // if (m_is_diffIP)
    // {
    //     m_diffusion->DiffuseCoeffs(nvariables, m_fields, inarray,
    //     outarrayDiff,
    //                                m_bndEvaluateTime, pFwd, pBwd);
    //     for (int i = 0; i < nvariables; ++i)
    //     {
    //         Vmath::Vadd(ncoeffs, outarrayDiff[i], 1, outarray[i], 1,
    //                     outarray[i], 1);
    //     }
    // }
    // else
    // {
    ASSERTL1(false, "LDGNS not yet validated for implicit compressible "
                    "flow solver");
    Array<OneD, Array<OneD, NekDouble>> inarrayDiff{nvariables};
    Array<OneD, Array<OneD, NekDouble>> inFwd{nvariables};
    Array<OneD, Array<OneD, NekDouble>> inBwd{nvariables};

    for (int i = 0; i < nvariables; ++i)
    {
        inarrayDiff[i] = Array<OneD, NekDouble>{npoints};
        inFwd[i]       = Array<OneD, NekDouble>{nTracePts};
        inBwd[i]       = Array<OneD, NekDouble>{nTracePts};
    }

    // Extract temperature
    m_varConv->GetElectronTemperature(inarray, inarrayDiff[pe_idx]);
    int s = 0;
    for (const auto &[k, v] : this->particle_sys->get_species())
    {
        m_varConv->GetIonTemperature(s, v.mass, inarray,
                                     inarrayDiff[pi_idx[s]]);
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
        s = 0;
        for (const auto &[k, v] : this->particle_sys->get_species())
        {
            m_varConv->GetIonTemperature(s, v.mass, pFwd, inFwd[pi_idx[s]]);
            m_varConv->GetIonTemperature(s, v.mass, pBwd, inBwd[pi_idx[s]]);
            s++;
        }
    }

    // Diffusion term in coeff rhs form
    m_diffusion->DiffuseCoeffs(nvariables, m_indfields, inarrayDiff,
                               outarrayDiff, inFwd, inBwd);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(ncoeffs, outarrayDiff[i], 1, outarray[i], 1, outarray[i],
                    1);
    }
    //}
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void ReducedBraginskii::DoParticlesCoeff(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int npts   = GetNpoints();
    int ncoeff = GetNcoeffs();
    Array<OneD, NekDouble> tmp(ncoeff, 0.0);

    // Add contribution to electron energy
    m_indfields[pe_idx]->FwdTrans(this->src_fields[0]->GetPhys(), tmp);
    Vmath::Vadd(npts, outarray[pe_idx], 1, tmp, 1, outarray[pe_idx], 1);

    int s = 0;
    for (const auto &[k, v] : this->particle_sys->get_species())
    {
        //  Add contribution to ion density
        m_indfields[ni_idx[s]]->FwdTrans(
            this->src_fields[ni_src_idx[s]]->GetPhys(), tmp);
        Vmath::Vadd(npts, outarray[ni_idx[s]], 1, tmp, 1, outarray[ni_idx[s]],
                    1);

        // Add contribution to ion energy
        m_indfields[pi_idx[s]]->FwdTrans(
            this->src_fields[pi_src_idx[s]]->GetPhys(), tmp);
        Vmath::Vadd(npts, outarray[pi_idx[s]], 1, tmp, 1, outarray[pi_idx[s]],
                    1);

        // Add number density source contribution to ion energy
        Array<OneD, NekDouble> dynamic_energy(npts);
        m_varConv->GetIonDynamicEnergy(s, v.mass, inarray, dynamic_energy);
        Vmath::Vmul(npts, dynamic_energy, 1,
                    this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                    dynamic_energy, 1);
        m_fields[pi_idx[s]]->FwdTrans(dynamic_energy, tmp);
        Vmath::Vadd(npts, outarray[pi_idx[s]], 1, tmp, 1, outarray[pi_idx[s]],
                    1);

        Vmath::Zero(npts, dynamic_energy, 1);
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(npts, this->b_unit[d], 1,
                         this->src_fields[vi_src_idx[s] + d]->GetPhys(), 1,
                         dynamic_energy, 1, dynamic_energy, 1);
        }
        m_fields[vi_idx[s]]->FwdTrans(dynamic_energy, tmp);
        Vmath::Vadd(npts, outarray[vi_idx[s]], 1, tmp, 1, outarray[vi_idx[s]],
                    1);
        s++;
    }
}

/**
 * @brief After reading ICs, calculate phi and grad(phi)
 */
void ReducedBraginskii::v_SetInitialConditions(NekDouble init_time,
                                               bool dump_ICs, const int domain)
{
    TokamakSystem::v_SetInitialConditions(init_time, dump_ICs, domain);
    Checkpoint_Output(0);
}

// void ReducedBraginskii::v_SetBoundaryConditions(
//     Array<OneD, Array<OneD, NekDouble>> &physarray, NekDouble time)
// {
//     size_t nTracePts  = GetTraceTotPoints();
//     size_t nvariables = physarray.size() - 1;
//     Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
//     for (size_t i = 0; i < nvariables; ++i)
//     {
//         Fwd[i] = Array<OneD, NekDouble>(nTracePts);
//         int p  = i < n_species * n_fields_per_species
//                      ? i % n_fields_per_species
//                      : i - n_species * n_fields_per_species;
//         m_fields[p]->ExtractTracePhys(physarray[i], Fwd[i]);
//     }
//     for (auto &bc : m_bndConds)
//     {
//         bc->Apply(Fwd, physarray, time);
//     }
// }

void ReducedBraginskii::v_ExtraFldOutput(
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
            i += (2 + m_spacedim);
        }
    }
}
} // namespace NESO::Solvers::tokamak