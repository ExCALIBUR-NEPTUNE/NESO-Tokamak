#include "ElectrostaticTurbulence.hpp"
#include "../RiemannSolvers/PlasmaSolver.hpp"
#include <SolverUtils/Advection/AdvectionNonConservative.h>
#include <SolverUtils/Advection/AdvectionWeakDG.h>

namespace PENKNIFE
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
    : PlasmaSystem(session, graph)
{
    this->n_indep_fields = 3; // p_e, w, phi
}

/**
 * @brief Initialise the class.
 */
void ElectrostaticTurbulence::v_InitObject(bool DeclareFields)
{
    PlasmaSystem::v_InitObject(DeclareFields);
    ee_idx    = m_indfields.size() - this->n_indep_fields;
    omega_idx = m_indfields.size() - this->n_indep_fields + 1;
    phi_idx   = m_indfields.size() - this->n_indep_fields + 2;

    m_varConv->ee_idx = this->ee_idx;

    m_temps = Array<OneD, MR::ExpListSharedPtr>(n_species + 1);
    int f   = 0;
    for (const auto &[s, v] : this->GetSpecies())
    {
        if (v.fields.find(field_to_index["e"]) != v.fields.end())
        {
            int ei_idx = v.fields.at(field_to_index["e"]);

            m_temps[f++] = MemoryManager<MR::DisContField>::AllocateSharedPtr(
                *std::dynamic_pointer_cast<MR::DisContField>(
                    m_indfields[ei_idx]),
                m_graph, "e");
        }
    }
    m_temps[f] = MemoryManager<MR::DisContField>::AllocateSharedPtr(
        *std::dynamic_pointer_cast<MR::DisContField>(m_indfields[ee_idx]),
        m_graph, "e");

    std::string diffName;
    m_session->LoadSolverInfo("DiffusionType", diffName, "LDGET");
    // m_diffusion =
    //     SolverUtils::GetDiffusionFactory().CreateInstance(diffName,
    //     diffName);
    // m_diffusion->SetFluxVector(&ElectrostaticTurbulence::GetFluxVectorDiff,
    //                            this);

    // workaround for bug in DiffusionLDG

    m_difffields =
        Array<OneD, MR::ExpListSharedPtr>(m_indfields.size() + n_species);
    int j = 0;
    for (int i = 0; i < m_indfields.size() - 1; ++i, ++j)
    {
        m_difffields[j] = m_indfields[i];
    }
    for (int i = 0; i < n_species + 1; ++i, ++j)
    {
        m_difffields[j] = m_temps[i];
    }

    // m_diffusion->InitObject(m_session, m_difffields);

    // Create storage for velocities

    // ExB velocity
    this->v_ExB = Array<OneD, Array<OneD, NekDouble>>(3);
    this->j_par = Array<OneD, NekDouble>(this->n_pts, 0.0);

    for (int d = 0; d < 3; ++d)
    {
        this->v_ExB[d] = Array<OneD, NekDouble>(this->n_pts, 0.0);
    }

    InitAdvection();

    // Since we are starting from a setup where each field is defined to be a
    // discontinuous field (and thus support DG), the first thing we do is to
    // recreate the phi field so that it is continuous, in order to support the
    // Poisson solve. Note that you can still perform a Poisson solve using a
    // discontinuous field, which is done via the hybridisable discontinuous
    // Galerkin (HDG) approach.

    this->phi   = m_fields[4];
    m_fields[4] = MemoryManager<MR::ContField>::AllocateSharedPtr(
        m_session, m_graph, "phi", true, false);
    m_indfields[phi_idx] = m_fields[4];

    m_ode.DefineOdeRhs(&ElectrostaticTurbulence::DoOdeRhs, this);

    if (this->particles_enabled)
    {
        std::vector<Sym<REAL>> src_syms;
        std::vector<int> src_components;

        int cnt = 0;
        for (const auto &[s, v] : this->GetIons())
        {
            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
            src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_DENSITY"));
            src_components.push_back(0);
            ni_src_idx.push_back(cnt++);

            if (v.fields.find(field_to_index["v"]) != v.fields.end())
            {
                for (int d = 0; d < this->m_spacedim; ++d)
                {
                    this->src_fields.emplace_back(
                        MemoryManager<MR::DisContField>::AllocateSharedPtr(
                            *std::dynamic_pointer_cast<MR::DisContField>(
                                m_fields[0])));
                    src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_MOMENTUM"));
                    src_components.push_back(d);
                }
                vi_src_idx.push_back(cnt);
                cnt += m_spacedim;
            }
            if (v.fields.find(field_to_index["e"]) != v.fields.end())
            {
                this->src_fields.emplace_back(
                    MemoryManager<MR::DisContField>::AllocateSharedPtr(
                        *std::dynamic_pointer_cast<MR::DisContField>(
                            m_fields[0])));

                src_syms.push_back(Sym<REAL>(v.name + "_SOURCE_ENERGY"));
                src_components.push_back(0);
                ei_src_idx.push_back(cnt++);
            }
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

/**
 * @brief Initialise the advection object.
 */
void ElectrostaticTurbulence::InitAdvection()
{
    for (const auto &[s, v] : this->GetSpecies())
    {
        if (v.fields.find(field_to_index["v"]) != v.fields.end())
        {
            int ni_idx = v.fields.at(field_to_index["n"]);
            int vi_idx = v.fields.at(field_to_index["v"]);
            this->advected_fields.push_back(ni_idx);
            this->advected_fields.push_back(vi_idx);

            if (v.fields.find(field_to_index["e"]) != v.fields.end())
            {
                int ei_idx = v.fields.at(field_to_index["e"]);
                this->advected_fields.push_back(ei_idx);
            }
        }
    }

    this->advected_fields.push_back(ee_idx);
    this->advected_fields.push_back(omega_idx);

    m_advfields = Array<OneD, MR::ExpListSharedPtr>(advected_fields.size());
    for (int a : this->advected_fields)
    {
        m_advfields[a] = m_indfields[a];
    }
    // Per-field advection velocities (phi not advected, omega is calculates
    // separately)
    this->adv_vel =
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(this->n_species + 2);
    this->dia_v =
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(this->n_species + 1);
    for (int i = 0; i < this->adv_vel.size(); ++i)
    {
        this->adv_vel[i] = Array<OneD, Array<OneD, NekDouble>>(3);
        for (int d = 0; d < 3; ++d)
        {
            this->adv_vel[i][d] = Array<OneD, NekDouble>(this->n_pts, 0.0);
        }
    }
    for (int i = 0; i < this->dia_v.size(); ++i)
    {
        this->dia_v[i] = Array<OneD, Array<OneD, NekDouble>>(3);
        for (int d = 0; d < 3; ++d)
        {
            this->dia_v[i][d] = Array<OneD, NekDouble>(this->n_pts, 0.0);
        }
    }
    this->omega_flux = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);

    for (int d = 0; d < m_spacedim; ++d)
    {
        this->omega_flux[d] = Array<OneD, NekDouble>(this->n_pts, 0.0);
    }

    if (m_indfields[0]->GetTrace())
    {
        auto nTrace = GetTraceNpoints();

        this->trace_vel_norm =
            Array<OneD, Array<OneD, NekDouble>>(m_advfields.size());
        this->trace_b_norm = Array<OneD, NekDouble>(nTrace, 0.0);

        for (int i = 0; i < this->trace_vel_norm.size(); ++i)
        {
            this->trace_vel_norm[i] = Array<OneD, NekDouble>(nTrace, 0.0);
        }
        this->omega_flux_norm = Array<OneD, NekDouble>(nTrace, 0.0);
        this->omega_flux_trace =
            Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int d = 0; d < m_spacedim; ++d)
        {
            this->omega_flux_trace[d] = Array<OneD, NekDouble>(nTrace, 0.0);
        }
    }

    // Create Riemann solver and set normal velocity
    // callback functions
    m_session->LoadSolverInfo("UpwindType", this->riemann_solver_type,
                              "VorticityUpwind");
    this->riemann_solver = SU::GetRiemannSolverFactory().CreateInstance(
        this->riemann_solver_type, m_session);
    auto t = std::dynamic_pointer_cast<PlasmaSolver>(this->riemann_solver);
    t->omega_idx = omega_idx;
    t->m_system  = as<PlasmaSystem>();

    this->riemann_solver->SetVector(
        "Vn", &ElectrostaticTurbulence::GetAdvVelNorm, this);
    this->riemann_solver->SetScalar(
        "wf", &ElectrostaticTurbulence::GetOmegaFlux, this);

    this->dia_riemann_solver = SU::GetRiemannSolverFactory().CreateInstance(
        "VorticityAverage", m_session);
    t = std::dynamic_pointer_cast<PlasmaSolver>(this->dia_riemann_solver);
    t->omega_idx = omega_idx;
    t->m_system  = as<PlasmaSystem>();

    this->dia_riemann_solver->SetVector(
        "Vn", &ElectrostaticTurbulence::GetAdvVelNorm, this);
    this->dia_riemann_solver->SetScalar(
        "wf", &ElectrostaticTurbulence::GetOmegaFlux, this);

    // Setup advection object
    m_session->LoadSolverInfo("AdvectionType", this->adv_type, "WeakDG");
    m_advection = SU::GetAdvectionFactory().CreateInstance(this->adv_type,
                                                           this->adv_type);
    m_advection->SetFluxVector(&ElectrostaticTurbulence::GetFluxVector, this);
    m_advection->SetRiemannSolver(this->riemann_solver);
    m_advection->InitObject(m_session, m_indfields);

    m_omega_advection = std::make_shared<OmegaAdvection>();
    m_omega_advection->SetFluxVector(&ElectrostaticTurbulence::CalcOmegaFlux,
                                     this);

    // workaround for bug in DiffusionLDG
    // m_difffields = Array<OneD, MR::ExpListSharedPtr>(m_indfields.size() - 1);
    // for (int f = 0; f < m_difffields.size(); ++f)
    // {
    //     m_difffields[f] = m_indfields[f];
    // }
    m_omega_advection->omega_idx = omega_idx;
    m_omega_advection->InitObject(m_session, m_indfields[omega_idx]);
}

bool ElectrostaticTurbulence::v_PostIntegrate(int step)
{
    m_fields[0]->FwdTrans(m_fields[0]->GetPhys(), m_fields[0]->UpdateCoeffs());
    m_fields[1]->FwdTrans(m_fields[1]->GetPhys(), m_fields[1]->UpdateCoeffs());

    // Writes a step of the particle trajectory.

    return PlasmaSystem::v_PostIntegrate(step);
}

/**
 * @brief Populate rhs array ( @p outarray )
 *
 * @param inarray physical values of all fields
 * @param[out] outarray output array (RHSs of time integration equations)
 */
void ElectrostaticTurbulence::DoOdeRhs(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time)
{
    // Get field indices
    int nTracePts = GetTraceTotPoints();
    for (int f = 0; f < outarray.size(); ++f)
    {
        Vmath::Zero(this->n_pts, outarray[f], 1);
    }

    int nvariables = inarray.size();

    m_varConv->GetElectronDensity(inarray, m_fields[0]->UpdatePhys());

    // Store forwards/backwards space along trace space
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables - 1);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvariables - 1);

    for (int i = 0; i < nvariables - 1; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        m_indfields[i]->GetFwdBwdTracePhys(inarray[i], Fwd[i], Bwd[i]);
    }

    // CalcOmegaFlux(inarray, this->omega_flux, this->omega_flux_norm);

    // Perform advection
    DoAdvection(inarray, outarray, time, Fwd, Bwd);

    AddForces(inarray, outarray);

    m_bndConds->Update(inarray, time);

    // CalcKappaTensor();

    // Perform Diffusion
    // DoDiffusion(inarray, outarray, Fwd, Bwd);

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

void ElectrostaticTurbulence::ApplyOmegaBC(
    const Array<OneD, Array<OneD, NekDouble>> &inarray, const NekDouble time)
{
    auto BndExps = m_indfields[omega_idx]->GetBndCondExpansions();

    for (size_t r = 0; r < BndExps.size(); ++r)
    {
        MultiRegions::ExpListSharedPtr bndExp = BndExps[r];

        int nEdgePts    = bndExp->GetTotPoints();
        int nEdgeCoeffs = bndExp->GetNcoeffs();

        MultiRegions::ExpListSharedPtr bndElmtExp;
        m_indfields[omega_idx]->GetBndElmtExpansion(r, bndElmtExp, false);

        Array<OneD, Array<OneD, NekDouble>> B_bndelmt(3);

        for (int d = 0; d < 3; ++d)
        {
            this->B[d]->ExtractPhysToBndElmt(r, this->B[d]->GetPhys(),
                                             B_bndelmt[d]);
        }

        int bnd_elmt_pts = B_bndelmt[0].size();
        Array<OneD, Array<OneD, NekDouble>> w(3);
        Array<OneD, Array<OneD, NekDouble>> v_bndelmt(3);
        for (int d = 0; d < 3; ++d)
        {
            w[d] = Array<OneD, NekDouble>(bnd_elmt_pts, 0.0);
        }
        Array<OneD, NekDouble> ni(bnd_elmt_pts);

        for (const auto &[s, v] : this->GetIons())
        {
            int ni_idx = v.fields.at(field_to_index["n"]);
            m_indfields[ni_idx]->ExtractPhysToBndElmt(r, inarray[ni_idx], ni);

            for (int d = 0; d < 3; ++d)
            {
                m_indfields[omega_idx]->ExtractPhysToBndElmt(
                    r, this->adv_vel[s][d], v_bndelmt[d]);
            }

            // Calculate sum of w = (n_i) v0 X b/|B|

            for (int p = 0; p < bnd_elmt_pts; ++p)
            {
                NekDouble mag_B = B_bndelmt[0][p] * B_bndelmt[0][p] +
                                  B_bndelmt[1][p] * B_bndelmt[1][p] +
                                  B_bndelmt[2][p] * B_bndelmt[2][p];
                w[0][p] += v.mass * ni[p] *
                           (v_bndelmt[1][p] * B_bndelmt[2][p] -
                            v_bndelmt[2][p] * B_bndelmt[1][p]) /
                           mag_B;
                w[1][p] += v.mass * ni[p] *
                           (v_bndelmt[2][p] * B_bndelmt[0][p] -
                            v_bndelmt[0][p] * B_bndelmt[2][p]) /
                           mag_B;
                w[2][p] += v.mass * ni[p] *
                           (v_bndelmt[0][p] * B_bndelmt[1][p] -
                            v_bndelmt[1][p] * B_bndelmt[0][p]) /
                           mag_B;
            }
        }
        Array<OneD, NekDouble> Omega(bnd_elmt_pts, 0.0);

        for (int d = 0; d < m_spacedim; ++d)
        {
            bndElmtExp->PhysDeriv(d, w[d], w[d]);
            Vmath::Vadd(bnd_elmt_pts, w[d], 1, Omega, 1, Omega, 1);
        }

        Array<OneD, NekDouble> Omega_bnd(nEdgePts, 0.0);

        m_fields[omega_idx]->ExtractElmtToBndPhys(r, Omega, Omega_bnd);
        bndExp->UpdatePhys() = Omega_bnd;

        bndExp->FwdTransBndConstrained(Omega_bnd, bndExp->UpdateCoeffs());
    }
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void ElectrostaticTurbulence::DoAdvection(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = this->advected_fields.size();
    int nTracePts  = GetTraceTotPoints();

    Array<OneD, Array<OneD, NekDouble>> outarrayAdv(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        outarrayAdv[i] = Array<OneD, NekDouble>(this->n_pts, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> inarrayAdv(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inFwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inBwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        inarrayAdv[i] = inarray[advected_fields[i]];
        inFwd[i]      = pFwd[advected_fields[i]];
        inBwd[i]      = pBwd[advected_fields[i]];
    }

    Array<OneD, Array<OneD, NekDouble>> advVel(m_spacedim);
    // Helmholtz Solve for electrostatic potential
    SolvePhi(inarray, m_fields[0]->GetPhys());

    // Calculate E
    ComputeE();
    // // Calculate ExB, parallel and diamagnetic velocities
    ComputevExB();
    CalcVelocities(inarray, outarray);
    AddDriftVelocities(inarray, outarray);

    m_advection->SetRiemannSolver(this->riemann_solver);

    m_advection->Advect(nvariables, m_advfields, advVel, inarrayAdv,
                        outarrayAdv, time, inFwd, inBwd);

    // m_omega_advection->Advect(m_indfields[omega_idx], inarrayAdv,
    // outarrayAdv,
    //                           inFwd, inBwd);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vsub(this->n_pts, outarray[this->advected_fields[i]], 1,
                    outarrayAdv[i], 1, outarray[this->advected_fields[i]], 1);
    }

    // AddDriftVelocities(inarray, outarray);
    // m_advection->SetRiemannSolver(this->dia_riemann_solver);

    // m_advection->Advect(nvariables, m_advfields, advVel, inarrayAdv,
    //                     outarrayAdv, time, inFwd, inBwd);
    // for (int i = 0; i < nvariables; ++i)
    // {
    //     Vmath::Vsub(this->n_pts, outarray[this->advected_fields[i]], 1,
    //                 outarrayAdv[i], 1, outarray[this->advected_fields[i]],
    //                 1);
    // }
}

/**
 * @brief Add particle sources to the rhs
 */
void ElectrostaticTurbulence::DoParticles(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    // Add contribution to electron energy
    Vmath::Vadd(this->n_pts, outarray[ee_idx], 1,
                this->src_fields[0]->GetPhys(), 1, outarray[ee_idx], 1);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        //  Add contribution to ion density
        Vmath::Vadd(this->n_pts, outarray[ni_idx], 1,
                    this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                    outarray[ni_idx], 1);

        if (v.fields.find(field_to_index["v"]) != v.fields.end())
        {
            int vi_idx = v.fields.at(field_to_index["v"]);

            for (int d = 0; d < m_spacedim; ++d)
            {
                Vmath::Vvtvp(this->n_pts, this->b_unit[d], 1,
                             this->src_fields[vi_src_idx[s] + d]->GetPhys(), 1,
                             outarray[vi_idx], 1, outarray[vi_idx], 1);
            }
        }

        if (v.fields.find(field_to_index["e"]) != v.fields.end())
        {
            int ei_idx = v.fields.at(field_to_index["e"]);

            // Add contribution to ion energy
            Vmath::Vadd(this->n_pts, outarray[ei_idx], 1,
                        this->src_fields[ei_src_idx[s]]->GetPhys(), 1,
                        outarray[ei_idx], 1);

            // Add number density source contribution to ion energy
            Array<OneD, NekDouble> dynamic_energy(this->n_pts);
            m_varConv->GetIonDynamicEnergy(s, v.mass, inarray, dynamic_energy);
            Vmath::Vvtvp(this->n_pts, dynamic_energy, 1,
                         this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                         outarray[ei_idx], 1, outarray[ei_idx], 1);
        }
    }
}

/**
 * @brief Calls HelmSolve to solve for the electric potential
 *
 * @param inarray Array of physical field values
 */
void ElectrostaticTurbulence::SolvePhi(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    [[maybe_unused]] const Array<OneD, NekDouble> &ne)
{
    StdRegions::ConstFactorMap factors;
    // Helmholtz => Poisson (lambda = 0)
    factors[StdRegions::eFactorLambda] = 0.0;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Array<OneD, NekDouble> D(this->n_pts, 0.0);
            for (int p = 0; p < this->n_pts; p++)
            {
                D[p] = -b_unit[i][p] * b_unit[j][p];
                if (i == j)
                {
                    D[p] += 1;
                }
                D[p] /= mag_B[p];
            }

            m_phi_varcoeff[vc[i][j]] = D;
        }
    }

    // Solve for phi. Output of this routine is in coefficient (spectral)
    // space, so backwards transform to physical space since we'll need that
    // for the advection step & computing drift velocity.

    m_indfields[phi_idx]->HelmSolve(
        inarray[omega_idx], this->phi->UpdateCoeffs(), factors, m_phi_varcoeff);

    m_indfields[phi_idx]->BwdTrans(this->phi->GetCoeffs(),
                                   this->phi->UpdatePhys());

    Array<OneD, NekDouble> tmp(this->n_pts, 0.0);
    Array<OneD, NekDouble> tmp2(this->n_pts, 0.0);
    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int ei_idx = v.fields.at(field_to_index["e"]);
        for (int p = 0; p < this->n_pts; ++p)
        {
            tmp2[p] += (2. / 3) * inarray[ei_idx][p] / v.charge;
            tmp[p] += v.mass * inarray[ni_idx][p];
        }
    }
    Array<OneD, NekDouble> &phi = this->phi->UpdatePhys();
    for (int p = 0; p < this->n_pts; ++p)
    {
        phi[p] = (phi[p] - tmp2[p]) / tmp[p];
    }
    m_indfields[phi_idx]->FwdTrans(m_indfields[phi_idx]->GetPhys(),
                                   m_indfields[phi_idx]->UpdateCoeffs());
}

/**
 * @brief Calculates initial potential and gradient
 */
void ElectrostaticTurbulence::CalcInitPhi()
{
    Array<OneD, Array<OneD, NekDouble>> inarray(m_indfields.size());
    Array<OneD, NekDouble> ne = m_fields[0]->UpdatePhys();
    for (int i = 0; i < m_indfields.size(); i++)
    {
        inarray[i] = m_indfields[i]->GetPhys();
    }
    m_varConv->GetElectronDensity(inarray, ne);
    m_fields[0]->FwdTrans(ne, m_fields[0]->UpdateCoeffs());
    CalcVelocities(inarray);
    AddDriftVelocities(inarray);
    CalcInitOmega();

    SolvePhi(inarray, ne);
    ComputeE();
}

/**
 * @brief Compute the gradient of phi for evaluation at the particle positions.
 */
void ElectrostaticTurbulence::ComputeE()
{
    m_indfields[phi_idx]->PhysDeriv(
        this->phi->GetPhys(), this->E[0]->UpdatePhys(),
        this->E[1]->UpdatePhys(), this->E[2]->UpdatePhys());

    Vmath::Neg(this->n_pts, this->E[0]->UpdatePhys(), 1);
    Vmath::Neg(this->n_pts, this->E[1]->UpdatePhys(), 1);
    Vmath::Neg(this->n_pts, this->E[2]->UpdatePhys(), 1);

    this->E[0]->FwdTrans(this->E[0]->GetPhys(), this->E[0]->UpdateCoeffs());
    this->E[1]->FwdTrans(this->E[1]->GetPhys(), this->E[1]->UpdateCoeffs());
    this->E[2]->FwdTrans(this->E[2]->GetPhys(), this->E[2]->UpdateCoeffs());
}

/**
 * @brief Calculate ExB velocity
 */
void ElectrostaticTurbulence::ComputevExB()
{
    const Array<OneD, NekDouble> &Ex = this->E[0]->GetPhys();
    const Array<OneD, NekDouble> &Ey = this->E[1]->GetPhys();
    const Array<OneD, NekDouble> &Ez = this->E[2]->GetPhys();
    const Array<OneD, NekDouble> &Bx = this->B[0]->GetPhys();
    const Array<OneD, NekDouble> &By = this->B[1]->GetPhys();
    const Array<OneD, NekDouble> &Bz = this->B[2]->GetPhys();

    for (int p = 0; p < this->n_pts; ++p)
    {
        this->v_ExB[0][p] = (Ey[p] * Bz[p] - Ez[p] * By[p]) / this->mag_B[p];
        this->v_ExB[1][p] = (Ez[p] * Bx[p] - Ex[p] * Bz[p]) / this->mag_B[p];
        this->v_ExB[2][p] = (Ex[p] * By[p] - Ey[p] * Bx[p]) / this->mag_B[p];
    }
}

/**
 * @brief Calculate advection velocities
 */
void ElectrostaticTurbulence::CalcVelocities(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    for (int f = 0; f < this->adv_vel.size(); ++f)
    {
        for (int d = 0; d < this->adv_vel[f].size(); ++d)
        {
            Vmath::Zero(this->n_pts, adv_vel[f][d], 1);
        }
    }
    const Array<OneD, NekDouble> &ne = m_fields[0]->GetPhys();

    Array<OneD, Array<OneD, NekDouble>> gradv(m_spacedim);
    for (int d = 0; d < m_spacedim; ++d)
    {
        gradv[d] = Array<OneD, NekDouble>(this->n_pts, 0.0);
    }

    // Zero Electron velocity
    Array<OneD, NekDouble> &j_i = m_fields[1]->UpdatePhys();
    Vmath::Zero(this->n_pts, j_i, 1);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int ei_idx = v.fields.at(field_to_index["e"]);
        // Calculate Ion parallel velocities

        for (int p = 0; p < this->n_pts; ++p)
        {
            j_i[p] += v.charge * inarray[vi_idx][p] / v.mass;
            double v_i_par = inarray[vi_idx][p] / (v.mass * inarray[ni_idx][p]);
            for (int d = 0; d < 3; ++d)
            {
                this->adv_vel[s][d][p] =
                    this->v_ExB[d][p] + v_i_par * this->b_unit[d][p];
            }
        }
    }

    Array<OneD, Array<OneD, NekDouble>> gradp(3);
    for (int d = 0; d < 3; ++d)
    {
        gradp[d] = Array<OneD, NekDouble>(this->n_pts, 0.0);
    }
    if (m_spacedim == 3)
        m_indfields[ee_idx]->PhysDeriv(inarray[ee_idx], gradp[0], gradp[1],
                                       gradp[2]);
    else
        m_indfields[ee_idx]->PhysDeriv(inarray[ee_idx], gradp[0], gradp[1]);

    // TODO calculate conductivity
    double sigma = 1;
    for (int p = 0; p < this->n_pts; ++p)
    {
        this->j_par[p] = 0.0;
        for (int d = 0; d < 3; ++d)
        {
            this->j_par[p] += sigma *
                              (this->E[d]->GetPhys()[p] + gradp[d][p] / ne[p]) *
                              this->b_unit[d][p];
        }
    }
    for (int p = 0; p < this->n_pts; ++p)
    {
        for (int d = 0; d < 3; ++d)
        {
            this->adv_vel[n_species][d][p] =
                this->v_ExB[d][p] +
                this->b_unit[d][p] * (j_i[p] - this->j_par[p]) / ne[p];
            this->adv_vel[n_species + 1][d][p] = 0.5 * this->v_ExB[d][p];
        }
    }

    for (const auto &[s, v] : GetNeutrals())
    {
        int nn_idx = v.fields.at(field_to_index["n"]);
        int vn_idx = v.fields.at(field_to_index["v"]);
        int pn_idx = v.fields.at(field_to_index["e"]);
        // Calculate Neutral parallel velocities
        for (int p = 0; p < this->n_pts; ++p)
        {
            double v_n_par = inarray[vn_idx][p] / (v.mass * inarray[nn_idx][p]);
            for (int d = 0; d < m_spacedim; ++d)
            {
                this->adv_vel[s][d][p] = v_n_par * this->b_unit[d][p];
            }
        }
    }
}

void ElectrostaticTurbulence::AddForces(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{

    Array<OneD, NekDouble> gradp(this->n_pts, 0.0);
    Array<OneD, NekDouble> gradv(this->n_pts, 0.0);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int ei_idx = v.fields.at(field_to_index["e"]);

        for (int d = 0; d < m_spacedim; ++d)
        {
            m_indfields[ei_idx]->PhysDeriv(d, this->adv_vel[s][d], gradv);
            m_indfields[ei_idx]->PhysDeriv(inarray[ei_idx], gradp);
            for (int p = 0; p < this->n_pts; ++p)
            {
                outarray[vi_idx][p] +=
                    b_unit[d][p] *
                    (v.charge * inarray[ni_idx][p] * this->E[d]->GetPhys()[p] -
                     (2.0 / 3.0) * gradp[p]);
                outarray[ei_idx][p] -=
                    (2.0 / 3.0) * inarray[ei_idx][p] * gradv[p];
            }
        }
    }
    for (int d = 0; d < m_spacedim; ++d)
    {
        m_indfields[ee_idx]->PhysDeriv(d, this->adv_vel[n_species][d],
                                       gradv);
        for (int p = 0; p < this->n_pts; ++p)
        {
            outarray[ee_idx][p] -= (2.0 / 3.0) * inarray[ee_idx][p] * gradv[p];
        }
    }
}

/**
 * @brief Add drift velocities to the advection velocities
 */
void ElectrostaticTurbulence::AddDriftVelocities(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    const Array<OneD, NekDouble> &ne = m_fields[0]->GetPhys();

    const Array<OneD, NekDouble> &Bx = this->B[0]->GetPhys();
    const Array<OneD, NekDouble> &By = this->B[1]->GetPhys();
    const Array<OneD, NekDouble> &Bz = this->B[2]->GetPhys();

    Array<OneD, Array<OneD, NekDouble>> boverB(3);
    Array<OneD, Array<OneD, NekDouble>> curlb(3);

    for (int d = 0; d < 3; ++d)
    {
        boverB[d] = Array<OneD, NekDouble>(this->n_pts, 0.0);
        curlb[d]  = Array<OneD, NekDouble>(this->n_pts, 0.0);
        for (int p = 0; p < this->n_pts; ++p)
        {
            boverB[d][p] = this->B[d]->GetPhys()[p] / this->mag_B[p];
        }
    }
    Array<OneD, NekDouble> Dummy(this->n_pts);
    Array<OneD, NekDouble> Vx(this->n_pts, 0.0);
    Array<OneD, NekDouble> Uy(this->n_pts, 0.0);
    Array<OneD, NekDouble> Vz(this->n_pts, 0.0);
    Array<OneD, NekDouble> Uz(this->n_pts, 0.0);
    Array<OneD, NekDouble> Wx(this->n_pts, 0.0);
    Array<OneD, NekDouble> Wy(this->n_pts, 0.0);

    this->B[0]->PhysDeriv(boverB[0], Dummy, Uy, Uz);
    this->B[1]->PhysDeriv(boverB[1], Vx, Dummy, Vz);
    this->B[2]->PhysDeriv(boverB[2], Wx, Wy, Dummy);

    for (int p = 0; p < this->n_pts; ++p)
    {
        curlb[0][p] = Wy[p] - Vz[p];
        curlb[1][p] = Uz[p] - Wx[p];
        curlb[2][p] = Vx[p] - Uy[p];
    }

    Array<OneD, Array<OneD, NekDouble>> gradp(3);

    for (int d = 0; d < 3; ++d)
    {
        gradp[d] = Array<OneD, NekDouble>(this->n_pts, 0.0);
    }

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int ei_idx = v.fields.at(field_to_index["e"]);
        for (int p = 0; p < this->n_pts; ++p)
        {
            this->dia_v[s][0][p] = (2.0 / 3.0) * inarray[ei_idx][p] * curlb[0][p] /
                             (v.charge * inarray[ni_idx][p]);
            this->dia_v[s][1][p] = (2.0 / 3.0) * inarray[ei_idx][p] * curlb[1][p] /
                             (v.charge * inarray[ni_idx][p]);
            this->dia_v[s][2][p] = (2.0 / 3.0) * inarray[ei_idx][p] * curlb[2][p] /
                             (v.charge * inarray[ni_idx][p]);
        }
    }

    for (const auto &[s, v] : this->GetNeutrals())
    {
        int nn_idx = v.fields.at(field_to_index["n"]);
        int vn_idx = v.fields.at(field_to_index["v"]);
        int en_idx = v.fields.at(field_to_index["e"]);

        if (m_spacedim == 2)
            m_indfields[en_idx]->PhysDeriv(inarray[en_idx], gradp[0], gradp[1]);
        else if (m_spacedim == 3)
            m_indfields[en_idx]->PhysDeriv(inarray[en_idx], gradp[0], gradp[1],
                                           gradp[2]);

        if (outarray != NullNekDoubleArrayOfArray)
        {
            for (int p = 0; p < this->n_pts; ++p)
            {
                for (int d = 0; d < m_spacedim; ++d)
                {
                    outarray[vn_idx][p] +=
                        b_unit[d][p] * ((2.0 / 3.0) * gradp[d][p]);
                }
            }
        }
    }
    for (int p = 0; p < this->n_pts; ++p)
    {
        this->dia_v[n_species][0][p] = -(2.0 / 3.0) * inarray[ee_idx][p] * curlb[0][p] / ne[p];
        this->dia_v[n_species][1][p] = -(2.0 / 3.0) * inarray[ee_idx][p] * curlb[1][p] / ne[p];
        this->dia_v[n_species][2][p] = -(2.0 / 3.0) * inarray[ee_idx][p] * curlb[2][p] / ne[p];
    }
    
}

/**
 * @brief Calculate initial vorticity
 */
void ElectrostaticTurbulence::CalcInitOmega()
{
    const Array<OneD, NekDouble> &Bx = this->B[0]->GetPhys();
    const Array<OneD, NekDouble> &By = this->B[1]->GetPhys();
    const Array<OneD, NekDouble> &Bz = this->B[2]->GetPhys();
    Vmath::Zero(this->n_pts, m_indfields[omega_idx]->UpdatePhys(), 1);

    Array<OneD, Array<OneD, NekDouble>> w(3);
    for (int d = 0; d < 3; ++d)
    {
        w[d] = Array<OneD, NekDouble>(this->n_pts, 0.0);
    }

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx                       = v.fields.at(field_to_index["n"]);
        const Array<OneD, NekDouble> &ni = m_indfields[ni_idx]->GetPhys();
        // Calculate sum of Zw = (n_i) v0 X b/|B|

        for (int p = 0; p < this->n_pts; ++p)
        {
            double vx = this->adv_vel[s][0][p] + this->dia_v[s][0][p];
            double vy = this->adv_vel[s][1][p] + this->dia_v[s][1][p];
            double vz = this->adv_vel[s][2][p] + this->dia_v[s][2][p];
            w[0][p] +=
                v.mass * ni[p] * (vy * Bz[p] - vz * By[p]) / this->mag_B[p];
            w[1][p] +=
                v.mass * ni[p] * (vz * Bx[p] - vx * Bz[p]) / this->mag_B[p];
            w[2][p] +=
                v.mass * ni[p] * (vx * By[p] - vy * Bx[p]) / this->mag_B[p];
        }
    }
    for (int d = 0; d < m_spacedim; ++d)
    {
        m_indfields[omega_idx]->PhysDeriv(d, w[d], w[d]);

        Vmath::Vadd(this->n_pts, w[d], 1, m_indfields[omega_idx]->GetPhys(), 1,
                    m_indfields[omega_idx]->UpdatePhys(), 1);
    }
    // Vmath::Smul(this->n_pts, 1.0 / this->omega_c,
    //             m_indfields[omega_idx]->GetPhys(), 1,
    //             m_indfields[omega_idx]->UpdatePhys(), 1);

    m_indfields[omega_idx]->FwdTransLocalElmt(
        m_indfields[omega_idx]->GetPhys(),
        m_indfields[omega_idx]->UpdateCoeffs());
}

/**
 * @brief Calculate vorticity flux
 */
void ElectrostaticTurbulence::CalcOmegaFlux(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &volume_flux,
    Array<OneD, NekDouble> &trace_flux)
{
    size_t n_tracepoints = GetTraceNpoints();
    size_t ncoeffs       = GetNcoeffs();

    const Array<OneD, NekDouble> &Bx = this->B[0]->GetPhys();
    const Array<OneD, NekDouble> &By = this->B[1]->GetPhys();
    const Array<OneD, NekDouble> &Bz = this->B[2]->GetPhys();
    for (int d = 0; d < m_spacedim; ++d)
    {
        Vmath::Zero(this->n_pts, volume_flux[d], 1);
        // Vmath::Zero(ncoeffs, omega_flux[d], 1);
    }
    Vmath::Zero(n_tracepoints, trace_flux, 1);

    Array<OneD, Array<OneD, NekDouble>> normals(m_spacedim);

    for (int d = 0; d < m_spacedim; ++d)
    {
        normals[d] = Array<OneD, NekDouble>(n_tracepoints);
    }

    m_indfields[omega_idx]->GetTrace()->GetNormals(normals);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);

        // Calculate w = (n_i) v0 X b/|B|
        Array<OneD, Array<OneD, NekDouble>> w(m_spacedim);
        if (m_spacedim == 3)
        {
            w[0] = Array<OneD, NekDouble>(this->n_pts);
            w[1] = Array<OneD, NekDouble>(this->n_pts);
            w[2] = Array<OneD, NekDouble>(this->n_pts);
            for (int p = 0; p < this->n_pts; ++p)
            {
                w[0][p] = inarray[ni_idx][p] *
                          (this->adv_vel[s][1][p] * Bz[p] -
                           this->adv_vel[s][2][p] * By[p]) /
                          this->mag_B[p];
                w[1][p] = inarray[ni_idx][p] *
                          (this->adv_vel[s][2][p] * Bx[p] -
                           this->adv_vel[s][0][p] * Bz[p]) /
                          this->mag_B[p];
                w[2][p] = inarray[ni_idx][p] *
                          (this->adv_vel[s][0][p] * By[p] -
                           this->adv_vel[s][1][p] * Bx[p]) /
                          this->mag_B[p];
            }
        }
        else if (m_spacedim == 2)
        {
            w[0] = Array<OneD, NekDouble>(this->n_pts);
            w[1] = Array<OneD, NekDouble>(this->n_pts);
            for (int p = 0; p < this->n_pts; ++p)
            {
                w[0][p] = inarray[ni_idx][p] *
                          (this->adv_vel[s][1][p] * Bz[p] -
                           this->adv_vel[s][2][p] * By[p]) /
                          this->mag_B[p];
                w[1][p] = inarray[ni_idx][p] *
                          (this->adv_vel[s][2][p] * Bx[p] -
                           this->adv_vel[s][0][p] * Bz[p]) /
                          this->mag_B[p];
            }
        }

        for (int d = 0; d < m_spacedim; ++d)
        {
            Array<OneD, NekDouble> tmp(this->n_pts, 0.0);

            for (int d2 = 0; d2 < m_spacedim; ++d2)
            {
                Array<OneD, NekDouble> vw(this->n_pts, 0.0);
                Array<OneD, NekDouble> tmp_c(ncoeffs, 0.0);
                Array<OneD, NekDouble> Fwd(n_tracepoints, 0.0);
                Array<OneD, NekDouble> Bwd(n_tracepoints, 0.0);

                // Calculate ∇⋅(v0⊗w)
                Vmath::Vmul(this->n_pts, this->adv_vel[s][d], 1, w[d2], 1, vw,
                            1);
                m_indfields[omega_idx]->PhysDeriv(d2, vw, vw);

                // m_indfields[omega_idx]->GetFwdBwdTracePhys(vw, Fwd, Bwd);
                // Vmath::Vmul(n_tracepoints, normals[d2], 1, Fwd, 1, Fwd, 1);
                // m_indfields[omega_idx]->IProductWRTDerivBase(d2, vw, tmp_c);
                // Vmath::Neg(ncoeffs, tmp_c, 1);
                // m_indfields[omega_idx]->AddTraceIntegral(Fwd, tmp_c);
                // m_indfields[omega_idx]->MultiplyByElmtInvMass(tmp_c, tmp_c);
                // m_indfields[omega_idx]->BwdTrans(tmp_c, vw);

                Vmath::Vadd(this->n_pts, vw, 1, tmp, 1, tmp, 1);
            }
            // Vorticity Flux
            Vmath::Smul(this->n_pts, v.mass, tmp, 1, tmp, 1);
            Vmath::Vadd(this->n_pts, tmp, 1, volume_flux[d], 1, volume_flux[d],
                        1);
        }

        Array<OneD, NekDouble> tmp(this->n_pts, 0.0);
        Array<OneD, NekDouble> divw(this->n_pts, 0.0);

        // // Calculate ∇⋅w
        // for (int d = 0; d < m_spacedim; ++d)
        // {
        //     m_indfields[omega_idx]->PhysDeriv(d, w[d], tmp);
        //     Vmath::Vadd(this->n_pts, tmp, 1, divw, 1, divw, 1);
        // }
        // Vmath::Smul(this->n_pts, v.mass, divw, 1, divw, 1);

        // for (int d = 0; d < m_spacedim; ++d)
        // {
        //     Vmath::Vmul(this->n_pts, this->adv_vel[ni_idx][d], 1, divw, 1,
        //     tmp,
        //                 1);

        //     // Vorticity Flux
        //     Vmath::Vadd(this->n_pts, tmp, 1, volume_flux[d], 1,
        //     volume_flux[d],
        //                 1);
        // }
        Array<OneD, NekDouble> trace_vel(n_tracepoints, 0.0);
        Array<OneD, NekDouble> trace_vel_norm(n_tracepoints, 0.0);

        Array<OneD, NekDouble> Fwd(n_tracepoints, 0.0);
        Array<OneD, NekDouble> Bwd(n_tracepoints, 0.0);

        m_indfields[omega_idx]->GetFwdBwdTracePhys(inarray[omega_idx], Fwd,
                                                   Bwd);
        for (int d = 0; d < m_spacedim; ++d)
        {
            m_indfields[omega_idx]->ExtractTracePhys(this->v_ExB[d], trace_vel);
            for (int p = 0; p < n_tracepoints; ++p)
            {
                trace_vel_norm[p] += normals[d][p] * trace_vel[p];
            }
        }
        for (int p = 0; p < n_tracepoints; ++p)
        {
            trace_flux[p] +=
                (trace_vel_norm[p] > 0 ? 0.5 * Fwd[p] * trace_vel_norm[p]
                                       : 0.5 * Bwd[p] * trace_vel_norm[p]);
        }

        // Array<OneD, NekDouble> trace_flux_d(n_tracepoints, 0.0);
        // Array<OneD, NekDouble> Fwd(n_tracepoints, 0.0);
        // Array<OneD, NekDouble> Bwd(n_tracepoints, 0.0);
        // Array<OneD, NekDouble> qFwd(n_tracepoints, 0.0);
        // Array<OneD, NekDouble> qBwd(n_tracepoints, 0.0);
        // Array<OneD, NekDouble> uterm(n_tracepoints, 0.0);

        // m_indfields[omega_idx]->GetFwdBwdTracePhys(inarray[omega_idx], Fwd,
        //                                            Bwd);
        // Vmath::Vsub(n_tracepoints, Fwd, 1, Bwd, 1, uterm, 1);
        // Vmath::Smul(n_tracepoints, 100.0, uterm, 1, uterm, 1);

        // for (int d = 0; d < m_spacedim; ++d)
        // {
        //     // m_indfields[omega_idx]->GetFwdBwdTracePhys(volume_flux[d],
        //     qFwd,
        //     //                                            qBwd);
        //     m_indfields[omega_idx]->ExtractTracePhys(volume_flux[d],
        //                                              trace_flux_d);
        //     for (int p = 0; p < n_tracepoints; ++p)
        //     {
        //         // trace_flux[p] += normals[d][p] * qFwd[p];
        //         trace_flux[p] += normals[d][p] * trace_flux_d[p];
        //     }
        //     // Vmath::Vadd(n_tracepoints, uterm, 1, trace_flux, 1,
        //     trace_flux,
        //     // 1);
        // }
    }
}

/**
 *  @brief Compute components of advection velocities normal to trace
 * elements (faces, in 3D).
 */
Array<OneD, Array<OneD, NekDouble>> &ElectrostaticTurbulence::GetAdvVelNorm()
{
    // Number of trace (interface) points
    int num_trace_pts = GetTraceNpoints();
    Array<OneD, Array<OneD, NekDouble>> normals(m_spacedim);
    Array<OneD, NekDouble> tmp(num_trace_pts);

    for (int d = 0; d < m_spacedim; ++d)
    {
        normals[d] = Array<OneD, NekDouble>(num_trace_pts);
    }
    // Compute advection vel dot trace normals and store
    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int ei_idx = v.fields.at(field_to_index["e"]);
        m_indfields[ni_idx]->GetTrace()->GetNormals(normals);
        // Ensure output array is zeroed
        Vmath::Zero(num_trace_pts, this->trace_vel_norm[ni_idx], 1);
        for (int d = 0; d < m_spacedim; ++d)
        {
            m_indfields[ni_idx]->ExtractTracePhys(this->adv_vel[s][d], tmp);
            for (int p = 0; p < num_trace_pts; ++p)
            {
                this->trace_vel_norm[ni_idx][p] += normals[d][p] * tmp[p];
            }
        }

        m_indfields[vi_idx]->GetTrace()->GetNormals(normals);
        // Ensure output array is zeroed
        Vmath::Zero(num_trace_pts, this->trace_vel_norm[vi_idx], 1);
        for (int d = 0; d < m_spacedim; ++d)
        {
            m_indfields[vi_idx]->ExtractTracePhys(this->adv_vel[s][d], tmp);
            for (int p = 0; p < num_trace_pts; ++p)
            {
                this->trace_vel_norm[vi_idx][p] += normals[d][p] * tmp[p];
            }
        }

        m_indfields[ei_idx]->GetTrace()->GetNormals(normals);
        // Ensure output array is zeroed
        Vmath::Zero(num_trace_pts, this->trace_vel_norm[ei_idx], 1);
        for (int d = 0; d < m_spacedim; ++d)
        {
            m_indfields[ei_idx]->ExtractTracePhys(this->adv_vel[s][d], tmp);
            for (int p = 0; p < num_trace_pts; ++p)
            {
                this->trace_vel_norm[ei_idx][p] += normals[d][p] * tmp[p];
            }
        }
    }
    m_indfields[ee_idx]->GetTrace()->GetNormals(normals);
    // Ensure output array is zeroed
    Vmath::Zero(num_trace_pts, this->trace_vel_norm[ee_idx], 1);
    for (int d = 0; d < m_spacedim; ++d)
    {
        m_indfields[ee_idx]->ExtractTracePhys(this->adv_vel[n_species][d], tmp);
        for (int p = 0; p < num_trace_pts; ++p)
        {
            this->trace_vel_norm[ee_idx][p] += normals[d][p] * tmp[p];
        }
    }
    m_indfields[omega_idx]->GetTrace()->GetNormals(normals);
    // Ensure output array is zeroed
    Vmath::Zero(num_trace_pts, this->trace_vel_norm[omega_idx], 1);
    for (int d = 0; d < m_spacedim; ++d)
    {
        m_indfields[omega_idx]->ExtractTracePhys(
            this->adv_vel[n_species + 1][d], tmp);
        for (int p = 0; p < num_trace_pts; ++p)
        {
            this->trace_vel_norm[omega_idx][p] += normals[d][p] * tmp[p];
        }
    }
    return this->trace_vel_norm;
}

/**
 * @brief Fetch the flux of vorticity
 */
Array<OneD, NekDouble> &ElectrostaticTurbulence::GetOmegaFlux()
{
    // int num_trace_pts = GetTraceNpoints();
    // Array<OneD, Array<OneD, NekDouble>> normals(m_spacedim);

    // for (int d = 0; d < m_spacedim; ++d)
    // {
    //     normals[d] = Array<OneD, NekDouble>(num_trace_pts);
    // }

    // m_indfields[omega_idx]->GetTrace()->GetNormals(normals);
    // Vmath::Zero(num_trace_pts, this->omega_flux_norm, 1);

    // for (int d = 0; d < m_spacedim; ++d)
    // {
    //     m_indfields[omega_idx]->ExtractTracePhys(this->omega_flux[d],
    //                                              this->omega_flux_trace[d]);
    //     for (int p = 0; p < num_trace_pts; ++p)
    //     {
    //         this->omega_flux_norm[p] +=
    //             normals[d][p] * this->omega_flux_trace[d][p];
    //     }
    // }

    return this->omega_flux_norm;
}

/**
 *  @brief Construct flux array.
 *
 * @param  field_vals Physical values for each advection field
 * @param[out] flux       Flux array
 */
void ElectrostaticTurbulence::GetFluxVector(
    const Array<OneD, Array<OneD, NekDouble>> &field_vals,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    for (const auto &[s, v] : this->GetSpecies())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int ei_idx = v.fields.at(field_to_index["e"]);

        for (int d = 0; d < m_spacedim; ++d)
        {
            for (int p = 0; p < this->n_pts; ++p)
            {
                double v  = this->adv_vel[s][d][p];
                double vd = this->dia_v[s][d][p];

                fluxes[ni_idx][d][p] = (v + vd) * field_vals[ni_idx][p];
                fluxes[vi_idx][d][p] = (v + vd) * field_vals[vi_idx][p];
                fluxes[ei_idx][d][p] =
                    (v + (5.0 / 3.0) * vd) * field_vals[ei_idx][p];
            }
        }
    }
    for (int d = 0; d < m_spacedim; ++d)
    {
        for (int p = 0; p < this->n_pts; ++p)
        {
            double v  = this->adv_vel[n_species][d][p];
            double vd = this->dia_v[n_species][d][p];

            fluxes[ee_idx][d][p] =
                (v + (5.0 / 3.0) * vd) * field_vals[ee_idx][p];
            fluxes[omega_idx][d][p] =
                this->adv_vel[n_species + 1][d][p] * field_vals[omega_idx][p];
        }
    }
    // for (int d = 0; d < m_spacedim; ++d)
    // {
    //     for (int d2 = 0; d2 < m_spacedim; ++d2)
    //     {
    //         auto &Ek = this->E[d2]->GetPhys();
    //         for (int p = 0; p < this->n_pts; ++p)
    //         {
    //             fluxes[omega_idx][d][p] += Ek[p] * b_unit[d2][p] *
    //             b_unit[d][p];
    //         }
    //     }
    // }
    // Omega flux
    // for (int d = 0; d < m_spacedim; ++d)
    // {
    //     Vmath::Vcopy(this->n_pts, this->omega_flux[d], 1,
    //     fluxes[omega_idx][d],
    //                  1);
    // }
}

/**
 * @brief Add diffusion to the rhs
 */
void ElectrostaticTurbulence::DoDiffusion(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = inarray.size() - 1;
    int nTracePts  = GetTraceTotPoints();

    // this should be preallocated
    Array<OneD, Array<OneD, NekDouble>> outarrayDiff(nvariables);
    for (int i = 0; i < nvariables; ++i)
    {
        outarrayDiff[i] = Array<OneD, NekDouble>(this->n_pts, 0.0);
    }

    Array<OneD, Array<OneD, NekDouble>> inarrayprim(nvariables + n_species + 1);
    Array<OneD, Array<OneD, NekDouble>> inFwdprim(nvariables + n_species + 1);
    Array<OneD, Array<OneD, NekDouble>> inBwdprim(nvariables + n_species + 1);

    int t = 0;
    for (const auto &[s, v] : this->GetSpecies())
    {
        int ei_idx          = v.fields.at(field_to_index["e"]);
        int ni_idx          = v.fields.at(field_to_index["n"]);
        int vi_idx          = v.fields.at(field_to_index["v"]);
        int ti_idx          = nvariables + t++;
        inarrayprim[ni_idx] = Array<OneD, NekDouble>(inarray[ni_idx]);
        inFwdprim[ni_idx]   = Array<OneD, NekDouble>(pFwd[ni_idx]);
        inBwdprim[ni_idx]   = Array<OneD, NekDouble>(pBwd[ni_idx]);

        inarrayprim[ei_idx] = Array<OneD, NekDouble>(inarray[ei_idx]);
        inFwdprim[ei_idx]   = Array<OneD, NekDouble>(pFwd[ei_idx]);
        inBwdprim[ei_idx]   = Array<OneD, NekDouble>(pBwd[ei_idx]);

        inarrayprim[ti_idx] = Array<OneD, NekDouble>(m_temps[t]->GetPhys());
        inFwdprim[ti_idx]   = Array<OneD, NekDouble>(nTracePts, 0.0);
        inBwdprim[ti_idx]   = Array<OneD, NekDouble>(nTracePts, 0.0);
        m_varConv->GetIonTemperature(s, v.mass, inarray, inarrayprim[ti_idx]);
        m_varConv->GetIonTemperature(s, v.mass, pFwd, inFwdprim[ti_idx]);
        m_varConv->GetIonTemperature(s, v.mass, pBwd, inBwdprim[ti_idx]);

        inarrayprim[vi_idx] = Array<OneD, NekDouble>(this->n_pts, 0.0);
        inFwdprim[vi_idx]   = Array<OneD, NekDouble>(nTracePts, 0.0);
        inBwdprim[vi_idx]   = Array<OneD, NekDouble>(nTracePts, 0.0);

        for (int p = 0; p < nTracePts; ++p)
        {
            inFwdprim[vi_idx][p] = pFwd[vi_idx][p] / (v.mass * pFwd[vi_idx][p]);
            inBwdprim[vi_idx][p] = pBwd[vi_idx][p] / (v.mass * pBwd[vi_idx][p]);
        }
    }

    inarrayprim[ee_idx] = Array<OneD, NekDouble>(inarray[ee_idx]);
    inFwdprim[ee_idx]   = Array<OneD, NekDouble>(pFwd[ee_idx]);
    inBwdprim[ee_idx]   = Array<OneD, NekDouble>(pBwd[ee_idx]);
    int te_idx          = nvariables + t;

    inarrayprim[te_idx] = Array<OneD, NekDouble>(m_temps[t]->GetPhys());
    inFwdprim[te_idx]   = Array<OneD, NekDouble>(nTracePts, 0.0);
    inBwdprim[te_idx]   = Array<OneD, NekDouble>(nTracePts, 0.0);
    t++;
    // Extract temperature

    m_varConv->GetElectronTemperature(inarray, inarrayprim[te_idx]);
    m_varConv->GetElectronTemperature(pFwd, inFwdprim[te_idx]);
    m_varConv->GetElectronTemperature(pBwd, inBwdprim[te_idx]);

    // CalcKappaTensor();
    m_diffusion->Diffuse(nvariables + t, m_difffields, inarrayprim,
                         outarrayDiff, inFwdprim, inBwdprim);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(this->n_pts, outarrayDiff[i], 1, outarray[i], 1,
                    outarray[i], 1);
    }
}

void ElectrostaticTurbulence::CalcKPar()
{
    // Change to fn of fields
    NekDouble k_par;
    m_session->LoadParameter("k_par", k_par, 100.0);
    m_kpar = Array<OneD, NekDouble>(this->n_pts, k_par);
}

void ElectrostaticTurbulence::CalcKPerp()
{
    // Change to fn of fields
    NekDouble k_perp;
    m_session->LoadParameter("k_perp", k_perp, 1.0);
    m_kperp = Array<OneD, NekDouble>(this->n_pts, k_perp);
}

void ElectrostaticTurbulence::CalcDiffTensor()
{
    CalcKPar();
    CalcKPerp();
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Array<OneD, NekDouble> d(this->n_pts, 0.0);
            for (int k = 0; k < this->n_pts; k++)
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
    NekDouble kappa_par;
    m_session->LoadParameter("kappa_par", kappa_par, 0.);
    m_kappapar = Array<OneD, NekDouble>(this->n_pts, kappa_par);
}

void ElectrostaticTurbulence::CalcKappaPerp()
{
    // Change to fn of T
    NekDouble kappa_perp;
    m_session->LoadParameter("kappa_perp", kappa_perp, 0.01);
    m_kappaperp = Array<OneD, NekDouble>(this->n_pts, kappa_perp);
}

void ElectrostaticTurbulence::CalcKappaTensor()
{
    CalcKappaPar();
    CalcKappaPerp();
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Array<OneD, NekDouble> kappa(this->n_pts, 0.0);
            for (int k = 0; k < this->n_pts; k++)
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
void ElectrostaticTurbulence::GetFluxVectorDiff(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    const auto &Bx = this->B[0]->GetPhys();
    const auto &By = this->B[1]->GetPhys();
    const auto &Bz = this->B[2]->GetPhys();

    int t = 0;
    for (const auto &[s, v] : this->GetIons())
    {
        int vi_idx = v.fields.at(field_to_index["v"]);
        int ei_idx = v.fields.at(field_to_index["e"]);
        int ni_idx = v.fields.at(field_to_index["n"]);
        int ti_idx = inarray.size() - n_species - 1 + t++;

        if (m_spacedim == 3)
        {
            for (int p = 0; p < this->n_pts; ++p)
            {
                double qx         = (2. / 3) * qfield[0][ei_idx][p] / v.charge;
                double qy         = (2. / 3) * qfield[1][ei_idx][p] / v.charge;
                double qz         = (2. / 3) * qfield[2][ei_idx][p] / v.charge;
                double dia_flux_x = (-Bz[p] * qy + By[p] * qz) / this->mag_B[p];
                double dia_flux_y = (Bz[p] * qx - Bx[p] * qz) / this->mag_B[p];
                double dia_flux_z = (-By[p] * qx + Bx[p] * qy) / this->mag_B[p];

                fluxes[0][ni_idx][p] += dia_flux_x;
                fluxes[1][ni_idx][p] += dia_flux_y;
                fluxes[2][ni_idx][p] += dia_flux_z;

                fluxes[0][vi_idx][p] +=
                    v.mass * inarray[vi_idx][p] * dia_flux_x;
                fluxes[1][vi_idx][p] +=
                    v.mass * inarray[vi_idx][p] * dia_flux_y;
                fluxes[2][vi_idx][p] +=
                    v.mass * inarray[vi_idx][p] * dia_flux_z;

                fluxes[0][ei_idx][p] += 2.5 * inarray[ti_idx][p] * dia_flux_x;
                fluxes[1][ei_idx][p] += 2.5 * inarray[ti_idx][p] * dia_flux_y;
                fluxes[2][ei_idx][p] += 2.5 * inarray[ti_idx][p] * dia_flux_z;
            }
        }
        else
        {
            for (int p = 0; p < this->n_pts; ++p)
            {
                double qx         = (2. / 3) * qfield[0][ei_idx][p] / v.charge;
                double qy         = (2. / 3) * qfield[1][ei_idx][p] / v.charge;
                double dia_flux_x = (-Bz[p] * qy) / this->mag_B[p];
                double dia_flux_y = (Bz[p] * qx) / this->mag_B[p];

                // fluxes[0][ni_idx][p] += dia_flux_x;
                // fluxes[1][ni_idx][p] += dia_flux_y;

                // fluxes[0][vi_idx][p] +=
                //     v.mass * inarray[vi_idx][p] * dia_flux_x;
                // fluxes[1][vi_idx][p] +=
                //     v.mass * inarray[vi_idx][p] * dia_flux_y;

                // fluxes[0][ei_idx][p] += 2.5 * inarray[ti_idx][p] *
                // dia_flux_x; fluxes[1][ei_idx][p] += 2.5 * inarray[ti_idx][p]
                // * dia_flux_y;
            }
        }
    }
    int te_idx = inarray.size() - 1;

    // if (m_spacedim == 3)
    // {
    //     for (int p = 0; p < this->n_pts; ++p)
    //     {
    //         double qx = (2. / 3) * qfield[0][ee_idx][p];
    //         double qy = (2. / 3) * qfield[1][ee_idx][p];
    //         double qz = (2. / 3) * qfield[2][ee_idx][p];
    //         double dia_flux_x = (-Bz[p] * qy + By[p] * qz) / this->mag_B[p];
    //         double dia_flux_y = (Bz[p] * qx - Bx[p] * qz) / this->mag_B[p];
    //         double dia_flux_z = (-By[p] * qx + Bx[p] * qy) / this->mag_B[p];

    //         fluxes[0][ee_idx][p] += 2.5 * inarray[te_idx][p] * dia_flux_x;
    //         fluxes[1][ee_idx][p] += 2.5 * inarray[te_idx][p] * dia_flux_y;
    //         fluxes[2][ee_idx][p] += 2.5 * inarray[te_idx][p] * dia_flux_z;
    //     }
    // }
    // else
    // {
    //     for (int p = 0; p < this->n_pts; ++p)
    //     {
    //         double qx = (2. / 3) * qfield[0][ee_idx][p];
    //         double qy = (2. / 3) * qfield[1][ee_idx][p];

    //         double dia_flux_x = (-Bz[p] * qy) / this->mag_B[p];
    //         double dia_flux_y = (Bz[p] * qx) / this->mag_B[p];

    //         fluxes[0][ee_idx][p] += 2.5 * inarray[te_idx][p] * dia_flux_x;
    //         fluxes[1][ee_idx][p] += 2.5 * inarray[te_idx][p] * dia_flux_y;
    //     }
    // }

    // for (int j = 0; j < m_spacedim; ++j)
    // {
    //     for (int k = 0; k < m_spacedim; k++)
    //     {
    //         for (int p = 0; p < this->n_pts; ++p)
    //         {
    //             fluxes[j][omega_idx][p] += m_zeta *
    //                                        (1. - b_unit[k][p] * b_unit[j][p])
    //                                        * qfield[k][omega_idx][p];
    //         }
    //     }
    // }
}

// void ElectrostaticTurbulence::GetFluxPenalty(
//     const Array<OneD, const Array<OneD, NekDouble>> &uFwd,
//     const Array<OneD, const Array<OneD, NekDouble>> &uBwd,
//     Array<OneD, Array<OneD, NekDouble>> &penaltyCoeff)
// {
//     size_t nTracePts = uFwd[0].size();

//     // Compute average temperature
//     size_t nVariables = uFwd.size();
//     Array<OneD, NekDouble> tAve{nTracePts, 0.0};
//     Vmath::Svtsvtp(nTracePts, 0.5, uFwd[nVariables - 1], 1, 0.5,
//                    uBwd[nVariables - 1], 1, tAve, 1);

//     // Get average viscosity and thermal conductivity
//     Array<OneD, NekDouble> muAve{nTracePts, 0.0};
//     Array<OneD, NekDouble> tcAve{nTracePts, 0.0};

//     GetViscosityAndThermalCondFromTemp(tAve, muAve, tcAve);

//     // Compute penalty term
//     for (size_t i = 0; i < nVariables; ++i)
//     {
//         // Get jump of u variables
//         Vmath::Vsub(nTracePts, uFwd[i], 1, uBwd[i], 1, penaltyCoeff[i], 1);
//         // Multiply by variable coefficient = {coeff} ( u^+ - u^- )
//         if (i < nVariables - 1)
//         {
//             Vmath::Vmul(nTracePts, muAve, 1, penaltyCoeff[i], 1,
//                         penaltyCoeff[i], 1);
//         }
//         else
//         {
//             Vmath::Vmul(nTracePts, tcAve, 1, penaltyCoeff[i], 1,
//                         penaltyCoeff[i], 1);
//         }
//     }
// }

void ElectrostaticTurbulence::AddNeutralSources(
    const Array<OneD, Array<OneD, NekDouble>> &in_arr,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
}

/**
 * @brief Populate rhs array ( @p outarray )
 *
 * @param inarray physical values of all fields
 * @param[out] outarray output array (RHSs of time integration equations)
 */
void ElectrostaticTurbulence::DoOdeImplicitRhs(
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
void ElectrostaticTurbulence::DoOdeRhsCoeff(
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
void ElectrostaticTurbulence::DoAdvectionCoeff(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = inarray.size() - 1;

    std::dynamic_pointer_cast<SU::AdvectionWeakDG>(m_advection)
        ->AdvectCoeffs(nvariables, m_indfields, this->v_ExB, inarray, outarray,
                       time, pFwd, pBwd);
}

void ElectrostaticTurbulence::DoDiffusionCoeff(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    const Array<OneD, const Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, const Array<OneD, NekDouble>> &pBwd)
{
    size_t nvariables = inarray.size();
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
    Array<OneD, Array<OneD, NekDouble>> inarrayDiff{nvariables - 1};
    Array<OneD, Array<OneD, NekDouble>> inFwd{nvariables - 1};
    Array<OneD, Array<OneD, NekDouble>> inBwd{nvariables - 1};

    for (int i = 0; i < nvariables; ++i)
    {
        inarrayDiff[i] = Array<OneD, NekDouble>{this->n_pts};
        inFwd[i]       = Array<OneD, NekDouble>{nTracePts};
        inBwd[i]       = Array<OneD, NekDouble>{nTracePts};
    }

    // Extract temperature
    m_varConv->GetElectronTemperature(inarray, inarrayDiff[ee_idx]);
    for (const auto &[s, v] : this->GetIons())
    {
        int ei_idx = v.fields.at(field_to_index["e"]);
        m_varConv->GetIonTemperature(s, v.mass, inarray, inarrayDiff[ei_idx]);
    }

    // Repeat calculation for trace space
    if (pFwd == NullNekDoubleArrayOfArray || pBwd == NullNekDoubleArrayOfArray)
    {
        inFwd = NullNekDoubleArrayOfArray;
        inBwd = NullNekDoubleArrayOfArray;
    }
    else
    {
        m_varConv->GetElectronTemperature(pFwd, inFwd[ee_idx]);
        for (const auto &[s, v] : this->GetIons())
        {
            int ei_idx = v.fields.at(field_to_index["e"]);
            m_varConv->GetIonTemperature(s, v.mass, pFwd, inFwd[ei_idx]);
            m_varConv->GetIonTemperature(s, v.mass, pBwd, inBwd[ei_idx]);
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
void ElectrostaticTurbulence::DoParticlesCoeff(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int ncoeff = GetNcoeffs();
    Array<OneD, NekDouble> tmp(ncoeff, 0.0);

    // Add contribution to electron energy
    m_indfields[ee_idx]->FwdTrans(this->src_fields[0]->GetPhys(), tmp);
    Vmath::Vadd(this->n_pts, outarray[ee_idx], 1, tmp, 1, outarray[ee_idx], 1);

    for (const auto &[s, v] : this->GetIons())
    {
        int ni_idx = v.fields.at(field_to_index["n"]);
        int vi_idx = v.fields.at(field_to_index["v"]);
        int ei_idx = v.fields.at(field_to_index["e"]);
        //  Add contribution to ion density
        m_indfields[ni_idx]->FwdTrans(
            this->src_fields[ni_src_idx[s]]->GetPhys(), tmp);
        Vmath::Vadd(this->n_pts, outarray[ni_idx], 1, tmp, 1, outarray[ni_idx],
                    1);

        // Add contribution to ion energy
        m_indfields[ei_idx]->FwdTrans(
            this->src_fields[ei_src_idx[s]]->GetPhys(), tmp);
        Vmath::Vadd(this->n_pts, outarray[ei_idx], 1, tmp, 1, outarray[ei_idx],
                    1);

        // Add number density source contribution to ion energy
        Array<OneD, NekDouble> dynamic_energy(this->n_pts);
        m_varConv->GetIonDynamicEnergy(s, v.mass, inarray, dynamic_energy);
        Vmath::Vmul(this->n_pts, dynamic_energy, 1,
                    this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                    dynamic_energy, 1);
        m_fields[ei_idx]->FwdTrans(dynamic_energy, tmp);
        Vmath::Vadd(this->n_pts, outarray[ei_idx], 1, tmp, 1, outarray[ei_idx],
                    1);

        Vmath::Zero(this->n_pts, dynamic_energy, 1);
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(this->n_pts, this->b_unit[d], 1,
                         this->src_fields[vi_src_idx[s] + d]->GetPhys(), 1,
                         dynamic_energy, 1, dynamic_energy, 1);
        }
        m_fields[vi_idx]->FwdTrans(dynamic_energy, tmp);
        Vmath::Vadd(this->n_pts, outarray[vi_idx], 1, tmp, 1, outarray[vi_idx],
                    1);
    }
}

/**
 * @brief After reading ICs, calculate phi and grad(phi)
 */
void ElectrostaticTurbulence::v_SetInitialConditions(NekDouble init_time,
                                                     bool dump_ICs,
                                                     const int domain)
{
    PlasmaSystem::v_SetInitialConditions(init_time, dump_ICs, domain);
    CalcInitPhi();
    Checkpoint_Output(0);
}

void ElectrostaticTurbulence::SetBoundaryConditions(
    Array<OneD, Array<OneD, NekDouble>> &physarray, NekDouble time)
{
    PlasmaSystem::SetBoundaryConditions(physarray, time);
    // for (auto &bc : m_bndConds)
    // {
    //     bc->Apply(Fwd, physarray, time);
    // }
    // ApplyOmegaBC(physarray, time);
}

void ElectrostaticTurbulence::load_params()
{
    PlasmaSystem::load_params();

    std::string boussinesq_str;
    m_session->LoadSolverInfo("Boussinesq Approximation", boussinesq_str, "On");
    this->m_boussinesq = (boussinesq_str == "On");
}

void ElectrostaticTurbulence::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    PlasmaSystem::v_ExtraFldOutput(fieldcoeffs, variables);
    const int nCoeffs = m_fields[0]->GetNcoeffs();

    m_fields[0]->FwdTransLocalElmt(this->phi->GetPhys(), fieldcoeffs[4]);

    variables.emplace_back("J_par");
    Array<OneD, NekDouble> SrcFwd(nCoeffs);
    m_fields[0]->FwdTransLocalElmt(this->j_par, SrcFwd);
    fieldcoeffs.emplace_back(SrcFwd);

    for (int d = 0; d < 3; ++d)
    {
        variables.emplace_back("DIAV_" + std::to_string(d));
        Array<OneD, NekDouble> SrcFwd1(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->dia_v[0][d], SrcFwd1);
        fieldcoeffs.emplace_back(SrcFwd1);
    }

    int f = 0;
    for (const auto &[k, v] : this->GetSpecies())
    {
        variables.emplace_back("T_" + v.name);
        Array<OneD, NekDouble> Fwd(nCoeffs);
        this->m_temps[f]->FwdTransLocalElmt(this->m_temps[f]->GetPhys(), Fwd);
        fieldcoeffs.emplace_back(Fwd);
        f++;
    }
    variables.emplace_back("T_e");
    Array<OneD, NekDouble> Fwd(nCoeffs);
    this->m_temps[f]->FwdTransLocalElmt(this->m_temps[f]->GetPhys(), Fwd);
    fieldcoeffs.emplace_back(Fwd);

    if (this->particles_enabled)
    {
        int cnt = 0;
        for (auto &[k, v] : this->particle_sys->get_species())
        {
            variables.emplace_back(k + "_SOURCE_DENSITY");
            Array<OneD, NekDouble> SrcFwd1(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[cnt++]->GetPhys(),
                                           SrcFwd1);
            fieldcoeffs.emplace_back(SrcFwd1);

            for (int d = 0; d < this->m_spacedim; ++d)
            {
                variables.emplace_back(k + "_SOURCE_MOMENTUM" +
                                       std::to_string(d));
                Array<OneD, NekDouble> SrcFwd1(nCoeffs);
                m_fields[0]->FwdTransLocalElmt(
                    this->src_fields[cnt++]->GetPhys(), SrcFwd1);
                fieldcoeffs.emplace_back(SrcFwd1);
            }

            variables.emplace_back(k + "_SOURCE_ENERGY");
            Array<OneD, NekDouble> SrcFwd2(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[cnt++]->GetPhys(),
                                           SrcFwd2);
            fieldcoeffs.emplace_back(SrcFwd2);
        }
    }
}
} // namespace PENKNIFE