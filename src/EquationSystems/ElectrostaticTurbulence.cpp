#include "ElectrostaticTurbulence.hpp"
#include <SolverUtils/Advection/AdvectionWeakDG.h>

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
    : TokamakSystem(session, graph)
{
    this->n_indep_fields       = 2;
    this->n_fields_per_species = 3;
}

void ElectrostaticTurbulence::v_InitObject(bool DeclareFields)
{
    TokamakSystem::v_InitObject(DeclareFields);
    m_varConv = MemoryManager<VariableConverter>::AllocateSharedPtr(
        m_session, m_spacedim, m_graph);

    // Since we are starting from a setup where each field is defined to be a
    // discontinuous field (and thus support DG), the first thing we do is to
    // recreate the phi field so that it is continuous, in order to support the
    // Poisson solve. Note that you can still perform a Poisson solve using a
    // discontinuous field, which is done via the hybridisable discontinuous
    // Galerkin (HDG) approach.
    this->phi = MemoryManager<MR::ContField>::AllocateSharedPtr(
        m_session, m_graph, "phi", true, true);

    std::string diffName;
    m_session->LoadSolverInfo("DiffusionType", diffName, "LDG");
    m_diffusion =
        SolverUtils::GetDiffusionFactory().CreateInstance(diffName, diffName);
    m_diffusion->SetFluxVector(&ElectrostaticTurbulence::GetFluxVectorDiff,
                               this);
    m_diffusion->InitObject(m_session, m_indfields);

    // Create storage for velocities
    int npts = GetNpoints();

    auto tmp_fields = Array<OneD, MultiRegions::ExpListSharedPtr>(
        m_fields.size() + n_species * 3);

    for (int i = 0; i < m_fields.size(); ++i)
    {
        tmp_fields[i] = m_fields[i];
    }

    // ExB velocity
    this->v_ExB = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);

    // Parallel velocities
    // this->v_e_par = Array<OneD, NekDouble>(npts, 0.0);
    this->v_i_par = std::vector<Array<OneD, NekDouble>>(n_species);

    // Drift velocities
    this->v_de = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    this->v_di = std::vector<Array<OneD, Array<OneD, NekDouble>>>(n_species);

    // Per-field advection velocities
    this->adv_vel =
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(2 + 3 * n_species);

    for (int i = 0; i < this->adv_vel.size(); ++i)
    {
        this->adv_vel[i] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    }

    for (int d = 0; d < m_spacedim; ++d)
    {
        this->v_ExB[d] = Array<OneD, NekDouble>(npts, 0.0);
        this->v_de[d]  = Array<OneD, NekDouble>(npts, 0.0);
        // this->adv_vel[ne_idx][d] = Array<OneD, NekDouble>(npts, 0.0);
    }

    pe_idx    = 3 * this->n_species;
    omega_idx = 3 * this->n_species + 1;
    for (const auto &[s, v] : this->particle_sys->get_species())
    {
        ni_idx.push_back(3 * s);
        vi_idx.push_back(1 + 3 * s);
        pi_idx.push_back(2 + 3 * s);

        this->v_i_par[s] = Array<OneD, NekDouble>(npts, 0.0);
        this->v_di[s]    = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int d = 0; d < m_graph->GetSpaceDimension(); ++d)
        {
            this->v_di[s][d]            = Array<OneD, NekDouble>(npts, 0.0);
            this->adv_vel[ni_idx[s]][d] = Array<OneD, NekDouble>(npts, 0.0);
            this->adv_vel[vi_idx[s]][d] = Array<OneD, NekDouble>(npts, 0.0);
            this->adv_vel[pi_idx[s]][d] = Array<OneD, NekDouble>(npts, 0.0);
        }
        m_varConv->charge[s] = v.charge;
        m_varConv->mass[s]   = v.mass;
    }

    m_fields = tmp_fields;

    m_varConv->ni_idx = ni_idx;
    m_varConv->vi_idx = vi_idx;
    m_varConv->pi_idx = pi_idx;

    if (m_fields[0]->GetTrace())
    {
        auto nTrace = GetTraceNpoints();
        this->trace_vel_norm =
            Array<OneD, Array<OneD, NekDouble>>(adv_vel.size());
        for (int i = 0; i < this->trace_vel_norm.size(); ++i)
        {
            this->trace_vel_norm[i] = Array<OneD, NekDouble>(nTrace, 0.0);
        }
    }

    // Create Riemann solver and set normal velocity
    // callback functions
    this->riemann_solver = SU::GetRiemannSolverFactory().CreateInstance(
        this->riemann_solver_type, m_session);
    this->riemann_solver->SetVector(
        "Vn", &ElectrostaticTurbulence::GetAdvVelNorm, this);

    // Setup advection object
    m_advection = SU::GetAdvectionFactory().CreateInstance(this->adv_type,
                                                           this->adv_type);
    m_advection->SetFluxVector(&ElectrostaticTurbulence::GetFluxVector, this);
    m_advection->SetRiemannSolver(this->riemann_solver);
    m_advection->InitObject(m_session, m_indfields);

    m_ode.DefineOdeRhs(&ElectrostaticTurbulence::DoOdeRhs, this);

    if (this->particles_enabled)
    {
        std::vector<Sym<REAL>> src_syms;
        std::vector<int> src_components;

        int s = 0;
        for (auto &[k, v] : this->particle_sys->get_species())
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
            ++s;
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
    int npts      = GetNpoints();
    int nTracePts = GetTraceTotPoints();

    if (this->m_explicitAdvection)
    {
        zero_array_of_arrays(outarray);
    }
    int nvariables = inarray.size();

    Array<OneD, Array<OneD, NekDouble>> tmp(nvariables);
    tmp = inarray;

    Array<OneD, NekDouble> ne(npts);
    m_varConv->GetElectronDensity(inarray, ne);

    // Store forwards/backwards space along trace space
    Array<OneD, Array<OneD, NekDouble>> Fwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> Bwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        Fwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        Bwd[i] = Array<OneD, NekDouble>(nTracePts, 0.0);
        m_indfields[i]->GetFwdBwdTracePhys(tmp[i], Fwd[i], Bwd[i]);
    }

    // Helmholtz Solve for electrostatic potential
    SolvePhi(tmp, ne);
    // Calculate E
    ComputeE();
    // Calculate ExB, parallel and diamagnetic velocities
    CalcVelocities(tmp, ne);

    // Perform advection
    DoAdvection(tmp, outarray, time, Fwd, Bwd);

    // add extra terms to electron momentum
    Array<OneD, Array<OneD, NekDouble>> gradp(m_spacedim);
    Array<OneD, NekDouble> extra(npts);
    // m_fields[3]->PhysDeriv(tmp[pe_idx], gradp[0], gradp[1], gradp[2]);
    //
    // for (int d = 0; d < m_spacedim; ++d)
    //{
    //    Vmath::Vvtvp(npts, b_unit[d], 1, gradp[d], 1, extra, 1, extra, 1);
    //}
    // Vmath::Smul(npts, 2.0 / 3.0, extra, 1, extra, 1);
    // Vmath::Vadd(npts, extra, 1, outarray[ve_idx], 1, outarray[ve_idx], 1);
    // Vmath::Zero(npts, extra, 1);
    // for (int d = 0; d < m_spacedim; ++d)
    //{
    //    Vmath::Vvtvp(npts, b_unit[d], 1, this->E[d]->GetPhys(), 1, extra, 1,
    //                 extra, 1);
    //}
    // Vmath::Vmul(npts, outarray[ne_idx], 1, extra, 1, extra, 1);
    // Vmath::Vadd(npts, extra, 1, outarray[ve_idx], 1, outarray[ve_idx], 1);

    // Add extra terms to electron energy
    Array<OneD, Array<OneD, NekDouble>> vel(m_spacedim);

    for (int d = 0; d < m_spacedim; ++d)
    {
        Vmath::Vvtvp(npts, this->b_unit[d], 1, this->v_e_par, 1, this->v_ExB[d],
                     1, vel[d], 1);
        m_fields[0]->PhysDeriv(d, vel[d], gradp[d]);
        Vmath::Vadd(npts, gradp[d], 1, extra, 1, extra, 1);
    }
    Vmath::Smul(npts, 2.0 / 3.0, extra, 1, extra, 1);
    Vmath::Vadd(npts, extra, 1, outarray[pe_idx], 1, outarray[pe_idx], 1);
    Vmath::Zero(npts, extra, 1);

    for (const auto &[s, v] : this->particle_sys->get_species())
    {
        // Add extra terms to ion momentum
        m_fields[0]->PhysDeriv(tmp[pi_idx[s]], gradp[0], gradp[1], gradp[2]);

        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(npts, b_unit[d], 1, gradp[d], 1, extra, 1, extra, 1);
        }
        Vmath::Smul(npts, 2.0 / 3.0, extra, 1, extra, 1);
        Vmath::Vadd(npts, extra, 1, outarray[vi_idx[s]], 1, outarray[vi_idx[s]],
                    1);
        Vmath::Zero(npts, extra, 1);
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(npts, b_unit[d], 1, this->E[d]->GetPhys(), 1, extra, 1,
                         extra, 1);
        }
        Vmath::Vmul(npts, outarray[0], 1, extra, 1, extra, 1);
        Vmath::Smul(npts, v.charge, extra, 1, extra, 1);

        Vmath::Vadd(npts, extra, 1, outarray[vi_idx[s]], 1, outarray[vi_idx[s]],
                    1);
        // Add extra terms to ion energy
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(npts, this->b_unit[d], 1, this->v_i_par[s], 1,
                         this->v_ExB[d], 1, vel[d], 1);
            m_fields[0]->PhysDeriv(d, vel[d], gradp[d]);
            Vmath::Vadd(npts, gradp[d], 1, extra, 1, extra, 1);
        }
        Vmath::Smul(npts, 2.0 / 3.0, extra, 1, extra, 1);
        Vmath::Vadd(npts, extra, 1, outarray[pi_idx[s]], 1, outarray[pi_idx[s]],
                    1);
    }

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Neg(npts, outarray[i], 1);
    }
    CalcDiffTensor();
    CalcKappaTensor();
    // Perform Diffusion
    DoDiffusion(tmp, outarray, Fwd, Bwd);

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
void ElectrostaticTurbulence::DoAdvection(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    int nvariables = inarray.size();

    m_advection->Advect(nvariables, m_indfields, this->v_ExB, inarray, outarray,
                        time, pFwd, pBwd);
}

/**
 * @brief Compute the advection terms for the right-hand side
 */
void ElectrostaticTurbulence::DoParticles(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    int npts = GetNpoints();

    // Add contribution to electron energy
    Vmath::Vadd(npts, outarray[pe_idx], 1, this->src_fields[0]->GetPhys(), 1,
                outarray[pe_idx], 1);

    for (const auto &[s, v] : this->particle_sys->get_species())
    {
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
        m_varConv->GetIonDynamicEnergy(s, v.mass, inarray, dynamic_energy);
        Vmath::Vvtvp(npts, dynamic_energy, 1,
                     this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                     outarray[pi_idx[s]], 1, outarray[pi_idx[s]], 1);

        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(npts, this->b_unit[d], 1,
                         this->src_fields[vi_src_idx[s] + d]->GetPhys(), 1,
                         outarray[vi_idx[s]], 1, outarray[vi_idx[s]], 1);
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
    int npts = GetNpoints();

    StdRegions::ConstFactorMap factors;
    // Helmholtz => Poisson (lambda = 0)
    factors[StdRegions::eFactorLambda] = 0.0;

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Array<OneD, NekDouble> D(npts, 0.0);
            for (int k = 0; k < npts; k++)
            {
                D[k] = -b_unit[i][k] * b_unit[j][k];
                if (i == j)
                {
                    D[k] += 1;
                }
            }
            if (!m_boussinesq)
            {
                Vmath::Vmul(npts, D, 1, ne, 1, D, 1);
            }
            m_phi_varcoeff[vc[i][j]] = D;
        }
    }

    // Solve for phi. Output of this routine is in coefficient (spectral)
    // space, so backwards transform to physical space since we'll need that
    // for the advection step & computing drift velocity.
    this->phi->HelmSolve(inarray[omega_idx], this->phi->UpdateCoeffs(), factors,
                         m_phi_varcoeff);
    this->phi->BwdTrans(this->phi->GetCoeffs(), this->phi->UpdatePhys());
}

/**
 * @brief Calculates initial potential and gradient
 */
void ElectrostaticTurbulence::CalcInitPhi()
{
    int npts = GetNpoints();
    Array<OneD, Array<OneD, NekDouble>> inarray(m_fields.size());
    Array<OneD, NekDouble> ne(npts);
    for (int i = 0; i < m_fields.size(); i++)
    {
        inarray[i] = m_fields[i]->GetPhys();
    }
    m_varConv->GetElectronDensity(inarray, ne);
    SolvePhi(inarray, ne);
    ComputeE();
}

/**
 * @brief Compute the gradient of phi for evaluation at the particle positions.
 */
void ElectrostaticTurbulence::ComputeE()
{
    this->phi->PhysDeriv(this->phi->GetPhys(), this->E[0]->UpdatePhys(),
                         this->E[1]->UpdatePhys(), this->E[2]->UpdatePhys());
    int npts = GetNpoints();

    Vmath::Smul(npts, 1.0, this->E[0]->GetPhys(), 1, this->E[0]->UpdatePhys(),
                1);
    Vmath::Smul(npts, 1.0, this->E[1]->GetPhys(), 1, this->E[1]->UpdatePhys(),
                1);
    Vmath::Smul(npts, 1.0, this->E[2]->GetPhys(), 1, this->E[2]->UpdatePhys(),
                1);

    this->E[0]->FwdTrans(this->E[0]->GetPhys(), this->E[0]->UpdateCoeffs());
    this->E[1]->FwdTrans(this->E[1]->GetPhys(), this->E[1]->UpdateCoeffs());
    this->E[2]->FwdTrans(this->E[2]->GetPhys(), this->E[2]->UpdateCoeffs());
}

void ElectrostaticTurbulence::CalcVelocities(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, NekDouble> &ne)
{
    int npts = GetNpoints();

    // Calculate ExB velocity
    Vmath::Vvtvvtm(npts, this->E[1]->GetPhys(), 1, this->B[2]->GetPhys(), 1,
                   this->E[2]->GetPhys(), 1, this->B[1]->GetPhys(), 1,
                   this->v_ExB[0], 1);
    Vmath::Vdiv(npts, this->v_ExB[0], 1, this->mag_B, 1, this->v_ExB[0], 1);

    Vmath::Vvtvvtm(npts, this->E[2]->GetPhys(), 1, this->B[0]->GetPhys(), 1,
                   this->E[0]->GetPhys(), 1, this->B[2]->GetPhys(), 1,
                   this->v_ExB[1], 1);
    Vmath::Vdiv(npts, this->v_ExB[1], 1, this->mag_B, 1, this->v_ExB[1], 1);

    if (this->m_graph->GetMeshDimension() == 3)
    {
        Vmath::Vvtvvtm(npts, this->E[0]->GetPhys(), 1, this->B[1]->GetPhys(), 1,
                       this->E[1]->GetPhys(), 1, this->B[0]->GetPhys(), 1,
                       this->v_ExB[2], 1);
        Vmath::Vdiv(npts, this->v_ExB[2], 1, this->mag_B, 1, this->v_ExB[2], 1);
    }

    // Calculate Electron parallel velocity
    Vmath::Vdiv(npts, inarray[2], 1, inarray[0], 1, this->v_e_par, 1);

    // Calculate Ion parallel velocities
    for (const auto &[s, v] : this->particle_sys->get_species())
    {
        Vmath::Vdiv(npts, inarray[vi_idx[s]], 1, inarray[ni_idx[s]], 1,
                    this->v_i_par[s], 1);
        Vmath::Smul(npts, 1.0 / v.mass, this->v_i_par[s], 1, this->v_i_par[s],
                    1);
    }

    // Calculate electron diagmagnetic velocity
    Array<OneD, Array<OneD, NekDouble>> gradp(m_spacedim);
    m_fields[3]->PhysDeriv(inarray[pe_idx], gradp[0], gradp[1], gradp[2]);

    Vmath::Vvtvvtm(npts, gradp[1], 1, this->B[2]->GetPhys(), 1, gradp[2], 1,
                   this->B[1]->GetPhys(), 1, this->v_de[0], 1);
    Vmath::Vdiv(npts, this->v_de[0], 1, this->mag_B, 1, this->v_de[0], 1);
    Vmath::Vdiv(npts, this->v_de[0], 1, ne, 1, this->v_de[0], 1);
    Vmath::Smul(npts, 2.0 / 3.0, this->v_de[0], 1, this->v_de[0], 1);

    Vmath::Vvtvvtm(npts, gradp[2], 1, this->B[0]->GetPhys(), 1, gradp[0], 1,
                   this->B[2]->GetPhys(), 1, this->v_de[1], 1);
    Vmath::Vdiv(npts, this->v_de[1], 1, this->mag_B, 1, this->v_de[1], 1);
    Vmath::Vdiv(npts, this->v_de[1], 1, ne, 1, this->v_de[1], 1);
    Vmath::Smul(npts, 2.0 / 3.0, this->v_de[1], 1, this->v_de[1], 1);

    if (this->m_graph->GetMeshDimension() == 3)
    {
        Vmath::Vvtvvtm(npts, gradp[0], 1, this->B[1]->GetPhys(), 1, gradp[1], 1,
                       this->B[0]->GetPhys(), 1, this->v_de[2], 1);
        Vmath::Vdiv(npts, this->v_de[2], 1, this->mag_B, 1, this->v_de[2], 1);
        Vmath::Vdiv(npts, this->v_de[2], 1, ne, 1, this->v_de[2], 1);
        Vmath::Smul(npts, 2.0 / 3.0, this->v_de[2], 1, this->v_de[2], 1);
    }

    for (int d = 0; d < m_spacedim; ++d)
    {
        Vmath::Svtvp(npts, 5.0 / 3.0, this->v_de[d], 1, this->v_ExB[d], 1,
                     this->adv_vel[pe_idx][d], 1);
    }

    // Calculate ion diamagnetic velocity
    for (const auto &[s, v] : this->particle_sys->get_species())
    {
        m_fields[pi_idx[s]]->PhysDeriv(inarray[pi_idx[s]], gradp[0], gradp[1],
                                       gradp[2]);

        Vmath::Vvtvvtm(npts, gradp[1], 1, this->B[2]->GetPhys(), 1, gradp[2], 1,
                       this->B[1]->GetPhys(), 1, this->v_di[s][0], 1);
        Vmath::Vdiv(npts, this->v_di[s][0], 1, this->mag_B, 1, this->v_di[s][0],
                    1);
        Vmath::Vdiv(npts, this->v_di[s][0], 1, inarray[0], 1, this->v_di[s][0],
                    1);
        Vmath::Smul(npts, 2.0 / (3.0 * v.charge), this->v_di[s][0], 1,
                    this->v_di[s][0], 1);
        Vmath::Vvtvvtm(npts, gradp[2], 1, this->B[0]->GetPhys(), 1, gradp[0], 1,
                       this->B[2]->GetPhys(), 1, this->v_di[s][1], 1);
        Vmath::Vdiv(npts, this->v_di[s][1], 1, this->mag_B, 1, this->v_di[s][1],
                    1);
        Vmath::Vdiv(npts, this->v_di[s][1], 1, ne, 1, this->v_di[s][1], 1);
        Vmath::Smul(npts, 2.0 / (3.0 * v.charge), this->v_di[s][1], 1,
                    this->v_di[s][1], 1);

        if (this->m_graph->GetMeshDimension() == 3)
        {
            Vmath::Vvtvvtm(npts, gradp[0], 1, this->B[1]->GetPhys(), 1,
                           gradp[1], 1, this->B[0]->GetPhys(), 1,
                           this->v_di[s][2], 1);
            Vmath::Vdiv(npts, this->v_di[s][2], 1, this->mag_B, 1,
                        this->v_di[s][2], 1);
            Vmath::Vdiv(npts, this->v_di[s][2], 1, ne, 1, this->v_di[s][2], 1);
            Vmath::Smul(npts, 2.0 / (3.0 * v.charge), this->v_di[s][2], 1,
                        this->v_di[s][2], 1);
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
Array<OneD, Array<OneD, NekDouble>> &ElectrostaticTurbulence::GetAdvVelNorm()
{
    // Number of trace (interface) points
    int num_trace_pts = GetTraceNpoints();
    // Auxiliary variable to compute normal velocities
    Array<OneD, NekDouble> tmp(num_trace_pts);

    // Compute advection vel dot trace normals and store
    for (int j = 0; j < this->adv_vel.size();)
    {
        // Ensure output array is zeroed
        Vmath::Zero(num_trace_pts, this->trace_vel_norm[j], 1);
        for (int d = 0; d < this->adv_vel[j].size(); ++d)
        {
            m_fields[0]->ExtractTracePhys(this->adv_vel[j][d], tmp);
            Vmath::Vvtvp(num_trace_pts, m_traceNormals[d], 1, tmp, 1,
                         this->trace_vel_norm[j], 1, this->trace_vel_norm[j],
                         1);
        }
    }
    return trace_vel_norm;
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
    ASSERTL1(
        fluxes[0].size() == this->adv_vel.size(),
        "Dimension of flux array and advection velocity array do not match");
    int npts = field_vals[0].size();

    for (int d = 0; d < m_spacedim; ++d)
    {
        // Vorticity Flux
        Vmath::Vmul(npts, field_vals[omega_idx], 1, this->adv_vel[omega_idx][d],
                    1, fluxes[omega_idx][d], 1);

        // Electron Energy Flux
        Vmath::Vmul(npts, field_vals[pe_idx], 1, this->adv_vel[pe_idx][d], 1,
                    fluxes[pe_idx][d], 1);
    }

    for (const auto &[s, v] : this->particle_sys->get_species())
    {
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
    }
}

void ElectrostaticTurbulence::DoDiffusion(
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

    // Get primitive variables [u,v,w,T]
    Array<OneD, Array<OneD, NekDouble>> inarrayDiff(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inFwd(nvariables);
    Array<OneD, Array<OneD, NekDouble>> inBwd(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        inarrayDiff[i] = Array<OneD, NekDouble>(npointsIn);
        inFwd[i]       = Array<OneD, NekDouble>(nTracePts);
        inBwd[i]       = Array<OneD, NekDouble>(nTracePts);
    }

    // Extract temperature
    m_varConv->GetElectronTemperature(inarray, inarrayDiff[pe_idx]);

    for (const auto &[s, v] : this->particle_sys->get_species())
    {
        m_varConv->GetIonTemperature(s, v.mass, inarray,
                                     inarrayDiff[pi_idx[s]]);
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
        for (const auto &[s, v] : this->particle_sys->get_species())
        {
            m_varConv->GetIonTemperature(s, v.mass, pFwd, inFwd[pi_idx[s]]);
            m_varConv->GetIonTemperature(s, v.mass, pBwd, inBwd[pi_idx[s]]);
        }
    }

    m_diffusion->Diffuse(nvariables, m_fields, inarrayDiff, outarrayDiff, inFwd,
                         inBwd);

    for (int i = 0; i < nvariables; ++i)
    {
        Vmath::Vadd(npointsOut, outarrayDiff[i], 1, outarray[i], 1, outarray[i],
                    1);
    }
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
 * @brief Construct the flux vector for the anisotropic diffusion problem.
 */
void ElectrostaticTurbulence::GetFluxVectorDiff(
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &qfield,
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes)
{
    unsigned int nPts = inarray[0].size();

    for (int j = 0; j < m_spacedim; ++j)
    {
        // Calc diffusion of n with D tensor
        // Vmath::Vmul(nPts, m_D[vc[j][0]].GetValue(), 1, qfield[0][ne_idx], 1,
        //            fluxes[j][ne_idx], 1);
        Vmath::Vmul(nPts, m_kappa[vc[j][0]].GetValue(), 1, qfield[0][pe_idx], 1,
                    fluxes[j][pe_idx], 1);
        for (int k = 1; k < m_spacedim; ++k)
        {
            // Vmath::Vvtvp(nPts, m_D[vc[j][k]].GetValue(), 1,
            // qfield[k][ne_idx],
            //              1, fluxes[j][ne_idx], 1, fluxes[j][ne_idx], 1);
            Vmath::Vvtvp(nPts, m_kappa[vc[j][k]].GetValue(), 1,
                         qfield[k][pe_idx], 1, fluxes[j][pe_idx], 1,
                         fluxes[j][pe_idx], 1);
        }
    }
    for (const auto &[s, v] : this->particle_sys->get_species())
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
    }
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
    int nvariables = inarray.size();

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
    Array<OneD, Array<OneD, NekDouble>> inarrayDiff{nvariables - 1};
    Array<OneD, Array<OneD, NekDouble>> inFwd{nvariables - 1};
    Array<OneD, Array<OneD, NekDouble>> inBwd{nvariables - 1};

    for (int i = 0; i < nvariables - 1; ++i)
    {
        inarrayDiff[i] = Array<OneD, NekDouble>{npoints};
        inFwd[i]       = Array<OneD, NekDouble>{nTracePts};
        inBwd[i]       = Array<OneD, NekDouble>{nTracePts};
    }

    // Extract temperature
    m_varConv->GetElectronTemperature(inarray, inarrayDiff[pe_idx]);

    for (const auto &[s, v] : this->particle_sys->get_species())
    {
        m_varConv->GetIonTemperature(s, v.mass, inarray,
                                     inarrayDiff[pi_idx[s]]);
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
        for (const auto &[s, v] : this->particle_sys->get_species())
        {
            m_varConv->GetIonTemperature(s, v.mass, pFwd, inFwd[pi_idx[s]]);
            m_varConv->GetIonTemperature(s, v.mass, pBwd, inBwd[pi_idx[s]]);
        }
    }

    // Diffusion term in physical rhs form
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
    int npts   = GetNpoints();
    int ncoeff = GetNcoeffs();
    Array<OneD, NekDouble> tmp(ncoeff, 0.0);

    // Add contribution to electron energy
    m_fields[0]->FwdTrans(this->src_fields[0]->GetPhys(), tmp);
    Vmath::Vadd(npts, outarray[pe_idx], 1, tmp, 1, outarray[pe_idx], 1);

    for (const auto &[s, v] : this->particle_sys->get_species())
    {
        //  Add contribution to ion density
        m_fields[0]->FwdTrans(this->src_fields[ni_src_idx[s]]->GetPhys(), tmp);
        Vmath::Vadd(npts, outarray[ni_idx[s]], 1, tmp, 1, outarray[ni_idx[s]],
                    1);

        // Add contribution to ion energy
        m_fields[0]->FwdTrans(this->src_fields[pi_src_idx[s]]->GetPhys(), tmp);
        Vmath::Vadd(npts, outarray[pi_idx[s]], 1, tmp, 1, outarray[pi_idx[s]],
                    1);

        // Add number density source contribution to ion energy
        Array<OneD, NekDouble> dynamic_energy(npts);
        m_varConv->GetIonDynamicEnergy(s, v.mass, inarray, dynamic_energy);
        Vmath::Vmul(npts, dynamic_energy, 1,
                    this->src_fields[ni_src_idx[s]]->GetPhys(), 1,
                    dynamic_energy, 1);
        m_fields[0]->FwdTrans(dynamic_energy, tmp);
        Vmath::Vadd(npts, outarray[pi_idx[s]], 1, tmp, 1, outarray[pi_idx[s]],
                    1);

        Vmath::Zero(npts, dynamic_energy, 1);
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(npts, this->b_unit[d], 1,
                         this->src_fields[vi_src_idx[s] + d]->GetPhys(), 1,
                         dynamic_energy, 1, dynamic_energy, 1);
        }
        m_fields[0]->FwdTrans(dynamic_energy, tmp);
        Vmath::Vadd(npts, outarray[vi_idx[s]], 1, tmp, 1, outarray[vi_idx[s]],
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
    std::string boussinesq_str;
    m_session->LoadSolverInfo("Boussinesq Approximation", boussinesq_str,
                              "Off");
    this->m_boussinesq = (boussinesq_str == "On");
}

void ElectrostaticTurbulence::v_ExtraFldOutput(
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

        for (const auto &[s, v] : this->particle_sys->get_species())
        {
            m_varConv->GetIonPressure(s, v.mass, tmp, pressure);
            m_varConv->GetIonTemperature(s, v.mass, tmp, temperature);

            Array<OneD, NekDouble> pFwd(nCoeffs), TFwd(nCoeffs);

            m_fields[0]->FwdTransLocalElmt(pressure, pFwd);
            m_fields[0]->FwdTransLocalElmt(temperature, TFwd);

            variables.push_back("pi");
            variables.push_back("Ti");

            fieldcoeffs.push_back(pFwd);
            fieldcoeffs.push_back(TFwd);
        }

        if (this->particles_enabled)
        {
            variables.push_back("ELECTRON_SOURCE_DENSITY");
            Array<OneD, NekDouble> SrcFwd0(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[0]->GetPhys(),
                                           SrcFwd0);
            fieldcoeffs.push_back(SrcFwd0);
            variables.push_back("ELECTRON_SOURCE_ENERGY");
            Array<OneD, NekDouble> SrcFwd1(nCoeffs);
            m_fields[0]->FwdTransLocalElmt(this->src_fields[1]->GetPhys(),
                                           SrcFwd1);
            fieldcoeffs.push_back(SrcFwd1);

            for (auto s : this->particle_sys->get_species())
            {
                variables.push_back(s.second.name + "_SOURCE_ENERGY");
                Array<OneD, NekDouble> SrcFwd2(nCoeffs);
                m_fields[0]->FwdTransLocalElmt(
                    this->src_fields[2 + s.first]->GetPhys(), SrcFwd2);
                fieldcoeffs.push_back(SrcFwd2);
            }
        }
    }
}
} // namespace NESO::Solvers::tokamak