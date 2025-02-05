#include "TokamakSystem.hpp"
#include <FieldUtils/Interpolator.h>
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
        "Tokamak equation system. Runs in either a 2D or 3D "
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
    : TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>(session, graph)
{
    //  Determine whether energy,enstrophy recording is enabled (output freq is
    //  set)
    this->energy_enstrophy_recording_enabled =
        session->DefinesParameter("growth_rates_recording_step");
}

std::shared_ptr<ParticleSystem> TokamakSystem::GetParticleSystem()
{
    return this->particle_sys;
}

/**
 * @brief Read the magnetic field from file.
 */
void TokamakSystem::ReadMagneticField()
{
    int d;
    int npoints = m_fields[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble>> B_in(3);
    for (d = 0; d < 3; ++d)
    {
        B_in[d] = Array<OneD, NekDouble>(npoints);
    }
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
            for (d = 0; d < 3; ++d)
            {
                this->B_pol[d] = B_in[d];
            }
        }
        else if (this->m_graph->GetMeshDimension() == 3)
        {
            LU::FunctionType vType =
                m_session->GetFunctionType("MagneticMeanField", "Bx");
            if (vType == LU::eFunctionTypeFile)
            {
                LU::PtsIO ptsIO(m_session->GetComm());
                std::string filename =
                    m_session->GetFunctionFilename("MagneticMeanField", "Bx");
                LU::PtsFieldSharedPtr inPts2D;
                ptsIO.Import(filename, inPts2D);
                // B2D is [x,y,Bx,By,Bz]
                Array<OneD, Array<OneD, NekDouble>> B2D;
                inPts2D->GetPts(B2D);
                unsigned int npoints_2d = inPts2D->GetNpoints();
                // B2D is [x,y,z,Bx,By,Bz]
                Array<OneD, Array<OneD, NekDouble>> B3D(6);
                unsigned int increments = 64;

                for (d = 0; d < 6; ++d)
                {
                    B3D[d] =
                        Array<OneD, NekDouble>(increments * npoints_2d, 0.0);
                    // B_pol[d] = Array<OneD, NekDouble>(npoints, 0.0);
                }
                Array<OneD, NekDouble> Bx(npoints_2d);
                Array<OneD, NekDouble> By(npoints_2d);
                Array<OneD, NekDouble> Bz(npoints_2d);
                Array<OneD, NekDouble> x(npoints_2d);
                Array<OneD, NekDouble> y(npoints_2d);
                Array<OneD, NekDouble> z(npoints_2d);
                for (int i = 0; i < increments; ++i)
                {
                    NekDouble theta = i * 2 * M_PI / increments;
                    for (int j = 0; j < npoints_2d; ++j)
                    {
                        x[j]  = cos(theta) * B2D[0][j];
                        y[j]  = B2D[1][j];
                        z[j]  = sin(theta) * B2D[0][j];
                        Bx[j] = cos(theta) * B2D[2][j] - sin(theta) * B2D[4][j];
                        By[j] = B2D[3][j];
                        Bz[j] = cos(theta) * B2D[4][j] + sin(theta) * B2D[2][j];
                    }

                    Vmath::Vcopy(npoints_2d, &x[0], 1, &B3D[0][i * npoints_2d],
                                 1);
                    Vmath::Vcopy(npoints_2d, &y[0], 1, &B3D[1][i * npoints_2d],
                                 1);
                    Vmath::Vcopy(npoints_2d, &z[0], 1, &B3D[2][i * npoints_2d],
                                 1);
                    Vmath::Vcopy(npoints_2d, &Bx[0], 1, &B3D[3][i * npoints_2d],
                                 1);
                    Vmath::Vcopy(npoints_2d, &By[0], 1, &B3D[4][i * npoints_2d],
                                 1);
                    Vmath::Vcopy(npoints_2d, &Bz[0], 1, &B3D[5][i * npoints_2d],
                                 1);

                    // Vmath::Vcopy(npoints_2d, &B2D[1][0], 1,
                    //              &B_pol[1][i * npoints_2d], 1);
                    // Vmath::Vcopy(npoints_2d, &B2D[0][0], 1,
                    //              &B_pol[0][i * npoints_2d], 1);
                    // Vmath::Vcopy(npoints_2d, &B2D[2][0], 1,
                    //              &B_pol[2][i * npoints_2d], 1);
                }
                LU::PtsFieldSharedPtr inPts =
                    MemoryManager<LU::PtsField>::AllocateSharedPtr(3, B3D);

                Array<OneD, Array<OneD, NekDouble>> pts(6);
                for (int i = 0; i < 6; ++i)
                {
                    pts[i] = Array<OneD, NekDouble>(npoints);
                }
                m_fields[0]->GetCoords(pts[0], pts[1], pts[2]);
                LU::PtsFieldSharedPtr outPts =
                    MemoryManager<LU::PtsField>::AllocateSharedPtr(3, Bstring,
                                                                   pts);
                FieldUtils::Interpolator<std::vector<MR::ExpListSharedPtr>>
                    interp;

                interp =
                    FieldUtils::Interpolator<std::vector<MR::ExpListSharedPtr>>(
                        LU::eShepard);
                interp.CalcWeights(inPts, outPts);
                if (m_session->GetComm()->GetRank() == 0)
                {
                    std::cout << std::endl;
                    if (m_session->DefinesCmdLineArgument("verbose"))
                    {
                        interp.PrintStatistics();
                    }
                }

                interp.Interpolate(inPts, outPts);
                for (d = 0; d < 3; ++d)
                {
                    // outPts->SetPts(d + 3, B_in[d]);
                    B_in[d] = outPts->GetPts(d + 3);
                }
            }
            else if (vType == LU::eFunctionTypeExpression)
            {
                unsigned int nq = m_fields[0]->GetNpoints();

                Array<OneD, NekDouble> x0(nq);
                Array<OneD, NekDouble> x1(nq);
                Array<OneD, NekDouble> x2(nq);

                // Get the coordinates (assuming all fields have the same
                // discretisation)
                m_fields[0]->GetCoords(x0, x1, x2);

                Array<OneD, NekDouble> r(nq);
                Array<OneD, NekDouble> phi(nq);
                for (int q = 0; q < nq; ++q)
                {
                    r[q]   = std::sqrt(x0[q] * x0[q] + x2[q] * x2[q]);
                    phi[q] = std::atan2(x2[q], x0[q]);
                }
                LibUtilities::EquationSharedPtr Bxfunc =
                    m_session->GetFunction("MagneticMeanField", "Bx");
                LibUtilities::EquationSharedPtr Byfunc =
                    m_session->GetFunction("MagneticMeanField", "By");
                LibUtilities::EquationSharedPtr Bzfunc =
                    m_session->GetFunction("MagneticMeanField", "Bz");
                Array<OneD, NekDouble> Br(nq);
                Array<OneD, NekDouble> Bphi(nq);

                Bxfunc->Evaluate(r, x1, phi, /*time=*/0, Br);
                Byfunc->Evaluate(r, x1, phi, /*time=*/0, B_in[1]);
                Bzfunc->Evaluate(r, x1, phi, /*time=*/0, Bphi);

                for (int q = 0; q < nq; ++q)
                {
                    B_in[0][q] = cos(phi[q]) * Br[q] - sin(phi[q]) * Bphi[q];
                    B_in[2][q] = cos(phi[q]) * Bphi[q] + sin(phi[q]) * Br[q];
                }
            }
        }
    }
    else if (m_session->DefinesFunction("MagneticField"))
    {
        GetFunction("MagneticField")->Evaluate(Bstring, B_in);
    }

    this->B = Array<OneD, MR::DisContFieldSharedPtr>(3);
    for (d = 0; d < 3; ++d)
    {
        B[d] = MemoryManager<MR::DisContField>::AllocateSharedPtr(
            *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]));
        B[d]->UpdatePhys() = B_in[d];
        B[d]->FwdTrans(B[d]->GetPhys(), B[d]->UpdateCoeffs());
    }

    this->mag_B = Array<OneD, NekDouble>(npoints, 0.0);
    for (d = 0; d < 3; ++d)
    {
        Vmath::Vvtvp(npoints, B[d]->GetPhys(), 1, B[d]->GetPhys(), 1, mag_B, 1,
                     mag_B, 1);
    }

    for (d = 0; d < 3; d++)
    {
        b_unit[d] = Array<OneD, NekDouble>(npoints, 0.0);
        for (int k = 0; k < npoints; ++k)
        {
            this->b_unit[d][k] =
                (this->mag_B[k] > 0)
                    ? B[d]->GetPhys()[k] / std::sqrt(this->mag_B[k])
                    : 0.0;
        }
    }
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

    ReadMagneticField();
    nSpecies = 2;

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
    int i;
    int num_vars = in_arr.size();
    int npoints  = GetNpoints();
    switch (m_projectionType)
    {
        case MultiRegions::eDiscontinuous:
        {
            // Just copy over array
            if (in_arr != out_arr)
            {
                int npoints = GetNpoints();

                for (i = 0; i < num_vars; ++i)
                {
                    Vmath::Vcopy(npoints, in_arr[i], 1, out_arr[i], 1);
                }
            }
            break;
        }
        case MultiRegions::eGalerkin:
        {
            Array<OneD, NekDouble> coeffs(m_fields[0]->GetNcoeffs());

            for (i = 0; i < num_vars; ++i)
            {
                m_fields[i]->FwdTrans(in_arr[i], coeffs);
                m_fields[i]->BwdTrans(coeffs, out_arr[i]);
            }
            break;
        }
        default:
        {
            ASSERTL0(false, "Unknown projection scheme");
            break;
        }
    }
    SetBoundaryConditions(out_arr, time);

}

void TokamakSystem::v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables)
{
    const int nPhys   = m_fields[0]->GetNpoints();
    const int nCoeffs = m_fields[0]->GetNcoeffs();

    variables.push_back("Bx");
    Array<OneD, NekDouble> BxFwd(nCoeffs);
    m_fields[0]->FwdTransLocalElmt(B[0]->GetPhys(), BxFwd);
    fieldcoeffs.push_back(BxFwd);

    variables.push_back("By");
    Array<OneD, NekDouble> ByFwd(nCoeffs);
    m_fields[0]->FwdTransLocalElmt(B[1]->GetPhys(), ByFwd);
    fieldcoeffs.push_back(ByFwd);

    variables.push_back("Bz");
    Array<OneD, NekDouble> BzFwd(nCoeffs);
    m_fields[0]->FwdTransLocalElmt(B[2]->GetPhys(), BzFwd);
    fieldcoeffs.push_back(BzFwd);

    variables.push_back("Ex");
    Array<OneD, NekDouble> ExFwd(nCoeffs);
    m_fields[0]->FwdTransLocalElmt(E[0]->GetPhys(), ExFwd);
    fieldcoeffs.push_back(ExFwd);

    variables.push_back("Ey");
    Array<OneD, NekDouble> EyFwd(nCoeffs);
    m_fields[0]->FwdTransLocalElmt(E[1]->GetPhys(), EyFwd);
    fieldcoeffs.push_back(EyFwd);

    variables.push_back("Ez");
    Array<OneD, NekDouble> EzFwd(nCoeffs);
    m_fields[0]->FwdTransLocalElmt(E[2]->GetPhys(), EzFwd);
    fieldcoeffs.push_back(EzFwd);

    for (auto s : this->particle_sys->get_species())
    {
        variables.push_back(s.second.name + "_SOURCE_DENSITY");
        Array<OneD, NekDouble> SrcFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(src_fields[s.first]->GetPhys(), SrcFwd);
        fieldcoeffs.push_back(SrcFwd);
    }

    variables.push_back("Rank");
    Array<OneD, NekDouble> Rank(nPhys, this->m_session->GetComm()->GetRank());
    Array<OneD, NekDouble> RankFwd(nCoeffs);
    m_fields[0]->FwdTransLocalElmt(Rank, RankFwd);
    fieldcoeffs.push_back(RankFwd);
}

void TokamakSystem::v_GenerateSummary(SU::SummaryList &s)
{
    TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>::v_GenerateSummary(s);
    SU::AddSummaryItem(s, "Riemann solver", this->riemann_solver_type);
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

    // Turn off forward-transform of initial conditions.
    m_homoInitialFwd = false;

    // Store FieldSharedPtr casts of fields in a map, indexed by name

    this->E = Array<OneD, MR::DisContFieldSharedPtr>(3);
    for (int d = 0; d < 3; ++d)
    {
        this->E[d] = MemoryManager<MR::DisContField>::AllocateSharedPtr(
            *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]));
    }

    // Bind projection function for time integration object
    m_ode.DefineProjection(&TokamakSystem::DoOdeProjection, this);
    int nConvectiveFields = m_fields.size();

    m_implHelper = std::make_shared<ImplicitHelper>(m_session, m_fields, m_ode,
                                                    nConvectiveFields);
    m_implHelper->InitialiseNonlinSysSolver();
    m_ode.DefineImplicitSolve(&ImplicitHelper::ImplicitTimeInt, m_implHelper);

    // Forcing terms for coupling to Reactions
    m_forcing = SU::Forcing::Load(m_session, shared_from_this(), m_fields,
                                  m_fields.size());

    Array<OneD, const SD::BoundaryConditionShPtr> BndConds =
        m_fields[0]->GetBndConditions();

    for (size_t n = 0; n < BndConds.size(); ++n)
    {
        std::string type = BndConds[n]->GetUserDefined();
        if (type.rfind("ObliqueOutflow", 0) == 0)
        {
            ASSERTL0(BndConds[n]->GetBoundaryConditionType() == SD::eNeumann,
                     "Oblique outflow boundary condition must be of type "
                     "Neumann <N>");
        }

        else if (type.rfind("Oblique", 0) == 0)
        {
            ASSERTL0(BndConds[n]->GetBoundaryConditionType() == SD::eNeumann,
                     "Oblique boundary condition must be of type Neumann <N>");
        }
        else if (type.rfind("Sheath", 0) == 0)
        {
            ASSERTL0(BndConds[n]->GetBoundaryConditionType() == SD::eDirichlet,
                     "Sheath boundary condition must be of type Dirichlet <D>");
        }
        if (BndConds[n]->GetBoundaryConditionType() == SD::ePeriodic)
        {
            continue;
        }

        if (!type.empty())
        {
            Array<OneD, Array<OneD, NekDouble>> magneticFieldBndElmt(3);
            for (int d = 0; d < 3; d++)
            {
                m_fields[0]->ExtractPhysToBnd(n, b_unit[d],
                                              magneticFieldBndElmt[d]);
            }
            m_bndConds.push_back(GetTokamakBndCondFactory().CreateInstance(
                type, m_session, m_fields, magneticFieldBndElmt, m_spacedim,
                n));
        }
    }

    SetBoundaryConditionsBwdWeight();

    if (this->particles_enabled)
    {
        // Set up object to evaluate density field
        this->particle_sys->setup_evaluate_E(this->E[0], this->E[1],
                                             this->E[2]);
        this->particle_sys->setup_evaluate_B(this->B[0], this->B[1],
                                             this->B[2]);
        this->particle_sys->setup_evaluate_ne(
            std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]));

        for (int i = 0; i < this->particle_sys->get_species().size(); ++i)
        {
            this->src_fields.emplace_back(
                MemoryManager<MR::DisContField>::AllocateSharedPtr(
                    *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0])));
        }

        this->particle_sys->setup_project(this->src_fields);
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
        this->particle_sys->write(step + 1);
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
        if (step == 0)
        {
            this->particle_sys->write(0);
        }

        this->particle_sys->integrate(m_time + m_timestep, this->part_timestep);
        this->particle_sys->project_source_terms();
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
