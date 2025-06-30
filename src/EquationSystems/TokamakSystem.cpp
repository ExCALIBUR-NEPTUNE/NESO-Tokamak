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
}

std::shared_ptr<ParticleSystem> TokamakSystem::GetParticleSystem()
{
    return this->particle_sys;
}

/**
 * @brief Read the magnetic field from file.
 */
void TokamakSystem::ReadMagneticField(NekDouble time)
{
    int d;
    int npoints = m_fields[0]->GetNpoints();
    Array<OneD, Array<OneD, NekDouble>> B_in(3);
    for (d = 0; d < 3; ++d)
    {
        B_in[d] = Array<OneD, NekDouble>(npoints);
    }

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

                Bxfunc->Evaluate(r, x1, phi, time, Br);
                Byfunc->Evaluate(r, x1, phi, time, B_in[1]);
                Bzfunc->Evaluate(r, x1, phi, time, Bphi);

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
        GetFunction("MagneticField")->Evaluate(Bstring, B_in, time);
    }

    this->B = Array<OneD, MR::DisContFieldSharedPtr>(3);
    for (d = 0; d < 3; ++d)
    {
        this->B[d] = MemoryManager<MR::DisContField>::AllocateSharedPtr(
            *std::dynamic_pointer_cast<MR::DisContField>(m_fields[0]));
        this->B[d]->UpdatePhys() = B_in[d];
        this->B[d]->FwdTrans(B[d]->GetPhys(), B[d]->UpdateCoeffs());
    }

    this->mag_B = Array<OneD, NekDouble>(npoints, 0.0);
    for (d = 0; d < 3; ++d)
    {
        Vmath::Vvtvp(npoints, B[d]->GetPhys(), 1, B[d]->GetPhys(), 1, mag_B, 1,
                     mag_B, 1);
    }

    for (d = 0; d < 3; d++)
    {
        this->b_unit[d] = Array<OneD, NekDouble>(npoints, 0.0);
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
    std::string transient_field_str;
    m_session->LoadSolverInfo("MagneticFieldEvolution", transient_field_str,
                              "Static");
    this->transient_field = (transient_field_str == "Transient");

    // Type of Riemann solver to use. Default = "Upwind"
    m_session->LoadSolverInfo("UpwindType", this->riemann_solver_type,
                              "MultiFieldUpwind");

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
                m_indfields[i]->FwdTrans(in_arr[i], coeffs);
                m_indfields[i]->BwdTrans(coeffs, out_arr[i]);
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
    bool extraFields;
    m_session->MatchSolverInfo("OutputSpeciesFields", "True", extraFields,
                               true);
    if (extraFields)
    {
        const int nPhys   = m_fields[0]->GetNpoints();
        const int nCoeffs = m_fields[0]->GetNcoeffs();
        for (const auto &[k, v] : this->neso_config->get_species())
        {
            std::string name = std::get<0>(v);

            for (int f = 0; f < this->n_fields_per_species; ++f)
            {
                int fi = f + k * this->n_fields_per_species;
                variables.push_back(m_session->GetVariable(f) + "_" + name);
                Array<OneD, NekDouble> Fwd(nCoeffs);
                this->m_indfields[fi]->FwdTransLocalElmt(
                    this->m_indfields[fi]->GetPhys(), Fwd);
                fieldcoeffs.push_back(Fwd);
            }
        }
        if (Te)
        {
            variables.push_back("Te");
            Array<OneD, NekDouble> Fwd(nCoeffs);
            this->Te->FwdTransLocalElmt(this->Te->GetPhys(), Fwd);
            fieldcoeffs.push_back(Fwd);
        }
    }

    m_session->MatchSolverInfo("OutputEMFields", "True", extraFields, true);
    if (extraFields)
    {
        const int nPhys   = m_fields[0]->GetNpoints();
        const int nCoeffs = m_fields[0]->GetNcoeffs();

        variables.push_back("Bx");
        Array<OneD, NekDouble> BxFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->B[0]->GetPhys(), BxFwd);
        fieldcoeffs.push_back(BxFwd);

        variables.push_back("By");
        Array<OneD, NekDouble> ByFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->B[1]->GetPhys(), ByFwd);
        fieldcoeffs.push_back(ByFwd);

        variables.push_back("Bz");
        Array<OneD, NekDouble> BzFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->B[2]->GetPhys(), BzFwd);
        fieldcoeffs.push_back(BzFwd);

        variables.push_back("Ex");
        Array<OneD, NekDouble> ExFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->E[0]->GetPhys(), ExFwd);
        fieldcoeffs.push_back(ExFwd);

        variables.push_back("Ey");
        Array<OneD, NekDouble> EyFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->E[1]->GetPhys(), EyFwd);
        fieldcoeffs.push_back(EyFwd);

        variables.push_back("Ez");
        Array<OneD, NekDouble> EzFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(this->E[2]->GetPhys(), EzFwd);
        fieldcoeffs.push_back(EzFwd);
    }
    m_session->MatchSolverInfo("OutputPartitions", "True", extraFields, false);
    if (extraFields)
    {
        const int nPhys   = m_fields[0]->GetNpoints();
        const int nCoeffs = m_fields[0]->GetNcoeffs();
        variables.push_back("Rank");
        Array<OneD, NekDouble> Rank(nPhys,
                                    this->m_session->GetComm()->GetRank());
        Array<OneD, NekDouble> RankFwd(nCoeffs);
        m_fields[0]->FwdTransLocalElmt(Rank, RankFwd);
        fieldcoeffs.push_back(RankFwd);
    }
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

    this->ve = Array<OneD, MR::DisContFieldSharedPtr>(3);
    ReadMagneticField(0);
    this->n_species = this->neso_config->get_species().size();
    m_allfields     = Array<OneD, MR::ExpListSharedPtr>(
        m_fields.size() + this->n_species * n_fields_per_species);
    m_indfields = Array<OneD, MR::ExpListSharedPtr>(
        n_indep_fields + this->n_species * n_fields_per_species);

    for (int i = 0; i < m_fields.size(); ++i)
    {
        m_allfields[i] = m_fields[i];
    }

    for (const auto &[k, v] : this->neso_config->get_species())
    {
        std::string name = std::get<0>(v);
        for (int f = 0; f < this->n_fields_per_species; ++f)
        {
            if (m_projectionType == MR::eGalerkin ||
                m_projectionType == MR::eMixed_CG_Discontinuous)
            {
                m_indfields[k * n_fields_per_species + f] =
                    MemoryManager<MR::ContField>::AllocateSharedPtr(
                        *std::dynamic_pointer_cast<MR::ContField>(m_fields[f]),
                        m_graph, m_session->GetVariable(f), true,
                        m_checkIfSystemSingular[f]);
            }
            else
            {
                m_indfields[k * n_fields_per_species + f] =
                    MemoryManager<MR::DisContField>::AllocateSharedPtr(
                        *std::dynamic_pointer_cast<MR::DisContField>(
                            m_fields[f]),
                        m_graph, m_session->GetVariable(f));
            }
        }
    }
    for (int i = 0; i < n_indep_fields; ++i)
    {
        m_indfields[n_fields_per_species * n_species + i] =
            m_fields[m_fields.size() - n_indep_fields + i];
    }

    // Bind projection function for time integration object
    m_ode.DefineProjection(&TokamakSystem::DoOdeProjection, this);
    if (m_projectionType == MR::eDiscontinuous)
    {
        m_implHelper = std::make_shared<ImplicitHelper>(
            m_session, m_indfields, m_ode, m_indfields.size());
        m_implHelper->InitialiseNonlinSysSolver();
        m_ode.DefineImplicitSolve(&ImplicitHelper::ImplicitTimeInt,
                                  m_implHelper);
    }

    // Forcing terms
    m_forcing = SU::Forcing::Load(m_session, shared_from_this(), m_indfields,
                                  m_indfields.size());

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
}

/**
 * @brief Initialises the time integration scheme (as specified in the
 * session file), and perform the time integration.
 */
void TokamakSystem::v_DoSolve()
{
    ASSERTL0(m_intScheme != nullptr, "No time integration scheme.");

    int i          = 1;
    int nvariables = 0;
    int nfields    = m_indfields.size();
    if (m_intVariables.empty())
    {
        for (i = 0; i < nfields; ++i)
        {
            m_intVariables.push_back(i);
        }
        nvariables = nfields;
    }
    else
    {
        nvariables = m_intVariables.size();
    }

    // Set up wrapper to fields data storage.
    Array<OneD, Array<OneD, NekDouble>> fields(nvariables);

    // Order storage to list time-integrated fields first.
    // @TODO: Refactor to take coeffs (FwdTrans) if boolean flag (in constructor
    // function) says to.
    for (i = 0; i < nvariables; ++i)
    {
        fields[i] = m_indfields[m_intVariables[i]]->UpdatePhys();
        m_indfields[m_intVariables[i]]->SetPhysState(false);
    }

    // @TODO: Virtual function that allows to transform the field space, embed
    // the MultiplyMassMatrix in here.
    // @TODO: Specify what the fields variables are physical or coefficient,
    // boolean in UnsteadySystem class...

    v_ALEPreMultiplyMass(fields);

    // Initialise time integration scheme.
    m_intScheme->InitializeScheme(m_timestep, fields, m_time, m_ode);

    LibUtilities::Timer timer;
    bool doCheckTime      = false;
    int step              = m_initialStep;
    int stepCounter       = 0;
    NekDouble intTime     = 0.0;
    NekDouble cpuTime     = 0.0;
    NekDouble cpuPrevious = 0.0;
    NekDouble elapsed     = 0.0;

    m_lastCheckTime = 0.0;

    Array<OneD, int> abortFlags(2, 0);
    std::string abortFile = "abort";
    if (m_session->DefinesSolverInfo("CheckAbortFile"))
    {
        abortFile = m_session->GetSolverInfo("CheckAbortFile");
    }

    NekDouble tmp_cflSafetyFactor = m_cflSafetyFactor;
    while ((step < m_steps || m_time < m_fintime - NekConstants::kNekZeroTol) &&
           abortFlags[1] == 0)
    {
        if (m_CFLGrowth > 1.0 && m_cflSafetyFactor < m_CFLEnd)
        {
            tmp_cflSafetyFactor =
                std::min(m_CFLEnd, m_CFLGrowth * tmp_cflSafetyFactor);
        }

        // Frozen preconditioner checks.
        if (!m_comm->IsParallelInTime())
        {
            if (v_UpdateTimeStepCheck())
            {
                m_cflSafetyFactor = tmp_cflSafetyFactor;

                if (m_cflSafetyFactor)
                {
                    m_timestep = GetTimeStep(fields);
                }

                // Ensure that the final timestep finishes at the final
                // time, or at a prescribed IO_CheckTime.
                if (m_time + m_timestep > m_fintime && m_fintime > 0.0)
                {
                    m_timestep = m_fintime - m_time;
                }
                else if (m_checktime &&
                         m_time + m_timestep - m_lastCheckTime >= m_checktime)
                {
                    m_lastCheckTime += m_checktime;
                    m_timestep  = m_lastCheckTime - m_time;
                    doCheckTime = true;
                }
            }
        }

        if (m_TimeIncrementFactor > 1.0)
        {
            NekDouble timeincrementFactor = m_TimeIncrementFactor;
            m_timestep *= timeincrementFactor;

            if (m_time + m_timestep > m_fintime && m_fintime > 0.0)
            {
                m_timestep = m_fintime - m_time;
            }
        }

        // Perform any solver-specific pre-integration steps.
        timer.Start();
        if (v_PreIntegrate(
                step)) // Could be possible to put a preintegrate step in the
                       // ALEHelper class, put in the Unsteady Advection class
        {
            break;
        }

        ASSERTL0(m_timestep > 0, "m_timestep < 0");

        fields = m_intScheme->TimeIntegrate(stepCounter, m_timestep);
        timer.Stop();

        m_time += m_timestep;
        elapsed = timer.TimePerTest(1);
        intTime += elapsed;
        cpuTime += elapsed;

        // Write out status information.
        v_PrintStatusInformation(step, cpuTime);
        if (m_infosteps &&
            m_session->GetComm()->GetSpaceComm()->GetRank() == 0 &&
            !((step + 1) % m_infosteps))
        {
            cpuPrevious = cpuTime;
            cpuTime     = 0.0;
        }

        // Transform data into coefficient space
        for (i = 0; i < nvariables; ++i)
        {
            // copy fields into ExpList::m_phys and assign the new
            // array to fields
            m_indfields[m_intVariables[i]]->SetPhys(fields[i]);
            fields[i] = m_indfields[m_intVariables[i]]->UpdatePhys();
            if (v_RequireFwdTrans())
            {
                if (m_comm->IsParallelInTime())
                {
                    m_indfields[m_intVariables[i]]->FwdTrans(
                        m_indfields[m_intVariables[i]]->GetPhys(),
                        m_indfields[m_intVariables[i]]->UpdateCoeffs());
                }
                else
                {
                    m_indfields[m_intVariables[i]]->FwdTransLocalElmt(
                        m_indfields[m_intVariables[i]]->GetPhys(),
                        m_indfields[m_intVariables[i]]->UpdateCoeffs());
                }
            }
            m_indfields[m_intVariables[i]]->SetPhysState(false);
        }

        // Perform any solver-specific post-integration steps.
        if (v_PostIntegrate(step))
        {
            break;
        }

        // Test for abort conditions (nan, or abort file).
        if (m_abortSteps && !((step + 1) % m_abortSteps))
        {
            abortFlags[0] = 0;
            for (i = 0; i < nvariables; ++i)
            {
                if (Vmath::Nnan(
                        m_indfields[m_intVariables[i]]->GetPhys().size(),
                        m_indfields[m_intVariables[i]]->GetPhys(), 1) > 0)
                {
                    abortFlags[0] = 1;
                }
            }

            // Rank zero looks for abort file and deletes it
            // if it exists. Communicates the abort.
            if (m_session->GetComm()->GetSpaceComm()->GetRank() == 0)
            {
                if (fs::exists(abortFile))
                {
                    fs::remove(abortFile);
                    abortFlags[1] = 1;
                }
            }

            m_session->GetComm()->GetSpaceComm()->AllReduce(
                abortFlags, LibUtilities::ReduceMax);

            ASSERTL0(!abortFlags[0], "NaN found during time integration.");
        }

        // Write out checkpoint files.
        if ((m_checksteps && !((step + 1) % m_checksteps)) || doCheckTime)
        {

            Checkpoint_Output(m_nchk);
            m_nchk++;

            doCheckTime = false;
        }

        // Step advance.
        ++step;
        ++stepCounter;
    }

    // Print out summary statistics.
    v_PrintSummaryStatistics(intTime);

    // If homogeneous, transform back into physical space if necessary.

    for (i = 0; i < nvariables; ++i)
    {
        m_indfields[m_intVariables[i]]->SetPhysState(true);
    }
    if (this->particle_sys)
    {
        this->particle_sys->free();
    }
}

/**
 * @brief Override v_PostIntegrate to do particle output and compute
 * diagnostics, if enabled
 * @param step Time step number
 */
bool TokamakSystem::v_PostIntegrate(int step)
{
    this->solver_callback_handler.call_post_integrate(this);

    // Writes a step of the particle trajectory.

    return TimeEvoEqnSysBase<SU::UnsteadySystem,
                             ParticleSystem>::v_PostIntegrate(step);
}

/**
 * @brief Override v_PreIntegrate to do particle system integration
 *
 * @param step Time step number
 */
bool TokamakSystem::v_PreIntegrate(int step)
{
    this->solver_callback_handler.call_pre_integrate(this);

    if (this->transient_field)
    {
        ReadMagneticField(m_time);
    }

    // for (int f = 0; f < this->n_fields_per_species; ++f)
    //{
    Vmath::Zero(n_pts, m_fields[0]->UpdatePhys(), 1);
    for (const auto &[k, v] : this->neso_config->get_species())
    {
        Vmath::Vadd(n_pts, m_fields[0]->GetPhys(), 1,
                    m_indfields[k * n_fields_per_species]->GetPhys(), 1,
                    m_fields[0]->UpdatePhys(), 1);
    }
    m_fields[0]->FwdTransLocalElmt(this->m_fields[0]->GetPhys(),
                                   m_fields[0]->UpdateCoeffs());

    if (this->Te)
    {
        Vmath::Vdiv(n_pts, m_indfields[m_indfields.size() - 1]->GetPhys(), 1,
                    m_fields[0]->GetPhys(), 1, Te->UpdatePhys(), 1);
        Te->FwdTransLocalElmt(this->Te->GetPhys(), Te->UpdateCoeffs());
    }
    //}

    if (this->particles_enabled)
    {
        // Integrate the particle system to the requested time.
        this->particle_sys->evaluate_fields(this->E, this->B, this->ne,
                                            this->Te, this->ve);
        if (particle_output_freq > 0 && (step % particle_output_freq) == 0)
        {
            this->particle_sys->write(step);
        }
        for(auto& fld : this->src_fields)
        {
            Vmath::Zero(fld->GetNpoints(), fld->UpdatePhys(), 1);
        }
        this->particle_sys->integrate(m_time + m_timestep, this->part_timestep);
    }

    return UnsteadySystem::v_PreIntegrate(step);
}

NESOSessionFunctionSharedPtr TokamakSystem::get_species_function(
    int s, std::string name, const MR::ExpListSharedPtr &field, bool cache)
{
    MR::ExpListSharedPtr vField = field;
    if (!field)
    {
        vField = m_fields[0];
    }

    if (cache)
    {
        if ((m_nesoSessionFunctions[s].find(name) ==
             m_nesoSessionFunctions[s].end()) ||
            (m_nesoSessionFunctions[s][name]->GetSession() !=
             this->neso_config) ||
            (m_nesoSessionFunctions[s][name]->GetExpansion() != vField))
        {
            m_nesoSessionFunctions[s][name] =
                MemoryManager<NESOSessionFunction>::AllocateSharedPtr(
                    s, this->neso_config, vField, name, cache);
        }

        return m_nesoSessionFunctions[s][name];
    }
    else
    {
        return std::make_shared<NESOSessionFunction>(s, this->neso_config,
                                                     vField, name, cache);
    }
}

/**
 * @brief After reading ICs, calculate phi and grad(phi)
 */
void TokamakSystem::v_SetInitialConditions(NekDouble init_time, bool dump_ICs,
                                           const int domain)
{
    if (m_session->DefinesFunction("InitialConditions"))
    {
        GetFunction("InitialConditions")
            ->Evaluate(m_session->GetVariables(), m_fields, m_time, domain);
        // Enforce C0 Continutiy of initial condiiton
        if ((m_projectionType == MultiRegions::eGalerkin) ||
            (m_projectionType == MultiRegions::eMixed_CG_Discontinuous))
        {
            for (int i = 0; i < m_fields.size(); ++i)
            {
                m_fields[i]->LocalToGlobal();
                m_fields[i]->GlobalToLocal();
                m_fields[i]->BwdTrans(m_fields[i]->GetCoeffs(),
                                      m_fields[i]->UpdatePhys());
            }
        }

        if (m_session->GetComm()->GetRank() == 0)
        {
            for (int i = 0; i < m_fields.size(); ++i)
            {
                std::string varName = m_session->GetVariable(i);
                std::cout << "  - Field " << varName << ": "
                          << GetFunction("InitialConditions")
                                 ->Describe(varName, domain)
                          << std::endl;
            }
        }
    }
    else
    {
        int nq = m_fields[0]->GetNpoints();
        for (int i = 0; i < m_fields.size(); i++)
        {
            Vmath::Zero(nq, m_fields[i]->UpdatePhys(), 1);
            m_fields[i]->SetPhysState(true);
            Vmath::Zero(m_fields[i]->GetNcoeffs(), m_fields[i]->UpdateCoeffs(),
                        1);
            if (m_session->GetComm()->GetRank() == 0)
            {
                std::cout << "  - Field " << m_session->GetVariable(i)
                          << ": 0 (default)" << std::endl;
            }
        }
    }
    for (int s = 0; s < this->n_species; ++s)
    {
        if (this->neso_config->defines_species_function(s, "InitialConditions"))
        {
            auto fn = this->get_species_function(s, "InitialConditions");
            for (int f = 0; f < n_fields_per_species; ++f)
            {
                int fi = f + s * n_fields_per_species;
                fn->Evaluate(m_session->GetVariables()[f],
                             m_indfields[fi]->UpdatePhys(), m_time, domain);
                if (m_indfields[fi]->GetWaveSpace())
                {
                    m_indfields[fi]->HomogeneousFwdTrans(
                        m_indfields[fi]->GetTotPoints(),
                        m_indfields[fi]->GetPhys(),
                        m_indfields[fi]->UpdatePhys());
                }
                m_indfields[fi]->FwdTransLocalElmt(
                    m_indfields[fi]->GetPhys(),
                    m_indfields[fi]->UpdateCoeffs());

                if ((m_projectionType == MultiRegions::eGalerkin) ||
                    (m_projectionType == MultiRegions::eMixed_CG_Discontinuous))
                {

                    m_indfields[fi]->LocalToGlobal();
                    m_indfields[fi]->GlobalToLocal();
                    m_indfields[fi]->BwdTrans(m_indfields[fi]->GetCoeffs(),
                                              m_indfields[fi]->UpdatePhys());
                }
            }

            // Enforce C0 Continutiy of initial condiiton
        }
        else
        {
            int nq = m_indfields[0]->GetNpoints();
            for (int f = 0; f < this->n_fields_per_species; f++)
            {
                int fi = f + s * n_fields_per_species;
                Vmath::Zero(nq, m_indfields[fi]->UpdatePhys(), 1);
                m_indfields[fi]->SetPhysState(true);
                Vmath::Zero(m_indfields[fi]->GetNcoeffs(),
                            m_indfields[fi]->UpdateCoeffs(), 1);
            }
        }
    }

    if (dump_ICs && m_checksteps && m_nchk == 0 && !m_comm->IsParallelInTime())
    {
        Checkpoint_Output(m_nchk);
    }
    ++m_nchk;

    if (this->particle_sys)
    {
        this->particle_sys->initialise_particles_from_fields(
            this->E, this->B, this->ne, this->Te, this->ve);
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
