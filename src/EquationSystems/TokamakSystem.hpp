#ifndef TOKAMAKSYSTEM_HPP
#define TOKAMAKSYSTEM_HPP

#include "../BoundaryConditions/TokamakBndCond.hpp"
#include "../ParticleSystems/ParticleSystem.hpp"

#include "../../NESO/include/nektar_interface/solver_base/neso_session_function.hpp"
#include "nektar_interface/solver_base/time_evolved_eqnsys_base.hpp"
#include "nektar_interface/utilities.hpp"

#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Core/Misc.h>
#include <SolverUtils/Diffusion/Diffusion.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Forcing/Forcing.h>

#include "ImplicitHelper.hpp"
#include <solvers/solver_callback_handler.hpp>

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;
namespace SU = Nektar::SolverUtils;

namespace NESO::Solvers::tokamak
{

/**
 * @brief Equation system for the Tokamak proxyapp
 *
 */
class TokamakSystem
    : public TimeEvoEqnSysBase<SU::UnsteadySystem, ParticleSystem>
{
public:
    friend class MemoryManager<TokamakSystem>;

    /// Name of class
    static std::string class_name;

    /**
     * @brief Create an instance of this class and initialise it.
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

    virtual std::shared_ptr<ParticleSystem> GetParticleSystem();

    /// Callback handler to call user-defined callbacks
    SolverCallbackHandler<TokamakSystem> solver_callback_handler;

protected:
    TokamakSystem(const LU::SessionReaderSharedPtr &session,
                  const SD::MeshGraphSharedPtr &graph);

    bool transient_field;
    /// Magnetic field vector
    Array<OneD, MR::DisContFieldSharedPtr> B;
    /// Normalised magnetic field vector
    Array<OneD, Array<OneD, NekDouble>> b_unit;

    /// Magnitude of the magnetic field
    Array<OneD, NekDouble> mag_B;
    /// Trace of magnetic field
    Array<OneD, Array<OneD, NekDouble>> m_magneticFieldTrace;

    /// Electric Field
    Array<OneD, MR::DisContFieldSharedPtr> E;

    /// Advection object used in the electron density equation
    SU::AdvectionSharedPtr m_advection;
    /// Advection type
    std::string adv_type;

    /// Diffusion object used in anisotropic diffusion
    SU::DiffusionSharedPtr m_diffusion;

    /// Sheath potential
    NekDouble lambda;

    MR::DisContFieldSharedPtr ne;
    MR::DisContFieldSharedPtr Te;
    Array<OneD, MR::DisContFieldSharedPtr> ve;

    Array<OneD, MR::ExpListSharedPtr> m_allfields;
    Array<OneD, MR::ExpListSharedPtr> m_saved;
    Array<OneD, MR::ExpListSharedPtr> m_indfields;
    int n_indep_fields;
    int n_species;
    int n_fields_per_species;

    /** Density source fields cast to DisContFieldSharedPtr for use in
     * particle evaluation/projection methods
     */
    std::vector<MR::DisContFieldSharedPtr> src_fields;
    std::vector<int> components;
    std::vector<Sym<REAL>> src_syms;

    /// Bool to enable/disable growth rate recordings
    bool energy_enstrophy_recording_enabled;

    /// Boundary Conditions
    std::vector<TokamakBndCondSharedPtr> m_bndConds;

    /// Forcing terms
    std::vector<SolverUtils::ForcingSharedPtr> m_forcing;

    // Convenience key for varcoeffmaps
    static constexpr StdRegions::VarCoeffType vc[3][3] = {
        {StdRegions::eVarCoeffD00, StdRegions::eVarCoeffD01,
         StdRegions::eVarCoeffD02},
        {StdRegions::eVarCoeffD01, StdRegions::eVarCoeffD11,
         StdRegions::eVarCoeffD12},
        {StdRegions::eVarCoeffD02, StdRegions::eVarCoeffD12,
         StdRegions::eVarCoeffD22}};

    /// Number of particle timesteps per fluid timestep.
    int num_part_substeps;
    /// Number of time steps between particle trajectory step writes.
    int particle_output_freq;
    /// Particle timestep size.
    double part_timestep;

    /// Riemann solver type (used for all advection terms)
    std::string riemann_solver_type;
    /// Riemann solver object used in electron advection
    SU::RiemannSolverSharedPtr riemann_solver;

    std::shared_ptr<ImplicitHelper> m_implHelper;

    virtual void load_params() override;
    void ReadMagneticField(NekDouble time = 0);

    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                  Array<OneD, Array<OneD, NekDouble>> &out_arr,
                  const NekDouble time);

    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
        Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time);

    void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble>> &physarray,
                               NekDouble time);
    void SetBoundaryConditionsBwdWeight();
    virtual void v_ExtraFldOutput(
        std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
        std::vector<std::string> &variables) override;
    void v_DoSolve() override;
    virtual void v_GenerateSummary(SU::SummaryList &s) override;
    virtual void v_InitObject(bool DeclareField) override;
    virtual bool v_PostIntegrate(int step) override;
    virtual bool v_PreIntegrate(int step) override;
    virtual void v_SetInitialConditions(NekDouble init_time, bool dump_ICs,
                                        const int domain) override;

    NESOSessionFunctionSharedPtr get_species_function(
        int s, std::string name,
        const MR::ExpListSharedPtr &field = MR::NullExpListSharedPtr,
        bool cache                        = false);
    std::vector<std::map<std::string, NESOSessionFunctionSharedPtr>>
        m_nesoSessionFunctions;
};

} // namespace NESO::Solvers::tokamak
#endif
