#ifndef TOKAMAKSYSTEM_HPP
#define TOKAMAKSYSTEM_HPP

#include "../BoundaryConditions/TokamakBndCond.hpp"
#include "../Diagnostics/GrowthRatesRecorder.hpp"
#include "../ParticleSystems/ParticleSystem.hpp"

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

    /// Object to facilitate allows optional recording of energy and enstrophy
    std::shared_ptr<GrowthRatesRecorder<MR::DisContField>>
        energy_enstrophy_recorder;
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
    // B in cylindrical polar (r, z, theta)
    Array<OneD, Array<OneD, NekDouble>> B_pol;
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
    virtual void v_GenerateSummary(SU::SummaryList &s) override;
    virtual void v_InitObject(bool DeclareField) override;
    virtual bool v_PostIntegrate(int step) override;
    virtual bool v_PreIntegrate(int step) override;
    virtual void v_SetInitialConditions(NekDouble init_time, bool dump_ICs,
                                        const int domain) override;
};

} // namespace NESO::Solvers::tokamak
#endif
