#ifndef TOKAMAKSYSTEM_H
#define TOKAMAKSYSTEM_H

#include "../BoundaryConditions/TokamakBndCond.h"
#include "../Diagnostics/GrowthRatesRecorder.hpp"
#include "../ParticleSystems/ParticleSystem.hpp"

#include "nektar_interface/solver_base/time_evolved_eqnsys_base.hpp"
#include "nektar_interface/utilities.hpp"

#include "ImplicitHelper.hpp"
#include <LibUtilities/Memory/NekMemoryManager.hpp>
#include <SolverUtils/AdvectionSystem.h>
#include <SolverUtils/Core/Misc.h>
#include <SolverUtils/Diffusion/Diffusion.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

#include <solvers/solver_callback_handler.hpp>

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

    /// Object to facilitate allows optional recording of energy and enstrophy
    std::shared_ptr<GrowthRatesRecorder<MR::DisContField>>
        energy_enstrophy_recorder;
    /// Callback handler to call user-defined callbacks
    SolverCallbackHandler<TokamakSystem> solver_callback_handler;

protected:
    TokamakSystem(const LU::SessionReaderSharedPtr &session,
                  const SD::MeshGraphSharedPtr &graph);

    /// Simulate in mean fluid or turbulent mode
    std::string m_mode;

    /// Advection object used in the electron density equation
    SU::AdvectionSharedPtr m_advection;
    /// Advection type
    std::string adv_type;

    /// Diffusion object used in anisotropic diffusion
    SU::DiffusionSharedPtr m_diffusion;

    /// Hasegawa-Wakatani α
    NekDouble alpha;
    /// Hasegawa-Wakatani κ
    NekDouble kappa;

    /// Sheath potential
    NekDouble lambda;

    /// Magnetic field vector
    Array<OneD, Array<OneD, NekDouble>> B;
    /// Normalised magnetic field vector
    Array<OneD, Array<OneD, NekDouble>> b_unit;
    // B in cylindrical polar (r, z, theta)
    Array<OneD, Array<OneD, NekDouble>> B_pol;

    /// Magnitude of the magnetic field
    Array<OneD, NekDouble> mag_B;

    /// Trace of magnetic field
    Array<OneD, Array<OneD, NekDouble>> m_magneticFieldTrace;

    /** Source fields cast to DisContFieldSharedPtr, indexed by name, for use in
     * particle evaluation/projection methods
     */
    std::map<std::string, MR::DisContFieldSharedPtr> discont_fields;
    /// Bool to enable/disable growth rate recordings
    bool energy_enstrophy_recording_enabled;
    /// Storage for ExB drift velocity
    Array<OneD, Array<OneD, NekDouble>> ExB_vel;



    /// Boundary Conditions
    std::vector<TokamakBndCondSharedPtr> m_bndConds;

    void CalcInitPhi();
    void SolvePhi(const Array<OneD, const Array<OneD, NekDouble>> &in_arr);
    void ComputeGradPhi();

    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                  Array<OneD, Array<OneD, NekDouble>> &out_arr,
                  const NekDouble time);
    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
        Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time);

    void do_null_precon(const Array<OneD, const NekDouble> &in_arr,
                        Array<OneD, NekDouble> &out_arr, const bool &flag);

    Array<OneD, NekDouble> &GetAdvVelNorm(
        Array<OneD, NekDouble> &trace_vel_norm,
        const Array<OneD, Array<OneD, NekDouble>> &adv_vel);

    void GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble>> &fields_vals,
        const Array<OneD, Array<OneD, NekDouble>> &adv_vel,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux);

    void SetBoundaryConditions(Array<OneD, Array<OneD, NekDouble>> &physarray,
                               NekDouble time);
    void SetBoundaryConditionsBwdWeight();

    virtual void v_GenerateSummary(SU::SummaryList &s) override;
    virtual void v_InitObject(bool DeclareField) override;
    virtual bool v_PostIntegrate(int step) override;
    virtual bool v_PreIntegrate(int step) override;
    virtual void v_SetInitialConditions(NekDouble init_time, bool dump_ICs,
                                        const int domain) override;

private:
    StdRegions::VarCoeffMap m_phi_varcoeff;

    // For Diffusion
    NekDouble m_kperp;
    NekDouble m_kpar;
    StdRegions::ConstFactorMap m_factors;
    StdRegions::VarCoeffMap m_varcoeff;

    /// Storage for component of ne advection velocity normal to trace elements
    Array<OneD, NekDouble> norm_vel_elec;


    /// Number of particle timesteps per fluid timestep.
    int num_part_substeps;
    /// Number of time steps between particle trajectory step writes.
    int particle_output_freq;
    /// Particle timestep size.
    double part_timestep;
    /// Riemann solver object used in electron advection
    SU::RiemannSolverSharedPtr riemann_solver;

    std::shared_ptr<ImplicitHelper> m_implHelper;

    Array<OneD, NekDouble> &GetAdvVelNormElec();

    void GetFluxVectorDiff(
        const Array<OneD, Array<OneD, NekDouble>> &in_arr,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &viscous_tensor);
    void GetFluxVectorElec(
        const Array<OneD, Array<OneD, NekDouble>> &fields_vals,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux);
};

} // namespace NESO::Solvers::tokamak
#endif
