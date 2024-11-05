#ifndef TOKAMAKSYSTEM_HPP
#define TOKAMAKSYSTEM_HPP

#include "../BoundaryConditions/TokamakBndCond.hpp"
#include "../Diagnostics/GrowthRatesRecorder.hpp"
#include "../Forcing/ReactionsCoupling.hpp"
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

    int nSpecies;
    /// Diffusion object used in anisotropic diffusion
    SU::DiffusionSharedPtr m_diffusion;

    /// Hasegawa-Wakatani α
    NekDouble alpha;
    /// Hasegawa-Wakatani κ
    NekDouble kappa;

    /// Sheath potential
    NekDouble lambda;

    /// Potential Gradient/-E
    Array<OneD, MR::DisContFieldSharedPtr> m_grad_phi;

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

    virtual void load_params() override;
    void ReadMagneticField();
    void CalcKPar();
    void CalcKPerp();
    void CalcDiffTensor();
    void CalcKappaPar();
    void CalcKappaPerp();
    void CalcKappaTensor();
    void CalcInitPhi();
    void SolvePhi(const Array<OneD, const Array<OneD, NekDouble>> &in_arr);
    void ComputeGradPhi();

    void DoOdeRhsMF(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                    Array<OneD, Array<OneD, NekDouble>> &out_arr,
                    const NekDouble time);
    void DoOdeRhsET(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                    Array<OneD, Array<OneD, NekDouble>> &out_arr,
                    const NekDouble time);

    void DoOdeProjection(
        const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
        Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time);

    Array<OneD, NekDouble> &GetAdvVelNorm(
        Array<OneD, NekDouble> &trace_vel_norm,
        const Array<OneD, Array<OneD, NekDouble>> &adv_vel);

    void GetFluxVector(const Array<OneD, Array<OneD, NekDouble>> &fields_vals,
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
    // Convenience key for varcoeffmaps
    static constexpr StdRegions::VarCoeffType vc[3][3] = {
        {StdRegions::eVarCoeffD00, StdRegions::eVarCoeffD01,
         StdRegions::eVarCoeffD02},
        {StdRegions::eVarCoeffD01, StdRegions::eVarCoeffD11,
         StdRegions::eVarCoeffD12},
        {StdRegions::eVarCoeffD02, StdRegions::eVarCoeffD12,
         StdRegions::eVarCoeffD22}};

    StdRegions::VarCoeffMap m_phi_varcoeff;

    // For Diffusion
    StdRegions::ConstFactorMap m_factors;

    Array<OneD, NekDouble> m_kperp;
    Array<OneD, NekDouble> m_kpar;
    StdRegions::VarCoeffMap m_D;

    Array<OneD, NekDouble> m_kappaperp;
    Array<OneD, NekDouble> m_kappapar;
    StdRegions::VarCoeffMap m_kappa;

    /// Storage for component of ne advection velocity normal to trace elements
    Array<OneD, NekDouble> norm_vel_elec;

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

    /// Forcing terms
    std::vector<SolverUtils::ForcingSharedPtr> m_forcing;

    Array<OneD, NekDouble> &GetAdvVelNormElec();

    // Diffusive Flux vector
    void GetFluxVectorDiff(
        const Array<OneD, Array<OneD, NekDouble>> &in_arr,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);

    // Advective Flux vector
    void GetFluxVectorElec(
        const Array<OneD, Array<OneD, NekDouble>> &fields_vals,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);
};

} // namespace NESO::Solvers::tokamak
#endif
