#ifndef SINGLEDIFFUSIVEFIELD_HPP
#define SINGLEDIFFUSIVEFIELD_HPP
#include "TokamakSystem.hpp"

namespace NESO::Solvers::tokamak
{
class ElectrostaticTurbulence : public TokamakSystem
{
public:
    friend class MemoryManager<ElectrostaticTurbulence>;

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
            MemoryManager<ElectrostaticTurbulence>::AllocateSharedPtr(session,
                                                                      graph);
        p->InitObject();
        return p;
    }
    ~ElectrostaticTurbulence() override = default;

protected:
    ElectrostaticTurbulence(const LU::SessionReaderSharedPtr &session,
                            const SD::MeshGraphSharedPtr &graph);
    void v_InitObject(bool DeclareFields = true) override;
    void v_SetInitialConditions(NekDouble init_time, bool dump_ICs,
                                const int domain) override;
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                  Array<OneD, Array<OneD, NekDouble>> &out_arr,
                  const NekDouble time);
    void CalcKPar();
    void CalcKPerp();
    void CalcDiffTensor();
    void CalcKappaPar();
    void CalcKappaPerp();
    void CalcKappaTensor();

    void CalcInitPhi();
    void SolvePhi(const Array<OneD, const Array<OneD, NekDouble>> &in_arr);
    void ComputeGradPhi();

    Array<OneD, NekDouble> &GetAdvVelNorm(
        Array<OneD, NekDouble> &trace_vel_norm,
        const Array<OneD, Array<OneD, NekDouble>> &adv_vel);

    Array<OneD, NekDouble> &GetAdvVelNormElec();

    // Advective Flux vector
    void GetFluxVectorElec(
        const Array<OneD, Array<OneD, NekDouble>> &fields_vals,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);

    void GetFluxVector(const Array<OneD, Array<OneD, NekDouble>> &fields_vals,
                       const Array<OneD, Array<OneD, NekDouble>> &adv_vel,
                       Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &flux);
    // Diffusive Flux vector
    void GetFluxVectorDiff(
        const Array<OneD, Array<OneD, NekDouble>> &in_arr,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);

    void load_params() override;

private:
    /// Hasegawa-Wakatani α
    NekDouble alpha;
    /// Hasegawa-Wakatani κ
    NekDouble kappa;
    /// Potential Gradient/-E
    Array<OneD, MR::DisContFieldSharedPtr> m_grad_phi;
    /// Storage for ExB drift velocity
    Array<OneD, Array<OneD, NekDouble>> ExB_vel;
    StdRegions::VarCoeffMap m_phi_varcoeff;
    /// Storage for component of ne advection velocity normal to trace elements
    Array<OneD, NekDouble> norm_vel_elec;

    // For Diffusion
    StdRegions::ConstFactorMap m_factors;

    Array<OneD, NekDouble> m_kperp;
    Array<OneD, NekDouble> m_kpar;
    StdRegions::VarCoeffMap m_D;

    Array<OneD, NekDouble> m_kappaperp;
    Array<OneD, NekDouble> m_kappapar;
    StdRegions::VarCoeffMap m_kappa;
};

} // namespace NESO::Solvers::tokamak
#endif