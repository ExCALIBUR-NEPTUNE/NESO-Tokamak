#ifndef ELECTROSTATICTURBULENCE_HPP
#define ELECTROSTATICTURBULENCE_HPP
#include "../Misc/VariableConverter.hpp"
#include "TokamakSystem.hpp"

namespace NESO::Solvers::tokamak
{
// [omega, {n, mnv, 1.5p}_e, {n, mnv, 1.5p}_i, ...]
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

    /// Advection functions
    void DoAdvection(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray,
                     const NekDouble time,
                     const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                     const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    void DoParticles(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &out_arr);

    void DoOdeImplicitRhs(
        const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
        Array<OneD, Array<OneD, NekDouble>> &out_arr, const NekDouble time);

    void DoOdeRhsCoeff(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                       Array<OneD, Array<OneD, NekDouble>> &out_arr,
                       const NekDouble time);
    void DoAdvectionCoeff(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                          Array<OneD, Array<OneD, NekDouble>> &outarray,
                          const NekDouble time,
                          const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                          const Array<OneD, Array<OneD, NekDouble>> &pBwd);
    void DoDiffusionCoeff(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const Array<OneD, const Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, const Array<OneD, NekDouble>> &pBwd);

    void CalcInitPhi();
    void SolvePhi(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                  [[maybe_unused]] const Array<OneD, NekDouble> &ne);
    void ComputeE();

    void CalcVelocities(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                        const Array<OneD, NekDouble> &ne);

    Array<OneD, Array<OneD, NekDouble>> &GetAdvVelNorm();

    // Advective Flux vector
    void GetFluxVector(
        const Array<OneD, Array<OneD, NekDouble>> &field_vals,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);

    /// Diffusion functions
    void DoDiffusion(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray,
                     const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                     const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    void CalcKPar();
    void CalcKPerp();
    void CalcDiffTensor();
    void CalcKappaPar();
    void CalcKappaPerp();
    void CalcKappaTensor();
    // Diffusive Flux vector
    void GetFluxVectorDiff(
        const Array<OneD, Array<OneD, NekDouble>> &in_arr,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);

    void load_params() override;

    void v_ExtraFldOutput(std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
                          std::vector<std::string> &variables) override;

private:
    int omega_idx;
    int pe_idx;

    std::vector<int> ni_idx;
    std::vector<int> vi_idx;
    std::vector<int> pi_idx;

    std::vector<int> ni_src_idx;
    std::vector<int> vi_src_idx;
    std::vector<int> pi_src_idx;

    /// Hasegawa-Wakatani α
    NekDouble alpha;
    /// Hasegawa-Wakatani κ
    NekDouble kappa;

    // Electric Potential
    MR::ContFieldSharedPtr phi;

    /// Velocities
    /// Storage for ExB drift velocity
    Array<OneD, Array<OneD, NekDouble>> v_ExB;
    // Electron parallel velocity
    Array<OneD, NekDouble> v_e_par;
    // Ion parallel velocities
    std::vector<Array<OneD, NekDouble>> v_i_par;
    // Electron diagmagnetic drift velocity
    Array<OneD, Array<OneD, NekDouble>> v_de;
    // Ion diagmagnetic drift velocities
    std::vector<Array<OneD, Array<OneD, NekDouble>>> v_di;
    // Per field advection velocities
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> adv_vel;
    // Per field advection velocities normal to trace elements
    Array<OneD, Array<OneD, NekDouble>> trace_vel_norm;

    StdRegions::VarCoeffMap m_phi_varcoeff;

    /// Whether the Boussinesq approximation is used for the vorticity
    bool m_boussinesq;

    // For Diffusion
    StdRegions::ConstFactorMap m_factors;

    Array<OneD, NekDouble> m_kperp;
    Array<OneD, NekDouble> m_kpar;
    StdRegions::VarCoeffMap m_D;

    Array<OneD, NekDouble> m_kappaperp;
    Array<OneD, NekDouble> m_kappapar;
    StdRegions::VarCoeffMap m_kappa;

    VariableConverterSharedPtr m_varConv;
};

} // namespace NESO::Solvers::tokamak
#endif