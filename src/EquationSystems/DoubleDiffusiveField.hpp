#ifndef DOUBLEDIFFUSIVEFIELD_HPP
#define DOUBLEDIFFUSIVEFIELD_HPP
#include "TokamakSystem.hpp"
#include "../Misc/VariableConverter.hpp"

namespace NESO::Solvers::tokamak
{
class DoubleDiffusiveField : public TokamakSystem
{
public:
    friend class MemoryManager<DoubleDiffusiveField>;

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
            MemoryManager<DoubleDiffusiveField>::AllocateSharedPtr(session,
                                                                   graph);
        p->InitObject();
        return p;
    }
    ~DoubleDiffusiveField() override = default;

protected:
    DoubleDiffusiveField(const LU::SessionReaderSharedPtr &session,
                         const SD::MeshGraphSharedPtr &graph);
    void v_InitObject(bool DeclareFields = true) override;

    void ImplicitTimeIntCG(
        const Array<OneD, const Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
        const NekDouble lambda);

    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                  Array<OneD, Array<OneD, NekDouble>> &out_arr,
                  const NekDouble time);

    void CalcKPar();
    void CalcKPerp();
    void CalcDiffTensor();

    void DoDiffusion(const Array<OneD, Array<OneD, NekDouble>> &inarray,
                     Array<OneD, Array<OneD, NekDouble>> &outarray,
                     const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                     const Array<OneD, Array<OneD, NekDouble>> &pBwd);
    // Diffusive Flux vector
    void GetFluxVectorDiff(
        const Array<OneD, Array<OneD, NekDouble>> &in_arr,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);

    void load_params() override;

    void v_ExtraFldOutput(std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
                          std::vector<std::string> &variables) override;

    // For Diffusion
    StdRegions::ConstFactorMap m_factors;
    NekDouble m_epsilon;

    bool m_useSpecVanVisc;
    NekDouble
        m_sVVCutoffRatio; // Cut-off ratio from which to start decaying modes
    NekDouble m_sVVDiffCoeff; // Diffusion coefficient of SVV modes

    NekDouble k_par;
    NekDouble k_perp;


    Array<OneD, NekDouble> m_kperp;
    Array<OneD, NekDouble> m_kpar;
    StdRegions::VarCoeffMap m_D;
    // std::vector<StdRegions::VarCoeffMap> m_D;

    NekDouble m_k_B;
};

} // namespace NESO::Solvers::tokamak
#endif