#ifndef SINGLEDIFFUSIVEFIELD_HPP
#define SINGLEDIFFUSIVEFIELD_HPP
#include "TokamakSystem.hpp"

namespace NESO::Solvers::tokamak
{
class SingleDiffusiveField : public TokamakSystem
{
public:
    friend class MemoryManager<SingleDiffusiveField>;

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
            MemoryManager<SingleDiffusiveField>::AllocateSharedPtr(session,
                                                                   graph);
        p->InitObject();
        return p;
    }
    ~SingleDiffusiveField() override = default;

protected:
    SingleDiffusiveField(const LU::SessionReaderSharedPtr &session,
                         const SD::MeshGraphSharedPtr &graph);
    void v_InitObject(bool DeclareFields = true) override;

    void ImplicitTimeIntCG(
    const Array<OneD, const Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble time,
    const NekDouble lambda);

    void v_ExtraFldOutput(
    std::vector<Array<OneD, NekDouble>> &fieldcoeffs,
    std::vector<std::string> &variables);
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                  Array<OneD, Array<OneD, NekDouble>> &out_arr,
                  const NekDouble time);
    void CalcKPar();
    void CalcKPerp();
    void CalcDiffTensor();

    // Diffusive Flux vector
    void GetFluxVectorDiff(
        const Array<OneD, Array<OneD, NekDouble>> &in_arr,
        const Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &q_field,
        Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &fluxes);

    void load_params() override;

    // For Diffusion
    StdRegions::ConstFactorMap m_factors;
    /// Weight for average calculation of diffusion term

    NekDouble m_k_B;
    NekDouble k_par;
    NekDouble k_perp;
    Array<OneD, NekDouble> m_kperp;
    Array<OneD, NekDouble> m_kpar;
    StdRegions::VarCoeffMap m_D;
};

} // namespace NESO::Solvers::tokamak
#endif