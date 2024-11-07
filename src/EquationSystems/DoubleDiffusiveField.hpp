#ifndef DOUBLEDIFFUSIVEFIELD_HPP
#define DOUBLEDIFFUSIVEFIELD_HPP
#include "TokamakSystem.hpp"

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
    void DoOdeRhs(const Array<OneD, const Array<OneD, NekDouble>> &in_arr,
                  Array<OneD, Array<OneD, NekDouble>> &out_arr,
                  const NekDouble time);
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