#ifndef TOKAMAK_OBLIQUEOUTFLOW_HPP
#define TOKAMAK_OBLIQUEOUTFLOW_HPP

#include "TokamakBndCond.hpp"

namespace NESO::Solvers::tokamak
{

class ObliqueOutflowBC : public TokamakBndCond
{
public:
    friend class MemoryManager<ObliqueOutflowBC>;

    static TokamakBndCondSharedPtr create(
        const LU::SessionReaderSharedPtr &pSession,
        const Array<OneD, MR::ExpListSharedPtr> &pFields,
        const Array<OneD, Array<OneD, NekDouble>> &pMagneticField,
        const int pSpaceDim, const int bcRegion)
    {
        TokamakBndCondSharedPtr p =
            MemoryManager<ObliqueOutflowBC>::AllocateSharedPtr(
                pSession, pFields, pMagneticField, pSpaceDim, bcRegion);
        return p;
    }

    static std::string className;

protected:
    void v_Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                 const NekDouble &time) override;

private:
    ObliqueOutflowBC(
        const LU::SessionReaderSharedPtr &pSession,
        const Array<OneD, MR::ExpListSharedPtr> &pFields,
        const Array<OneD, Array<OneD, NekDouble>> &pObliqueOutflowFields,
        const int pSpaceDim, const int bcRegion);
    ~ObliqueOutflowBC() override {};

    void CalcKPar();
    void CalcKPerp();
    void CalcDTensor();

    void CalcKappaPar();
    void CalcKappaPerp();
    void CalcKappaTensor();
    void AddRHS();
    Array<OneD, NekDouble> m_D[3][3];
    Array<OneD, NekDouble> m_kappa[3][3];
    Array<OneD, NekDouble> kpar;
    Array<OneD, NekDouble> kperp;
    Array<OneD, NekDouble> kappapar;
    Array<OneD, NekDouble> kappaperp;
    NekDouble gamma;
    NekDouble m_i;
    NekDouble k_B;
};

} // namespace NESO::Solvers::tokamak
#endif