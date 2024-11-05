#ifndef TOKAMAK_OBLIQUE_HPP
#define TOKAMAK_OBLIQUE_HPP

#include "TokamakBndCond.hpp"

namespace NESO::Solvers::tokamak
{

class ObliqueBC : public TokamakBndCond
{
public:
    friend class MemoryManager<ObliqueBC>;

    static TokamakBndCondSharedPtr create(
        const LU::SessionReaderSharedPtr &pSession,
        const Array<OneD, MR::ExpListSharedPtr> &pFields,
        const Array<OneD, Array<OneD, NekDouble>> &pMagneticField,
        const int pSpaceDim, const int bcRegion, const int index)

    {
        TokamakBndCondSharedPtr p = MemoryManager<ObliqueBC>::AllocateSharedPtr(
            pSession, pFields, pMagneticField, pSpaceDim, bcRegion, index);
        return p;
    }

    static std::string className;

protected:
    void v_Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                 const NekDouble &time) override;

private:
    ObliqueBC(const LU::SessionReaderSharedPtr &pSession,
              const Array<OneD, MR::ExpListSharedPtr> &pFields,
              const Array<OneD, Array<OneD, NekDouble>> &pObliqueFields,
              const int pSpaceDim, const int bcRegion, const int index);
    ~ObliqueBC(void) override {};

    void CalcKPar();
    void CalcKPerp();
    void CalcDTensor();
    void AddRHS();
    Array<OneD, NekDouble> m_D[3][3];
    Array<OneD, NekDouble> kpar;
    Array<OneD, NekDouble> kperp;
};

} // namespace NESO::Solvers::tokamak
#endif