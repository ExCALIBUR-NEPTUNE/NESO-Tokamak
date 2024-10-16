#ifndef TOKAMAK_SHEATH_H
#define TOKAMAK_SHEATH_H

#include "TokamakBndCond.h"

namespace NESO::Solvers::tokamak
{

class SheathBC : public TokamakBndCond
{
public:
    friend class MemoryManager<SheathBC>;

    static TokamakBndCondSharedPtr create(
        const LU::SessionReaderSharedPtr &pSession,
        const Array<OneD, MR::ExpListSharedPtr> &pFields,
        const Array<OneD, Array<OneD, NekDouble>> &pMagneticField,
        const int pSpaceDim, const int bcRegion, const int cnt)

    {
        TokamakBndCondSharedPtr p = MemoryManager<ObliqueBC>::AllocateSharedPtr(
            pSession, pFields, pMagneticField, pSpaceDim, bcRegion, cnt);
        return p;
    }

    static std::string className;

protected:
    void v_Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                 const NekDouble &time) override;

private:
    SheathBC(const LU::SessionReaderSharedPtr &pSession,
              const Array<OneD, MR::ExpListSharedPtr> &pFields,
              const Array<OneD, Array<OneD, NekDouble>> &pObliqueFields,
              const int pSpaceDim, const int bcRegion, const int cnt);
    ~SheathBC(void) override {};
};

} // namespace NESO::Solvers::tokamak
#endif