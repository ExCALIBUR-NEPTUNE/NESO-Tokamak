#ifndef TOKAMAK_OBLIQUE_H
#define TOKAMAK_OBLIQUE_H

#include "TokamakBndCond.h"

namespace NESO::Solvers::tokamak
{

class ObliqueBC : public TokamakBndCond
{
public:
    friend class MemoryManager<ObliqueBC>;

    static TokamakBndCondSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
        const Array<OneD, Array<OneD, NekDouble>> &pObliqueField,
        const int pSpaceDim, const int bcRegion, const int cnt)

    {
        TokamakBndCondSharedPtr p = MemoryManager<ObliqueBC>::AllocateSharedPtr(
            pSession, pFields, pTraceNormals, pObliqueField, pSpaceDim,
            bcRegion, cnt);
        return p;
    }

    static std::string className;

protected:
    void v_Apply(Array<OneD, Array<OneD, NekDouble>> &magnetic,
                 Array<OneD, Array<OneD, NekDouble>> &physarray,
                 const NekDouble &time) override;

private:
    ObliqueBC(const LibUtilities::SessionReaderSharedPtr &pSession,
              const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
              const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
              const Array<OneD, Array<OneD, NekDouble>> &pObliqueFields,
              const int pSpaceDim, const int bcRegion, const int cnt);
    ~ObliqueBC(void) override {};
};

} // namespace NESO::Solvers::tokamak
#endif