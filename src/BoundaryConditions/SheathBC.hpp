#ifndef TOKAMAK_SHEATH_HPP
#define TOKAMAK_SHEATH_HPP

#include "TokamakBndCond.hpp"

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
        const int pSpaceDim, const int bcRegion, const int index)

    {
        TokamakBndCondSharedPtr p = MemoryManager<SheathBC>::AllocateSharedPtr(
            pSession, pFields, pMagneticField, pSpaceDim, bcRegion, index);
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
              const int pSpaceDim, const int bcRegion, const int index);
    ~SheathBC(void) override {};

    // Hardcoded for now
    NekDouble adiabatic = 5/3;

    NekDouble lambda;
    NekDouble Ge;
    NekDouble Me;
    NekDouble Mi;
    NekDouble Zi;

    Array<OneD, NekDouble> sin_alpha;
};

} // namespace NESO::Solvers::tokamak
#endif