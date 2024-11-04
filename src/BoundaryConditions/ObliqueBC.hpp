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
    ObliqueBC(const LU::SessionReaderSharedPtr &pSession,
              const Array<OneD, MR::ExpListSharedPtr> &pFields,
              const Array<OneD, Array<OneD, NekDouble>> &pObliqueFields,
              const int pSpaceDim, const int bcRegion, const int cnt);
    ~ObliqueBC(void) override {};

    Array<OneD, NekDouble> m_diffusivity[3][3];
    Array<OneD, NekDouble> m_kpar;
    Array<OneD, NekDouble> m_kperp;
};

} // namespace NESO::Solvers::tokamak
#endif