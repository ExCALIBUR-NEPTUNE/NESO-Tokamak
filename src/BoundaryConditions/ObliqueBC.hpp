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
        const int pSpaceDim, const int bcRegion)
    {
        TokamakBndCondSharedPtr p = MemoryManager<ObliqueBC>::AllocateSharedPtr(
            pSession, pFields, pMagneticField, pSpaceDim, bcRegion);
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
              const int pSpaceDim, const int bcRegion);
    ~ObliqueBC() override {};

    void CalcKPar();
    void CalcKPerp();
    void CalcDTensor();

    Array<OneD, NekDouble> m_D[3][3];
    Array<OneD, NekDouble> kpar;
    Array<OneD, NekDouble> kperp;
    NekDouble T_bnd;
    NekDouble m_i;
    NekDouble k_B;
    NekDouble k_par;
    NekDouble k_perp;
};

} // namespace NESO::Solvers::tokamak
#endif