#ifndef BOHM_BC_HPP
#define BOHM_BC_HPP

#include "TokamakBndCond.hpp"

namespace NESO::Solvers::tokamak
{

class BohmBC : public TokamakBndCond
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
    BohmBC(const LU::SessionReaderSharedPtr &pSession,
           const Array<OneD, MR::ExpListSharedPtr> &pFields,
           const Array<OneD, Array<OneD, NekDouble>> &pObliqueFields,
           const int pSpaceDim, const int bcRegion);
    ~BohmBC() override {};

    // indices of fields to which this BC is applied
    std::vector<int> idx;
};

} // namespace NESO::Solvers::tokamak
#endif