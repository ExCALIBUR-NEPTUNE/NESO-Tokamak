#ifndef TOKAMAK_BOHM_HPP
#define TOKAMAK_BOHM_HPP

#include "TokamakBaseBndCond.hpp"

namespace NESO::Solvers::tokamak
{

class BohmBC : public TokamakBaseBndCond
{
public:
    friend class MemoryManager<BohmBC>;

    static TokamakBaseBndCondSharedPtr create(
        const LU::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<TokamakSystem> &pSystem,
        const Array<OneD, MR::ExpListSharedPtr> &pFields,
        const Array<OneD, MR::DisContFieldSharedPtr> &pB,
        const Array<OneD, MR::DisContFieldSharedPtr> &pE,
        Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
        Array<OneD, MultiRegions::ExpListSharedPtr> exp, const int pSpaceDim,
        const int bcRegion)

    {
        TokamakBaseBndCondSharedPtr p =
            MemoryManager<BohmBC>::AllocateSharedPtr(pSession, pSystem, pFields,
                                                     pB, cond, exp, pE,
                                                     pSpaceDim, bcRegion);
        return p;
    }

    static std::string className;

protected:
    void v_Apply(const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
                 const Array<OneD, const Array<OneD, NekDouble>> &physarray,
                 const NekDouble &time) override;

private:
    BohmBC(const LU::SessionReaderSharedPtr &pSession,
           const std::weak_ptr<TokamakSystem> &pSystem,
           const Array<OneD, MR::ExpListSharedPtr> &pFields,
           const Array<OneD, MR::DisContFieldSharedPtr> &pB,
           const Array<OneD, MR::DisContFieldSharedPtr> &pE,
           Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
           Array<OneD, MultiRegions::ExpListSharedPtr> exp, const int pSpaceDim,
           const int bcRegion);
    ~BohmBC() override {};

    // Hardcoded for now
    NekDouble gamma_i = 5 / 2;
    NekDouble gamma_e = 9 / 2;

    Array<OneD, Array<OneD, NekDouble>> v_ExB;

    NekDouble lambda;
    NekDouble Ge;
};

} // namespace NESO::Solvers::tokamak
#endif