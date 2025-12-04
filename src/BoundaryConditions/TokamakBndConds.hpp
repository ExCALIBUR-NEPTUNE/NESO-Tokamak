#ifndef TOKAMAK_BNDCONDS_HPP
#define TOKAMAK_BNDCONDS_HPP

#include "TokamakBaseBndCond.hpp"

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;

namespace NESO::Solvers::tokamak
{

class TokamakSystem;
class TokamakBoundaryConditions;
typedef std::shared_ptr<TokamakBoundaryConditions>
    IncBoundaryConditionsSharedPtr;

class TokamakBoundaryConditions
{
public:
    TokamakBoundaryConditions();

    void Initialize(const LibUtilities::SessionReaderSharedPtr pSession,
                    const std::weak_ptr<TokamakSystem> &pSystem,
                    Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
                    const Array<OneD, MR::DisContFieldSharedPtr> &pB,
                    const Array<OneD, MR::DisContFieldSharedPtr> &pE,
                    int pSpacedim);

    void Update(const Array<OneD, const Array<OneD, NekDouble>> &physarray,
                NekDouble time);

    std::map<int, TokamakBaseBndCondSharedPtr> &GetBounds()
    {
        return m_bounds;
    };

protected:
    std::map<int, TokamakBaseBndCondSharedPtr> m_bounds;
    static std::set<std::string> m_BndType;
    int m_spacedim;
    int m_bnd_dim;
};

} // namespace NESO::Solvers::tokamak
#endif