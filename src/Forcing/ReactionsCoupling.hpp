#ifndef REACTIONS_COUPLING_HPP
#define REACTIONS_COUPLING_HPP

#include "nektar_interface/utilities.hpp"
#include "nektar_interface/solver_base/time_evolved_eqnsys_base.hpp"
#include <SolverUtils/Forcing/Forcing.h>

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;

namespace NESO::Solvers::tokamak
{
class ReactionsCoupling : public SolverUtils::Forcing
{
public:
    friend class MemoryManager<ReactionsCoupling>;

    /// Creates an instance of this class
    static SolverUtils::ForcingSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem> &pEquation,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const unsigned int &pNumForcingFields, const TiXmlElement *pForce)
    {
        SolverUtils::ForcingSharedPtr p =
            MemoryManager<ReactionsCoupling>::AllocateSharedPtr(pSession,
                                                                pEquation);
        p->InitObject(pFields, pNumForcingFields, pForce);
        return p;
    }

    /// Name of the class
    static std::string className;

protected:
    void v_InitObject(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
        const unsigned int &pNumForcingFields,
        const TiXmlElement *pForce) override;

    void v_Apply(const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                 const Array<OneD, Array<OneD, NekDouble>> &inarray,
                 Array<OneD, Array<OneD, NekDouble>> &outarray,
                 const NekDouble &time) override;

private:
    NektarFieldIndexMap field_to_index;

    ReactionsCoupling(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const std::weak_ptr<SolverUtils::EquationSystem> &pEquation);
};

} // namespace NESO::Solvers::tokamak
#endif
