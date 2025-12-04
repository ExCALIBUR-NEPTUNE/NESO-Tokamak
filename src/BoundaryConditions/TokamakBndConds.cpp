#include "TokamakBndConds.hpp"
#include "../EquationSystems/TokamakSystem.hpp"

using namespace std;

namespace NESO::Solvers::tokamak
{

std::set<std::string> TokamakBoundaryConditions::m_BndType = {
    "Sheath", "Oblique", "ObliqueOutflow", "Bohm"};

TokamakBoundaryConditions::TokamakBoundaryConditions()
{
}

/**
 *  Initialize the partial slip boundary conditions
 *  The total points, unit normal vector, and the slip length are assigned
 */
void TokamakBoundaryConditions::Initialize(
    const LibUtilities::SessionReaderSharedPtr pSession,
    const std::weak_ptr<TokamakSystem>& pSystem,
    Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
    const Array<OneD, MR::DisContFieldSharedPtr> &pB,
    const Array<OneD, MR::DisContFieldSharedPtr> &pE, int pSpacedim)
{
    int nvariables = pFields.size();
    Array<OneD, Array<OneD, const SpatialDomains::BoundaryConditionShPtr>>
        BndConds(nvariables);
    Array<OneD, Array<OneD, MultiRegions::ExpListSharedPtr>> BndExp(nvariables);

    for (int i = 0; i < nvariables; ++i)
    {
        BndConds[i] = pFields[i]->GetBndConditions();
        BndExp[i]   = pFields[i]->GetBndCondExpansions();
    }

    for (size_t n = 0; n < BndExp[0].size(); ++n)
    {
        Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond(
            BndConds.size());
        Array<OneD, MultiRegions::ExpListSharedPtr> exp(BndConds.size());
        std::string bndtype;
        for (int k = 0; k < BndConds.size(); ++k)
        {
            cond[k] = BndConds[k][n];
            exp[k]  = BndExp[k][n];
            if (bndtype.size() == 0 &&
                m_BndType.find(cond[k]->GetUserDefined()) != m_BndType.end())
            {
                bndtype = cond[k]->GetUserDefined();
            }
        }
        if (bndtype.size())
        {
            m_bounds[n] = GetTokamakBaseBndCondFactory().CreateInstance(
                bndtype, pSession, pSystem, pFields, pB, pE, cond, exp,
                pSpacedim, n);
        }
    }
}

void TokamakBoundaryConditions::Update(
    const Array<OneD, const Array<OneD, NekDouble>> &physarray, NekDouble time)
{
    for (auto &it : m_bounds)
    {
        it.second->Apply(physarray, time);
    }
}

} // namespace NESO::Solvers::tokamak
