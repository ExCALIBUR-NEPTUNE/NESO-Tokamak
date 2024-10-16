#include "SheathBC.h"

using namespace std;

namespace NESO::Solvers::tokamak
{

std::string SheathBC::className =
    GetTokamakBndCondFactory().RegisterCreatorFunction(
        "Sheath", SheathBC::create, "Sheath boundary condition.");

SheathBC::SheathBC(const LU::SessionReaderSharedPtr &pSession,
                     const Array<OneD, MR::ExpListSharedPtr> &pFields,
                     const Array<OneD, Array<OneD, NekDouble>> &pObliqueField,
                     const int pSpaceDim, const int bcRegion, const int cnt)
    : TokamakBndCond(pSession, pFields, pObliqueField, pSpaceDim, bcRegion, cnt)
{
    int nBCEdgePts =
        m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetTotPoints();
}

void SheathBC::v_Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                        [[maybe_unused]] const NekDouble &time)
{
     int nBCEdgePts = m_bndExp[0]->GetTotPoints();
}

} // namespace NESO::Solvers::tokamak
