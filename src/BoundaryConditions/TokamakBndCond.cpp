#include "TokamakBndCond.h"

using namespace std;

namespace NESO::Solvers::tokamak
{
TokamakBndCondFactory &GetTokamakBndCondFactory()
{
    static TokamakBndCondFactory instance;
    return instance;
}

TokamakBndCond::TokamakBndCond(
    const LU::SessionReaderSharedPtr &pSession,
    const Array<OneD, MR::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble>> &pMagneticField,
    const int pSpaceDim, const int bcRegion, const int cnt)
    : m_session(pSession), m_fields(pFields),
      m_magneticFieldTrace(pMagneticField), m_spacedim(pSpaceDim),
      m_bcRegion(bcRegion), m_offset(cnt)
{
    m_bndExp = Array<OneD, MultiRegions::ExpListSharedPtr>(m_fields.size());
    for (int v = 0; v < m_fields.size(); ++v)
    {
        m_bndExp[v] = m_fields[0]->GetBndCondExpansions()[m_bcRegion];
    }
    m_fields[0]->GetBndElmtExpansion(m_bcRegion, m_bndElmtExp, false);

    m_fields[0]->GetBoundaryNormals(m_bcRegion, m_normals);
    m_diffusionAveWeight = 1.0;
}

/**
 * @param   bcRegion      id of the boundary region
 * @param   cnt
 * @param   Fwd
 * @param   physarray
 * @param   time
 */
void TokamakBndCond::Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                           const NekDouble &time)
{
    v_Apply(physarray, time);
}

/**
 * @ brief Newly added bc should specify this virtual function
 * if the Bwd/value in m_bndCondExpansions is the target value like Direchlet
 * bc weight should be 1.0.
 * if some average Fwd and Bwd/value in m_bndCondExpansions
 * is the target value like WallViscousBC weight should be 0.5.
 */
void TokamakBndCond::v_ApplyBwdWeight()
{
    size_t nVariables = m_fields.size();
    for (int i = 0; i < nVariables; ++i)
    {
        m_fields[i]->SetBndCondBwdWeight(m_bcRegion, m_diffusionAveWeight);
    }
}

} // namespace NESO::Solvers::tokamak
