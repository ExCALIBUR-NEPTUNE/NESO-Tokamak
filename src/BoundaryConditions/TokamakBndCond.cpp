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
    const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
    const Array<OneD, Array<OneD, NekDouble>> &pObliqueField,
    const int pSpaceDim, const int bcRegion, const int cnt)
    : m_session(pSession), m_fields(pFields), m_normals(pTraceNormals),
      m_obliqueField(pObliqueField), m_spacedim(pSpaceDim),
      m_bcRegion(bcRegion), m_offset(cnt)
{
    m_diffusionAveWeight = 1.0;
}

/**
 * @param   bcRegion      id of the boundary region
 * @param   cnt
 * @param   Fwd
 * @param   physarray
 * @param   time
 */
void TokamakBndCond::Apply(Array<OneD, Array<OneD, NekDouble>> &magnetic,
                           Array<OneD, Array<OneD, NekDouble>> &physarray,
                           const NekDouble &time)
{
    v_Apply(magnetic, physarray, time);
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
