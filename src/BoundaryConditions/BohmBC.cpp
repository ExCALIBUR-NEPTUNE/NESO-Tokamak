#include "BohmBC.hpp"

using namespace std;

namespace NESO::Solvers::tokamak
{

std::string BohmBC::className =
    GetTokamakBndCondFactory().RegisterCreatorFunction(
        "Bohm", BohmBC::create, "Bohm-Chodura velocity boundary condition.");

BohmBC::BohmBC(const LU::SessionReaderSharedPtr &pSession,
               const Array<OneD, MR::ExpListSharedPtr> &pFields,
               const Array<OneD, Array<OneD, NekDouble>> &pObliqueField,
               const int pSpaceDim, const int bcRegion)
    : TokamakBndCond(pSession, pFields, pObliqueField, pSpaceDim, bcRegion)
{
}

void BohmBC::v_Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                     [[maybe_unused]] const NekDouble &time)
{
    // Obtain field on boundary elmts
    Array<OneD, Array<OneD, NekDouble>> FwdBnd;
    for (int i = 0; i < idx.size(); ++i)
    {
        m_fields[0]->ExtractPhysToBnd(m_bcRegion, physarray[i], FwdBnd[i]);
    }
    Array<OneD, NekDouble> vbc(m_nEdgePts);
    m_varConv->GetSystemSoundSpeed(FwdBnd, vbc);

    Array<OneD, NekDouble> ndotb(m_nEdgePts, 0.0);
    for (int d = 0; d < m_spacedim; ++d)
    {
        Vmath::Vvtvp(m_nEdgePts, m_b[d], 1, m_normals[d], 1, ndotb, 1, ndotb,
                     1);
    }

    // negate when B field points inwards
    for (int p = 0; p < m_nEdgePts; ++p)
    {
        if (ndotb[p] < 0)
        {
            vbc[p] = -vbc[p];
        }
    }

    for (int i = 0; i < idx.size(); ++i)
    {
        m_bndExp[idx[i]]->UpdatePhys() = vbc;
    }
}

} // namespace NESO::Solvers::tokamak
