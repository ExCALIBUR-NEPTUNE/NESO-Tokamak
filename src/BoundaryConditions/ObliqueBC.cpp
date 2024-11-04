#include "ObliqueBC.hpp"

using namespace std;

namespace NESO::Solvers::tokamak
{

std::string ObliqueBC::className =
    GetTokamakBndCondFactory().RegisterCreatorFunction(
        "Oblique", ObliqueBC::create, "Oblique boundary condition.");

ObliqueBC::ObliqueBC(const LU::SessionReaderSharedPtr &pSession,
                     const Array<OneD, MR::ExpListSharedPtr> &pFields,
                     const Array<OneD, Array<OneD, NekDouble>> &pObliqueField,
                     const int pSpaceDim, const int bcRegion, const int cnt)
    : TokamakBndCond(pSession, pFields, pObliqueField, pSpaceDim, bcRegion, cnt)
{
    int nBCEdgePts =
        m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetTotPoints();
    Array<OneD, Array<OneD, NekDouble>> B(3);
    for (int i = 0; i < 3; ++i)
    {
        B[i] = Array<OneD, NekDouble>(nBCEdgePts, 0.0);
        m_fields[0]->ExtractElmtToBndPhys(m_bcRegion, m_magneticFieldTrace[i],
                                          B[i]);
        for (int j = 0; j < 3; ++j)
        {
            m_diffusivity[i][j] = Array<OneD, NekDouble>(nBCEdgePts, 0.0);
        }
    }
    NekDouble k_par, k_perp;
    m_session->LoadParameter("k_par", k_par, 100.0);
    m_session->LoadParameter("k_perp", k_perp, 1.0);
    m_kpar  = Array<OneD, NekDouble>(nBCEdgePts, k_par);
    m_kperp = Array<OneD, NekDouble>(nBCEdgePts, k_perp);

    for (int k = 0; k < nBCEdgePts; k++)
    {
        m_diffusivity[0][0][k] =
            (m_kpar[k] - m_kperp[k]) * B[0][k] * B[0][k] + m_kperp[k];
        m_diffusivity[0][1][k] = (m_kpar[k] - m_kperp[k]) * B[0][k] * B[1][k];
        m_diffusivity[0][2][k] = (m_kpar[k] - m_kperp[k]) * B[0][k] * B[2][k];
        m_diffusivity[1][0][k] = (m_kpar[k] - m_kperp[k]) * B[1][k] * B[0][k];
        m_diffusivity[1][1][k] =
            (m_kpar[k] - m_kperp[k]) * B[1][k] * B[1][k] + m_kperp[k];
        m_diffusivity[1][2][k] = (m_kpar[k] - m_kperp[k]) * B[1][k] * B[2][k];
        m_diffusivity[2][0][k] = (m_kpar[k] - m_kperp[k]) * B[2][k] * B[0][k];
        m_diffusivity[2][1][k] = (m_kpar[k] - m_kperp[k]) * B[2][k] * B[1][k];
        m_diffusivity[2][2][k] =
            (m_kpar[k] - m_kperp[k]) * B[2][k] * B[2][k] + m_kperp[k];
    }
}

void ObliqueBC::v_Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                        [[maybe_unused]] const NekDouble &time)
{
    int nBCEdgePts = m_bndExp[0]->GetTotPoints();
    for (int v = 0; v < physarray.size(); ++v)
    {
        // Obtain field on boundary elmts
        Array<OneD, NekDouble> Fwd;
        m_fields[0]->ExtractPhysToBndElmt(m_bcRegion, physarray[v], Fwd);

        // Obtain field derivative
        Array<OneD, Array<OneD, NekDouble>> grad(m_spacedim);
        for (int k = 0; k < m_spacedim; k++)
        {
            grad[k] = Array<OneD, NekDouble>(Fwd.size(), 0.0);
        }
        if (m_spacedim == 2)
        {
            m_bndElmtExp->PhysDeriv(Fwd, grad[0], grad[1]);
        }
        else if (m_spacedim == 3)
        {
            m_bndElmtExp->PhysDeriv(Fwd, grad[0], grad[1], grad[2]);
        }

        Array<OneD, Array<OneD, NekDouble>> bndGrad(m_spacedim);
        for (int k = 0; k < m_spacedim; k++)
        {
            bndGrad[k] = Array<OneD, NekDouble>(nBCEdgePts, 0.0);
            m_fields[0]->ExtractElmtToBndPhys(m_bcRegion, grad[k], bndGrad[k]);
        }
        // Obtain Flux and normal flux
        Array<OneD, Array<OneD, NekDouble>> Flux(m_spacedim);
        Array<OneD, NekDouble> Flux_n(nBCEdgePts, 0.0);
        for (int k = 0; k < m_spacedim; k++)
        {
            Flux[k] = Array<OneD, NekDouble>(nBCEdgePts, 0.0);

            for (int j = 0; j < m_spacedim; j++)
            {
                Vmath::Vvtvp(nBCEdgePts, m_diffusivity[k][j], 1, bndGrad[j], 1,
                             Flux[k], 1, Flux[k], 1);
            }
            Vmath::Vvtvp(nBCEdgePts, Flux[k], 1, m_normals[k], 1, Flux_n, 1,
                         Flux_n, 1);
        }

        m_bndExp[v]->FwdTrans(Flux_n, m_bndExp[v]->UpdateCoeffs());
    }
}

} // namespace NESO::Solvers::tokamak
