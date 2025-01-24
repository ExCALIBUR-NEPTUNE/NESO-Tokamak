#include "ObliqueOutflowBC.hpp"

using namespace std;

namespace NESO::Solvers::tokamak
{

std::string ObliqueOutflowBC::className =
    GetTokamakBndCondFactory().RegisterCreatorFunction(
        "ObliqueOutflow", ObliqueOutflowBC::create,
        "ObliqueOutflow boundary condition.");

ObliqueOutflowBC::ObliqueOutflowBC(
    const LU::SessionReaderSharedPtr &pSession,
    const Array<OneD, MR::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble>> &pObliqueOutflowField,
    const int pSpaceDim, const int bcRegion)
    : TokamakBndCond(pSession, pFields, pObliqueOutflowField, pSpaceDim,
                     bcRegion)
{
    m_session->LoadParameter("gamma", gamma, 7.0);
    m_session->LoadParameter("m_i", m_i);
    m_session->LoadParameter("k_B", k_B);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            m_D[i][j]     = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
            m_kappa[i][j] = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
        }
    }
}

void ObliqueOutflowBC::v_Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                               [[maybe_unused]] const NekDouble &time)
{
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> bndGrad(m_fields.size());

    for (int i = 0; i < m_fields.size(); ++i)
    {
        // Obtain field on boundary elmts
        Array<OneD, NekDouble> Fwd;
        m_fields[i]->ExtractPhysToBndElmt(m_bcRegion, physarray[i], Fwd);

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
        // Obtain gradient on boundary of elements
        bndGrad[i] = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (int k = 0; k < m_spacedim; k++)
        {
            bndGrad[i][k] = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
            m_fields[i]->ExtractElmtToBndPhys(m_bcRegion, grad[k],
                                              bndGrad[i][k]);
        }
    }
    // Obtain anisotropic tensors
    CalcDTensor();
    CalcKappaTensor();

    // Obtain Fluxes
    Array<OneD, Array<OneD, NekDouble>> density_flux(m_spacedim);
    Array<OneD, Array<OneD, NekDouble>> energy_flux(m_spacedim);

    for (int k = 0; k < m_spacedim; k++)
    {
        density_flux[k] = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
        energy_flux[k]  = Array<OneD, NekDouble>(m_nEdgePts, 0.0);

        for (int j = 0; j < m_spacedim; j++)
        {
            Vmath::Vvtvp(m_nEdgePts, m_D[k][j], 1, bndGrad[0][j], 1,
                         density_flux[k], 1, density_flux[k], 1);
            Vmath::Vvtvp(m_nEdgePts, m_kappa[k][j], 1, bndGrad[1][j], 1,
                         energy_flux[k], 1, energy_flux[k], 1);
        }
    }
    Array<OneD, NekDouble> nbndCoeffs(m_nEdgeCoeffs, 0.0);
    Array<OneD, NekDouble> TbndCoeffs(m_nEdgeCoeffs, 0.0);

    m_bndExp[0]->NormVectorIProductWRTBase(density_flux, nbndCoeffs);
    m_bndExp[2]->NormVectorIProductWRTBase(energy_flux, TbndCoeffs);

    // Obtain rhs of bcs
    Array<OneD, NekDouble> n = m_bndExp[0]->GetPhys();
    Array<OneD, NekDouble> T = m_bndExp[2]->GetPhys();

    Array<OneD, Array<OneD, NekDouble>> rhs(2);
    for (int p = 0; p < m_nEdgePts; ++p)
    {
        NekDouble cs = -std::sqrt(k_B * T[p] / m_i);
        rhs[0][p]    = n[p] * cs;
        rhs[1][p]    = gamma * n[p] * k_B * T[p] * cs;
    }

    Array<OneD, Array<OneD, NekDouble>> wk(2);
    wk[0] = Array<OneD, NekDouble>(m_nEdgeCoeffs);
    wk[1] = Array<OneD, NekDouble>(m_nEdgeCoeffs);
    m_bndExp[0]->FwdTrans(rhs[0], wk[0]);
    m_bndExp[2]->FwdTrans(rhs[1], wk[1]);

    // Combine fluxes and rhs
    Vmath::Vvtvp(m_nEdgeCoeffs, wk[0], 1, m_bndExp[0]->UpdateCoeffs(), 1,
                 nbndCoeffs, 1, m_bndExp[0]->UpdateCoeffs(), 1);
    Vmath::Vvtvp(m_nEdgeCoeffs, wk[1], 1, m_bndExp[2]->UpdateCoeffs(), 1,
                 TbndCoeffs, 1, m_bndExp[2]->UpdateCoeffs(), 1);

    // Finally calculate pressure on boundary

    m_bndExp[0]->BwdTrans(m_bndExp[0]->GetCoeffs(), n);
    m_bndExp[2]->BwdTrans(m_bndExp[2]->GetCoeffs(), T);
    Array<OneD, NekDouble> p(m_nEdgePts);
    Vmath::Vmul(m_nEdgePts, n, 1, T, 1, p, 1);
    m_bndExp[1]->FwdTrans(p, m_bndExp[1]->UpdateCoeffs());
}

void ObliqueOutflowBC::CalcKPar()
{
    NekDouble k_par;
    m_session->LoadParameter("k_par", k_par, 100.0);
    kpar = Array<OneD, NekDouble>(m_nEdgePts, k_par);
}

void ObliqueOutflowBC::CalcKPerp()
{
    NekDouble k_perp;
    m_session->LoadParameter("k_perp", k_perp, 1.0);
    kperp = Array<OneD, NekDouble>(m_nEdgePts, k_perp);
}

void ObliqueOutflowBC::CalcDTensor()
{
    CalcKPar();
    CalcKPerp();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < m_nEdgePts; k++)
            {
                m_D[i][j][k] = (kpar[k] - kperp[k]) * m_b[i][k] * m_b[j][k];
                if (i == j)
                {
                    m_D[i][j][k] += kperp[k];
                }
            }
        }
    }
}

void ObliqueOutflowBC::CalcKappaPar()
{
    NekDouble kappa_par;
    m_session->LoadParameter("k_par", kappa_par, 100.0);
    kpar = Array<OneD, NekDouble>(m_nEdgePts, kappa_par);
}

void ObliqueOutflowBC::CalcKappaPerp()
{
    NekDouble kappa_perp;
    m_session->LoadParameter("k_perp", kappa_perp, 1.0);
    kperp = Array<OneD, NekDouble>(m_nEdgePts, kappa_perp);
}

void ObliqueOutflowBC::CalcKappaTensor()
{
    CalcKappaPar();
    CalcKappaPerp();
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < m_nEdgePts; k++)
            {
                m_kappa[i][j][k] =
                    (kappapar[k] - kappaperp[k]) * m_b[i][k] * m_b[j][k];
                if (i == j)
                {
                    m_kappa[i][j][k] += kappaperp[k];
                }
            }
        }
    }
}

} // namespace NESO::Solvers::tokamak
