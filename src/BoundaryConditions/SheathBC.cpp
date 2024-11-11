#include "SheathBC.hpp"

using namespace std;

namespace NESO::Solvers::tokamak
{

std::string SheathBC::className =
    GetTokamakBndCondFactory().RegisterCreatorFunction(
        "Sheath", SheathBC::create, "Sheath boundary condition.");

SheathBC::SheathBC(const LU::SessionReaderSharedPtr &pSession,
                   const Array<OneD, MR::ExpListSharedPtr> &pFields,
                   const Array<OneD, Array<OneD, NekDouble>> &pObliqueField,
                   const int pSpaceDim, const int bcRegion)
    : TokamakBndCond(pSession, pFields, pObliqueField, pSpaceDim, bcRegion)
{
    m_session->LoadParameter("Ge", Ge, 0.0);
    m_session->LoadParameter("lambda", lambda, 0.0);
    m_session->LoadParameter("m_e", Me, 0.0);
    m_session->LoadParameter("m_i", Mi, 0.0);
    m_session->LoadParameter("Z", Zi, 0.0);

    int nBCEdgePts =
        m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetTotPoints();

    // Take dot product of unit b with normals to
    // get wall angle
    sin_alpha = Array<OneD, NekDouble>(nBCEdgePts, 0.0);
    for (int d = 0; d < m_spacedim; ++d)
    {
        // B.n = Bncos(pi-alpha)=sin(alpha)
        Vmath::Vvtvp(nBCEdgePts, m_b[d], 1, m_normals[d], 1, sin_alpha, 1,
                     sin_alpha, 1);
    }
}

void SheathBC::v_Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                       [[maybe_unused]] const NekDouble &time)
{
    int nBCEdgePts = m_bndExp[0]->GetTotPoints();
    int Te_idx     = 1;
    int ne_idx     = 2;
    int Ti_idx     = 3;
    int ni_idx     = 4;

    // Get electron density gradient
    Array<OneD, NekDouble> ne_bnd;
    m_fields[0]->ExtractPhysToBndElmt(m_bcRegion, physarray[ne_idx], ne_bnd);
    Array<OneD, Array<OneD, NekDouble>> ne_grad(m_spacedim);
    for (int d = 0; d < m_spacedim; d++)
    {
        ne_grad[d] = Array<OneD, NekDouble>(ne_bnd.size(), 0.0);
    }
    if (m_spacedim == 2)
    {
        m_bndElmtExp->PhysDeriv(ne_bnd, ne_grad[0], ne_grad[1]);
    }
    else if (m_spacedim == 3)
    {
        m_bndElmtExp->PhysDeriv(ne_bnd, ne_grad[0], ne_grad[1], ne_grad[2]);
    }
    Array<OneD, Array<OneD, NekDouble>> bndGrad(m_spacedim);
    for (int k = 0; k < m_spacedim; k++)
    {
        bndGrad[k] = Array<OneD, NekDouble>(nBCEdgePts, 0.0);
        m_fields[0]->ExtractElmtToBndPhys(m_bcRegion, ne_grad[k], bndGrad[k]);
    }
    Array<OneD, NekDouble> grad_ne(nBCEdgePts, 0.0);
    for (int d = 0; d < m_spacedim; ++d)
    {
        Vmath::Vvtvp(nBCEdgePts, bndGrad[d], 1, m_normals[d], 1, grad_ne, 1,
                     grad_ne, 1);
    }

    // // Get ion density gradients and calculate ion sum
    Array<OneD, NekDouble> ion_sum(nBCEdgePts, 0.0);
    // Add support for multiple ions later
    int n_species = 1;
    for (int s = 0; s < n_species; ++s)
    {
        Array<OneD, NekDouble> ni_bnd;

        m_fields[0]->ExtractPhysToBndElmt(m_bcRegion, physarray[ni_idx],
                                          ni_bnd);

        // Obtain field derivative
        Array<OneD, Array<OneD, NekDouble>> ni_grad(m_spacedim);
        for (int d = 0; d < m_spacedim; d++)
        {
            ni_grad[d] = Array<OneD, NekDouble>(ni_bnd.size(), 0.0);
        }
        if (m_spacedim == 2)
        {
            m_bndElmtExp->PhysDeriv(ni_bnd, ni_grad[0], ni_grad[1]);
        }
        else if (m_spacedim == 3)
        {
            m_bndElmtExp->PhysDeriv(ni_bnd, ni_grad[0], ni_grad[1], ni_grad[2]);
        }
        Array<OneD, Array<OneD, NekDouble>> bndGrad(m_spacedim);
        for (int k = 0; k < m_spacedim; k++)
        {
            bndGrad[k] = Array<OneD, NekDouble>(nBCEdgePts, 0.0);
            m_fields[0]->ExtractElmtToBndPhys(m_bcRegion, ni_grad[k],
                                              bndGrad[k]);
        }

        Array<OneD, NekDouble> grad_ni(nBCEdgePts, 0.0);
        for (int d = 0; d < m_spacedim; ++d)
        {
            Vmath::Vvtvp(nBCEdgePts, bndGrad[d], 1, m_normals[d], 1, grad_ni, 1,
                         grad_ni, 1);
        }

        Array<OneD, NekDouble> s_i(nBCEdgePts, 0.0);
        Vmath::Vdiv(nBCEdgePts, m_bndExp[ni_idx]->GetPhys(), 1,
                    m_bndExp[ne_idx]->GetPhys(), 1, s_i, 1);

        for (int i = 0; i < nBCEdgePts; ++i)
        {
            NekDouble C_i_sq = (adiabatic * m_bndExp[Ti_idx]->GetPhys()[i] +
                                Zi * s_i[i] * m_bndExp[Te_idx]->GetPhys()[i] *
                                    grad_ne[i] / grad_ni[i]) /
                               Mi;
            ion_sum[i] += s_i[i] * Zi * sin_alpha[i] * std::sqrt(C_i_sq);
        }
    }

    Array<OneD, NekDouble> phi(nBCEdgePts, 0.0);
    for (int i = 0; i < nBCEdgePts; ++i)
    {
        phi[i] = m_bndExp[Te_idx]->GetPhys()[i] *
                 log(sqrt(m_bndExp[Te_idx]->GetPhys()[i] / (Me * 2 * M_PI)) *
                     (1. - Ge) / ion_sum[i]);
    }
    Array<OneD, NekDouble> wall_potential(nBCEdgePts, lambda);
    Vmath::Vadd(nBCEdgePts, wall_potential, 1, phi, 1, phi, 1);

    m_bndExp[0]->FwdTrans(phi, m_bndExp[0]->UpdateCoeffs());
}

} // namespace NESO::Solvers::tokamak
