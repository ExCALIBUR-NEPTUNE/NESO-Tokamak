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
                     const int pSpaceDim, const int bcRegion)
    : TokamakBndCond(pSession, pFields, pObliqueField, pSpaceDim, bcRegion)
{
    m_session->LoadParameter("T_bnd", T_bnd);
    m_session->LoadParameter("m_i", m_i);
    m_session->LoadParameter("k_B", k_B);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            m_D[i][j] = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
        }
    }
}

void ObliqueBC::v_Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                        [[maybe_unused]] const NekDouble &time)
{
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> bndGrad(m_fields.size());
    for (int i = 0; i < m_fields.size(); ++i)
    {
        // Obtain field on boundary elmts
        Array<OneD, NekDouble> Fwd;
        m_fields[0]->ExtractPhysToBndElmt(m_bcRegion, physarray[i], Fwd);

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
    CalcDTensor();
    // Obtain Flux and normal flux
    Array<OneD, Array<OneD, NekDouble>> Flux(m_spacedim);
    Array<OneD, NekDouble> Flux_n(m_nEdgePts, 0.0);
    for (int k = 0; k < m_spacedim; k++)
    {
        Flux[k] = Array<OneD, NekDouble>(m_nEdgePts, 0.0);

        for (int j = 0; j < m_spacedim; j++)
        {
            Vmath::Vvtvp(m_nEdgePts, m_D[k][j], 1, bndGrad[0][j], 1, Flux[k], 1,
                         Flux[k], 1);
        }
        Vmath::Vvtvp(m_nEdgePts, Flux[k], 1, m_normals[k], 1, Flux_n, 1, Flux_n,
                     1);
    }

    m_bndExp[0]->IProductWRTBase(Flux_n, m_bndExp[0]->UpdateCoeffs());
    AddRHS();
}

void ObliqueBC::CalcKPar()
{
    NekDouble k_par;
    m_session->LoadParameter("k_par", k_par, 100.0);
    kpar = Array<OneD, NekDouble>(m_nEdgePts, k_par);
}

void ObliqueBC::CalcKPerp()
{
    NekDouble k_perp;
    m_session->LoadParameter("k_perp", k_perp, 1.0);
    kperp = Array<OneD, NekDouble>(m_nEdgePts, k_perp);
}

void ObliqueBC::CalcDTensor()
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

void ObliqueBC::AddRHS()
{
    Array<OneD, NekDouble> n = m_bndExp[0]->GetPhys();
    NekDouble cs             = std::sqrt(k_B * T_bnd / m_i);

    Array<OneD, NekDouble> rhs(m_fields.size());
    for (int p = 0; p < m_nEdgePts; ++p)
    {
        rhs[p] = n[p] * cs;
    }
    Array<OneD, NekDouble> wk(m_nEdgeCoeffs);
    m_bndExp[0]->FwdTrans(rhs, wk);

    Vmath::Vsub(m_nEdgeCoeffs, m_bndExp[0]->GetCoeffs(), 1, wk, 1,
                m_bndExp[0]->UpdateCoeffs(), 1);
}

} // namespace NESO::Solvers::tokamak
