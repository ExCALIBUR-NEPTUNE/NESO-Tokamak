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
                     const int pSpaceDim, const int bcRegion, const int index)
    : TokamakBndCond(pSession, pFields, pObliqueField, pSpaceDim, bcRegion,
                     index)
{
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
    // Obtain field on boundary elmts
    Array<OneD, NekDouble> Fwd;
    m_fields[m_index]->ExtractPhysToBndElmt(m_bcRegion, physarray[m_index],
                                            Fwd);

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
    Array<OneD, Array<OneD, NekDouble>> bndGrad(m_spacedim);
    for (int k = 0; k < m_spacedim; k++)
    {
        bndGrad[k] = Array<OneD, NekDouble>(m_nEdgePts, 0.0);
        m_fields[m_index]->ExtractElmtToBndPhys(m_bcRegion, grad[k],
                                                bndGrad[k]);
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
            Vmath::Vvtvp(m_nEdgePts, m_D[k][j], 1, bndGrad[j], 1, Flux[k], 1,
                         Flux[k], 1);
        }
        Vmath::Vvtvp(m_nEdgePts, Flux[k], 1, m_normals[k], 1, Flux_n, 1, Flux_n,
                     1);
    }

    m_bndExp[m_index]->IProductWRTBase(Flux_n,
                                       m_bndExp[m_index]->UpdateCoeffs());
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
}

} // namespace NESO::Solvers::tokamak
