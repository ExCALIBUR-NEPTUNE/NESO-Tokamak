#include "ObliqueBC.h"

using namespace std;

namespace NESO::Solvers::tokamak
{

std::string ObliqueBC::className =
    GetTokamakBndCondFactory().RegisterCreatorFunction(
        "Oblique", ObliqueBC::create, "Oblique boundary condition.");

ObliqueBC::ObliqueBC(const LU::SessionReaderSharedPtr &pSession,
                     const Array<OneD, MR::ExpListSharedPtr> &pFields,
                     const Array<OneD, Array<OneD, NekDouble>> &pNormals,
                     const Array<OneD, Array<OneD, NekDouble>> &pObliqueField,
                     const int pSpaceDim, const int bcRegion, const int cnt)
    : TokamakBndCond(pSession, pFields, pNormals, pObliqueField, pSpaceDim,
                     bcRegion, cnt)
{
}

void ObliqueBC::v_Apply(Array<OneD, Array<OneD, NekDouble>> &magnetic,
                        Array<OneD, Array<OneD, NekDouble>> &physarray,
                        [[maybe_unused]] const NekDouble &time)
{
    // Obtain boundary normals
    Array<OneD, Array<OneD, NekDouble>> Normals(m_spacedim);
    m_fields[0]->GetBoundaryNormals(m_bcRegion, Normals);

    int nBCEdgePts =
        m_fields[0]->GetBndCondExpansions()[m_bcRegion]->GetTotPoints();

    // Obtain magnetic field at boundary
    Array<OneD, Array<OneD, NekDouble>> Bndmagnetic(m_spacedim);

    for (int i = 0; i < m_spacedim; ++i)
    {
        Bndmagnetic[i] = Array<OneD, NekDouble>(nBCEdgePts, 0.0);
        m_fields[0]->ExtractElmtToBndPhys(m_bcRegion, magnetic[i],
                                          Bndmagnetic[i]);
    }

    for (int i = 0; i < physarray.size(); ++i)
    {
        // Obtain field derivative
        Array<OneD, NekDouble> Deriv(m_spacedim * physarray[i].size(), 0.0);
        m_fields[0]->GetBndCondExpansions()[m_bcRegion]->PhysDeriv(physarray[i],
                                                                   Deriv);

        // Obtain field derivative at boundary
        Array<OneD, NekDouble> BndDeriv(m_spacedim * nBCEdgePts, 0.0);
        m_fields[0]->ExtractElmtToBndPhys(m_bcRegion, Deriv, BndDeriv);

        // Calculate grad(T).B
        Array<OneD, NekDouble> Result(nBCEdgePts, 0.0);
        Array<OneD, NekDouble> NormDeriv(nBCEdgePts, 0.0);

        for (int j = 0; j < m_spacedim; ++j)
        {
            Vmath::Vvtvp(nBCEdgePts, &BndDeriv[j * nBCEdgePts], 1,
                         &Bndmagnetic[j][0], 1, &Result[0], 1, &Result[0], 1);
            Vmath::Vabs(nBCEdgePts, &Result[0], 1, &Result[0], 1);
        }
        Vmath::Smul(nBCEdgePts, -1.0, &Result[0], 1, &Result[0], 1);

        // Copy boundary adjusted values into the boundary expansion
        int nbcoeffs =
            m_fields[i]->GetBndCondExpansions()[m_bcRegion]->GetNcoeffs();
        Array<OneD, NekDouble> bndCoeffs(nbcoeffs, 0.0);

        m_fields[i]->GetBndCondExpansions()[m_bcRegion]->IProductWRTBase(
            Result, bndCoeffs);

        Vmath::Vadd(
            nbcoeffs, &bndCoeffs[0], 1,
            &(m_fields[i]->GetBndCondExpansions()[m_bcRegion]->GetCoeffs())[0],
            1,
            &(m_fields[i]
                  ->GetBndCondExpansions()[m_bcRegion]
                  ->UpdateCoeffs())[0],
            1);
    }
}

} // namespace NESO::Solvers::tokamak
