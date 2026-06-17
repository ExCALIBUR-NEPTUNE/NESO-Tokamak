///////////////////////////////////////////////////////////////////////////////
//
// File: OmegaAdvection.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: LDG diffusion class.
//
///////////////////////////////////////////////////////////////////////////////

#include <iomanip>
#include <iostream>

#include <boost/algorithm/string/predicate.hpp>

#include "OmegaAdvection.h"

namespace PENKNIFE
{

OmegaAdvection::OmegaAdvection()
{
}

void OmegaAdvection::InitObject(LibUtilities::SessionReaderSharedPtr pSession,
                                MultiRegions::ExpListSharedPtr pField)
{
    m_session = pSession;

    m_session->LoadSolverInfo("ShockCaptureType", m_shockCaptureType, "Off");

    // Set up penalty term for LDG
    m_session->LoadParameter("LDGc11", m_C11, 1.0);

    // Setting up the normals
    std::size_t nDim      = pField->GetCoordim(0);
    std::size_t nTracePts = pField->GetTrace()->GetTotPoints();

    m_traceNormals = Array<OneD, Array<OneD, NekDouble>>{nDim};
    for (std::size_t i = 0; i < nDim; ++i)
    {
        m_traceNormals[i] = Array<OneD, NekDouble>{nTracePts};
    }
    pField->GetTrace()->GetNormals(m_traceNormals);
}

void OmegaAdvection::Advect(const MultiRegions::ExpListSharedPtr &field,
                            const Array<OneD, Array<OneD, NekDouble>> &inarray,
                            Array<OneD, Array<OneD, NekDouble>> &outarray,
                            const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                            const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    std::size_t nCoeffs = field->GetNcoeffs();

    Array<OneD, NekDouble> tmp{nCoeffs};

    OmegaAdvection::AdvectCoeffs(field, inarray, tmp, pFwd, pBwd);

    field->BwdTrans(tmp, outarray[omega_idx]);
}

void OmegaAdvection::AdvectCoeffs(
    const MultiRegions::ExpListSharedPtr &field,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, NekDouble> &outarray,
    const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    if (field->GetGraph()->GetMovement()->GetMoveFlag()) // i.e. if
                                                         // m_ALESolver
    {
        field->GetTrace()->GetNormals(m_traceNormals);
    }

    std::size_t nDim      = field->GetCoordim(0);
    std::size_t nPts      = field->GetTotPoints();
    std::size_t nCoeffs   = field->GetNcoeffs();
    std::size_t nTracePts = field->GetTrace()->GetTotPoints();

    // Array<OneD, Array<OneD, NekDouble>> qfield{nDim};
    // for (std::size_t d = 0; d < nDim; ++d)
    // {
    //     qfield[d] = Array<OneD, NekDouble>{nPts, 0.0};
    // }

    // Initialize viscous tensor
    Array<OneD, Array<OneD, NekDouble>> viscTensor{nDim};
    for (std::size_t j = 0; j < nDim; ++j)
    {
        viscTensor[j] = Array<OneD, NekDouble>{nPts, 0.0};
    }
    Array<OneD, NekDouble> traceflux{nTracePts, 0.0};

    //AdvectCalcDerivative(field, inarray, viscTensor, pFwd, pBwd);
    AdvectFlux(field, inarray, viscTensor, traceflux);
    //AdvectVolumeFlux(field, inarray, viscTensor, traceflux);
    //DiffuseVolumeFlux(field, inarray, viscTensor);

    //AdvectTraceFlux(field, inarray, viscTensor, traceflux, pFwd, pBwd);
    //DiffuseTraceFlux(field, inarray, viscTensor, traceflux, pFwd, pBwd); 

    Array<OneD, Array<OneD, NekDouble>> qdbase{nDim};

    for (std::size_t j = 0; j < nDim; ++j)
    {
        qdbase[j] = viscTensor[j];
    }
    field->IProductWRTDerivBase(qdbase, outarray);

    Vmath::Neg(nCoeffs, outarray, 1);
    field->AddTraceIntegral(traceflux, outarray);
    field->SetPhysState(false);

    // If mesh is not distorted
    if (!field->GetGraph()->GetMovement()->GetMeshDistortedFlag())
    {
        field->MultiplyByElmtInvMass(outarray, outarray);
    }
}

// void OmegaAdvection::AdvectCalcDerivative(
//     const MultiRegions::ExpListSharedPtr &field,
//     const Array<OneD, Array<OneD, NekDouble>> &inarray,
//     Array<OneD, Array<OneD, NekDouble>> &qfield,
//     const Array<OneD, Array<OneD, NekDouble>> &pFwd,
//     const Array<OneD, Array<OneD, NekDouble>> &pBwd)
// {
//     std::size_t nDim      = field->GetCoordim(0);
//     std::size_t nCoeffs   = field->GetNcoeffs();
//     std::size_t nTracePts = field->GetTrace()->GetTotPoints();

//     Array<OneD, NekDouble> tmp{nCoeffs};
//     Array<OneD, Array<OneD, NekDouble>> traceflux{nDim};
//     Array<OneD, Array<OneD, NekDouble>> flux{nDim};
//     for (std::size_t j = 0; j < nDim; ++j)
//     {
//         traceflux[j] = Array<OneD, NekDouble>{nTracePts, 0.0};
//         flux[j] = Array<OneD, NekDouble>{nCoeffs, 0.0};
//     }

//     NumFluxforScalar(field, inarray[omega_idx], traceflux, pFwd, pBwd);

//     m_flux_vector(inarray, flux);
//     for (std::size_t j = 0; j < nDim; ++j)
//     {
//         //field->IProductWRTDerivBase(j, inarray[omega_idx], tmp);
//         //Vmath::Neg(nCoeffs, flux[j], 1);
//         field->AddTraceIntegral(traceflux[j], flux[j]);
//         field->SetPhysState(false);
//         field->MultiplyByElmtInvMass(flux[j], flux[j]);
//         field->BwdTrans(flux[j], qfield[j]);
//     }
// }

void OmegaAdvection::AdvectFlux(
    [[maybe_unused]] const MultiRegions::ExpListSharedPtr &field,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &viscTensor, Array<OneD, NekDouble> &traceflux)
{
    m_flux_vector(inarray, viscTensor, traceflux);
}

void OmegaAdvection::AdvectTraceFlux(
    const MultiRegions::ExpListSharedPtr &field,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &viscTensor,
    Array<OneD, NekDouble> &TraceFlux,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &pFwd,
    [[maybe_unused]] const Array<OneD, Array<OneD, NekDouble>> &pBwd)
{
    NumFluxforVector(field, inarray[omega_idx], viscTensor, TraceFlux);
}

void OmegaAdvection::NumFluxforScalar(
    const MultiRegions::ExpListSharedPtr &field,
    const Array<OneD, NekDouble> &ufield,
    Array<OneD, Array<OneD, NekDouble>> &uflux,
    const Array<OneD, Array<OneD, NekDouble>> &Fwd,
    const Array<OneD, Array<OneD, NekDouble>> &Bwd)
{
    std::size_t nTracePts = field->GetTrace()->GetTotPoints();
    std::size_t nDim      = field->GetCoordim(0);

    Array<OneD, NekDouble> fluxtemp{nTracePts, 0.0};

    // Get the sign of (v \cdot n), v = an arbitrary vector
    // Evaluate upwind flux:
    // uflux = \hat{u} \phi \cdot u = u^{(+,-)} n

    // Upwind
    Vmath::Vcopy(nTracePts, Fwd[omega_idx], 1, fluxtemp, 1);

    // Imposing weak boundary condition with flux
    if (field->GetBndCondExpansions().size())
    {
        ApplyScalarBCs(field, ufield, Fwd[omega_idx], Bwd[omega_idx], fluxtemp);
    }

    for (std::size_t j = 0; j < nDim; ++j)
    {
        Vmath::Vmul(nTracePts, m_traceNormals[j], 1, fluxtemp, 1, uflux[j], 1);
    }
}

void OmegaAdvection::ApplyScalarBCs(
    const MultiRegions::ExpListSharedPtr &field,
    [[maybe_unused]] const Array<OneD, const NekDouble> &ufield,
    const Array<OneD, const NekDouble> &Fwd,
    [[maybe_unused]] const Array<OneD, const NekDouble> &Bwd,
    Array<OneD, NekDouble> &penaltyflux)
{
    // Number of boundary regions
    std::size_t nBndRegions = field->GetBndCondExpansions().size();
    std::size_t cnt         = 0;
    for (std::size_t i = 0; i < nBndRegions; ++i)
    {
        if (field->GetBndConditions()[i]->GetBoundaryConditionType() ==
            SpatialDomains::ePeriodic)
        {
            continue;
        }

        // Number of boundary expansion related to that region
        std::size_t nBndEdges = field->GetBndCondExpansions()[i]->GetExpSize();

        // Weakly impose boundary conditions by modifying flux values
        for (std::size_t e = 0; e < nBndEdges; ++e)
        {
            std::size_t nBndEdgePts =
                field->GetBndCondExpansions()[i]->GetExp(e)->GetTotPoints();

            std::size_t id1 =
                field->GetBndCondExpansions()[i]->GetPhys_Offset(e);

            std::size_t id2 = field->GetTrace()->GetPhys_Offset(
                field->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt++));

            // AV boundary conditions
            if (boost::iequals(field->GetBndConditions()[i]->GetUserDefined(),
                               "Wall") ||
                boost::iequals(field->GetBndConditions()[i]->GetUserDefined(),
                               "Symmetry") ||
                boost::iequals(field->GetBndConditions()[i]->GetUserDefined(),
                               "WallViscous") ||
                boost::iequals(field->GetBndConditions()[i]->GetUserDefined(),
                               "WallAdiabatic") ||
                boost::iequals(field->GetBndConditions()[i]->GetUserDefined(),
                               "WallRotational"))
            {
                Vmath::Vcopy(nBndEdgePts, &Fwd[id2], 1, &penaltyflux[id2], 1);
            }
            // For Dirichlet boundary condition: uflux = g_D
            else if (field->GetBndConditions()[i]->GetBoundaryConditionType() ==
                     SpatialDomains::eDirichlet)
            {
                Vmath::Vcopy(
                    nBndEdgePts,
                    &(field->GetBndCondExpansions()[i]->GetPhys())[id1], 1,
                    &penaltyflux[id2], 1);
            }
            // For Neumann boundary condition: uflux = u+
            else if ((field->GetBndConditions()[i])
                         ->GetBoundaryConditionType() ==
                     SpatialDomains::eNeumann)
            {
                Vmath::Vcopy(nBndEdgePts, &Fwd[id2], 1, &penaltyflux[id2], 1);
            }
        }
    }
}

/**
 * @brief Build the numerical flux for the 2nd order derivatives
 * todo: add variable coeff and h dependence to penalty term
 */
void OmegaAdvection::NumFluxforVector(
    const MultiRegions::ExpListSharedPtr &field,
    const Array<OneD, NekDouble> &ufield,
    Array<OneD, Array<OneD, NekDouble>> &qfield, Array<OneD, NekDouble> &qflux)
{
    std::size_t nTracePts = field->GetTrace()->GetTotPoints();
    std::size_t nDim      = qfield.size();

    Array<OneD, NekDouble> Fwd{nTracePts};
    Array<OneD, NekDouble> Bwd{nTracePts};
    Array<OneD, NekDouble> qFwd{nTracePts};
    Array<OneD, NekDouble> qBwd{nTracePts};
    Array<OneD, NekDouble> qfluxtemp{nTracePts, 0.0};
    Array<OneD, NekDouble> uterm{nTracePts};

    Vmath::Zero(nTracePts, uterm, 1);

    // Evaulate upwind flux:
    // qflux = \hat{q} \cdot u = q \cdot n - C_(11)*(u^+ - u^-)

    // Generate Stability term = - C11 ( u- - u+ )
    field->GetFwdBwdTracePhys(ufield, Fwd, Bwd);
    Vmath::Vsub(nTracePts, Fwd, 1, Bwd, 1, uterm, 1);
    Vmath::Smul(nTracePts, m_C11, uterm, 1, uterm, 1);

    qflux = Array<OneD, NekDouble>{nTracePts, 0.0};
    for (std::size_t j = 0; j < nDim; ++j)
    {
        //  Compute Fwd and Bwd value of ufield of jth direction
        field->GetFwdBwdTracePhys(qfield[j], qFwd, qBwd);
        //field->ExtractTracePhys(qfield[j], qfluxtemp);

        // Downwind
        //Vmath::Vcopy(nTracePts, qFwd, 1, qfluxtemp, 1);
        for(size_t p = 0; p < nTracePts; ++p)
        {
            qfluxtemp[p] = 0.5 * (qFwd[p] + qBwd[p]);
        }

        Vmath::Vmul(nTracePts, m_traceNormals[j], 1, qfluxtemp, 1, qfluxtemp,
                    1);

        // Flux = {Fwd, Bwd} * (nx, ny, nz) + uterm * (nx, ny)
        // Vmath::Vadd(nTracePts, uterm, 1, qfluxtemp, 1, qfluxtemp, 1);

        // Imposing weak boundary condition with flux
        if (field->GetBndCondExpansions().size())
        {
            ApplyVectorBCs(field, j, qfield[j], qFwd, qBwd, qfluxtemp);
        }

        // q_hat \cdot n = (q_xi \cdot n_xi) or (q_eta \cdot n_eta)
        // n_xi = n_x * tan_xi_x + n_y * tan_xi_y + n_z * tan_xi_z
        // n_xi = n_x * tan_eta_x + n_y * tan_eta_y + n_z*tan_eta_z
        Vmath::Vadd(nTracePts, qfluxtemp, 1, qflux, 1, qflux, 1);
    }
}

/**
 * Diffusion: Imposing weak boundary condition for q with flux
 *  uflux = g_D  on Dirichlet boundary condition
 *  uflux = u_Fwd  on Neumann boundary condition
 */
void OmegaAdvection::ApplyVectorBCs(
    const MultiRegions::ExpListSharedPtr &field, const std::size_t dir,
    [[maybe_unused]] const Array<OneD, const NekDouble> &qfield,
    const Array<OneD, const NekDouble> &qFwd,
    [[maybe_unused]] const Array<OneD, const NekDouble> &qBwd,
    Array<OneD, NekDouble> &penaltyflux)
{
    std::size_t nBndRegions = field->GetBndCondExpansions().size();
    std::size_t cnt         = 0;

    for (std::size_t i = 0; i < nBndRegions; ++i)
    {
        if (field->GetBndConditions()[i]->GetBoundaryConditionType() ==
            SpatialDomains::ePeriodic)
        {
            continue;
        }
        std::size_t nBndEdges = field->GetBndCondExpansions()[i]->GetExpSize();

        // Weakly impose boundary conditions by modifying flux values
        for (std::size_t e = 0; e < nBndEdges; ++e)
        {
            std::size_t nBndEdgePts =
                field->GetBndCondExpansions()[i]->GetExp(e)->GetTotPoints();

            std::size_t id1 =
                field->GetBndCondExpansions()[i]->GetPhys_Offset(e);

            std::size_t id2 = field->GetTrace()->GetPhys_Offset(
                field->GetTraceMap()->GetBndCondIDToGlobalTraceID(cnt++));

            // AV boundary conditions
            if (boost::iequals(field->GetBndConditions()[i]->GetUserDefined(),
                               "Wall") ||
                boost::iequals(field->GetBndConditions()[i]->GetUserDefined(),
                               "Symmetry") ||
                boost::iequals(field->GetBndConditions()[i]->GetUserDefined(),
                               "WallViscous") ||
                boost::iequals(field->GetBndConditions()[i]->GetUserDefined(),
                               "WallAdiabatic") ||
                boost::iequals(field->GetBndConditions()[i]->GetUserDefined(),
                               "WallRotational"))
            {
                Vmath::Zero(nBndEdgePts, &penaltyflux[id2], 1);
            }
            // For Dirichlet boundary condition:
            // qflux = q+ - C_11 (u+ -    g_D) (nx, ny)
            else if (field->GetBndConditions()[i]->GetBoundaryConditionType() ==
                     SpatialDomains::eDirichlet)
            {
                Vmath::Vmul(nBndEdgePts, &m_traceNormals[dir][id2], 1,
                            &qFwd[id2], 1, &penaltyflux[id2], 1);
            }
            // For Neumann boundary condition: qflux = g_N
            else if ((field->GetBndConditions()[i])
                         ->GetBoundaryConditionType() ==
                     SpatialDomains::eNeumann)
            {
                Vmath::Vmul(nBndEdgePts, &m_traceNormals[dir][id2], 1,
                            &(field->GetBndCondExpansions()[i]->GetPhys())[id1],
                            1, &penaltyflux[id2], 1);
            }
        }
    }
}

} // namespace PENKNIFE
