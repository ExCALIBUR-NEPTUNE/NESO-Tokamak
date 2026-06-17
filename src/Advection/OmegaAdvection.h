///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLDG.h
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

#ifndef NEKTAR_SOLVERUTILS_DIFFUSIONWEAKDG
#define NEKTAR_SOLVERUTILS_DIFFUSIONWEAKDG

#include <SolverUtils/Diffusion/Diffusion.h>

using namespace Nektar::SolverUtils;
using namespace Nektar;
namespace PENKNIFE
{
class OmegaAdvection
{
public:
    int omega_idx;
    OmegaAdvection();

    void InitObject(LibUtilities::SessionReaderSharedPtr pSession,
                    MultiRegions::ExpListSharedPtr pFields);

    void Advect(const MultiRegions::ExpListSharedPtr &fields,
                const Array<OneD, Array<OneD, NekDouble>> &inarray,
                Array<OneD, Array<OneD, NekDouble>> &outarray,
                const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    template <typename FuncPointerT, typename ObjectPointerT>
    void SetFluxVector(FuncPointerT func, ObjectPointerT obj)
    {
        m_flux_vector = std::bind(func, obj, std::placeholders::_1,
                                  std::placeholders::_2, std::placeholders::_3);
    }

    void SetFluxVector(
        std::function<void(const Array<OneD, Array<OneD, NekDouble>> &,
                           Array<OneD, Array<OneD, NekDouble>> &,
                           Array<OneD, NekDouble> &)>
            fluxVector)
    {
        m_flux_vector = fluxVector;
    }

private:
    void AdvectCoeffs(const MultiRegions::ExpListSharedPtr &fields,
                      const Array<OneD, Array<OneD, NekDouble>> &inarray,
                      Array<OneD, NekDouble> &outarray,
                      const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                      const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    void AdvectCalcDerivative(
        const MultiRegions::ExpListSharedPtr &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &qfields,
        const Array<OneD, Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    void AdvectFlux(const MultiRegions::ExpListSharedPtr &fields,
                    const Array<OneD, Array<OneD, NekDouble>> &inarray,
                    Array<OneD, Array<OneD, NekDouble>> &VolumeFlux,
                    Array<OneD, NekDouble> &TraceFlux);

    // void AdvectVolumeFlux(const MultiRegions::ExpListSharedPtr &fields,
    //                       const Array<OneD, Array<OneD, NekDouble>> &inarray,
    //                       Array<OneD, Array<OneD, NekDouble>> &VolumeFlux);

    void AdvectTraceFlux(const MultiRegions::ExpListSharedPtr &fields,
                         const Array<OneD, Array<OneD, NekDouble>> &inarray,
                         Array<OneD, Array<OneD, NekDouble>> &VolumeFlux,
                         Array<OneD, NekDouble> &TraceFlux,
                         const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                         const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    std::string m_shockCaptureType;

    /// Coefficient of penalty term
    NekDouble m_C11;

    Array<OneD, Array<OneD, NekDouble>> m_traceNormals;
    LibUtilities::SessionReaderSharedPtr m_session;

    std::function<void(const Array<OneD, Array<OneD, NekDouble>> &,
                       Array<OneD, Array<OneD, NekDouble>> &,
                       Array<OneD, NekDouble> &)>
        m_flux_vector;

    void NumFluxforScalar(const MultiRegions::ExpListSharedPtr &fields,
                          const Array<OneD, NekDouble> &ufield,
                          Array<OneD, Array<OneD, NekDouble>> &uflux,
                          const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                          const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    void ApplyScalarBCs(const MultiRegions::ExpListSharedPtr &fields,
                        const Array<OneD, const NekDouble> &ufield,
                        const Array<OneD, const NekDouble> &Fwd,
                        const Array<OneD, const NekDouble> &Bwd,
                        Array<OneD, NekDouble> &penaltyflux);

    void NumFluxforVector(const MultiRegions::ExpListSharedPtr &field,
                          const Array<OneD, NekDouble> &ufield,
                          Array<OneD, Array<OneD, NekDouble>> &qfield,
                          Array<OneD, NekDouble> &qflux);

    void ApplyVectorBCs(const MultiRegions::ExpListSharedPtr &field,
                        const std::size_t dir,
                        const Array<OneD, const NekDouble> &qfield,
                        const Array<OneD, const NekDouble> &qFwd,
                        const Array<OneD, const NekDouble> &qBwd,
                        Array<OneD, NekDouble> &penaltyflux);
};
} // namespace PENKNIFE

#endif
