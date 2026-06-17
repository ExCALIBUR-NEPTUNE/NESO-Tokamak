#ifndef DIFFUSIONWEAKDGET
#define DIFFUSIONWEAKDGET

#include <SolverUtils/Diffusion/Diffusion.h>

using namespace Nektar::SolverUtils;
using namespace Nektar;

namespace PENKNIFE
{

class DiffusionLDGET : public Diffusion
{
public:
    static DiffusionSharedPtr create([[maybe_unused]] std::string diffType)
    {
        return DiffusionSharedPtr(new DiffusionLDGET());
    }

    static std::string type;

protected:
    DiffusionLDGET();

    void v_InitObject(
        LibUtilities::SessionReaderSharedPtr pSession,
        Array<OneD, MultiRegions::ExpListSharedPtr> pFields) override;

    void v_Diffuse(const std::size_t nConvective,
                   const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
                   const Array<OneD, Array<OneD, NekDouble>> &inarray,
                   Array<OneD, Array<OneD, NekDouble>> &outarray,
                   const Array<OneD, Array<OneD, NekDouble>> &pFwd,
                   const Array<OneD, Array<OneD, NekDouble>> &pBwd) override;

    void v_DiffuseCoeffs(
        const std::size_t nConvective,
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        Array<OneD, Array<OneD, NekDouble>> &outarray,
        const Array<OneD, Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, Array<OneD, NekDouble>> &pBwd) override;

    void v_DiffuseCalcDerivative(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        TensorOfArray3D<NekDouble> &qfields,
        const Array<OneD, Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, Array<OneD, NekDouble>> &pBwd) override;

    void v_DiffuseVolumeFlux(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        TensorOfArray3D<NekDouble> &qfields,
        TensorOfArray3D<NekDouble> &VolumeFlux,
        Array<OneD, int> &nonZeroIndex) override;

    void v_DiffuseTraceFlux(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &inarray,
        TensorOfArray3D<NekDouble> &qfields,
        TensorOfArray3D<NekDouble> &VolumeFlux,
        Array<OneD, Array<OneD, NekDouble>> &TraceFlux,
        const Array<OneD, Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, Array<OneD, NekDouble>> &pBwd,
        Array<OneD, int> &nonZeroIndex) override;

private:
    std::string m_shockCaptureType;

    /// Coefficient of penalty term
    NekDouble m_C11;

    Array<OneD, Array<OneD, NekDouble>> m_traceNormals;
    LibUtilities::SessionReaderSharedPtr m_session;

    void NumFluxforScalar(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &ufield,
        TensorOfArray3D<NekDouble> &uflux,
        const Array<OneD, Array<OneD, NekDouble>> &pFwd,
        const Array<OneD, Array<OneD, NekDouble>> &pBwd);

    void ApplyScalarBCs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const std::size_t var, const Array<OneD, const NekDouble> &ufield,
        const Array<OneD, const NekDouble> &Fwd,
        const Array<OneD, const NekDouble> &Bwd,
        Array<OneD, NekDouble> &penaltyflux);

    void NumFluxforVector(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const Array<OneD, Array<OneD, NekDouble>> &ufield,
        TensorOfArray3D<NekDouble> &qfield,
        Array<OneD, Array<OneD, NekDouble>> &qflux,
        const Array<OneD, Array<OneD, NekDouble>> &uFwd,
        const Array<OneD, Array<OneD, NekDouble>> &uBwd);

    void ApplyVectorBCs(
        const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
        const std::size_t var, const std::size_t dir,
        const Array<OneD, const NekDouble> &qfield,
        const Array<OneD, const NekDouble> &qFwd,
        const Array<OneD, const NekDouble> &qBwd,
        Array<OneD, NekDouble> &penaltyflux);
};

} // namespace PENKNIFE

#endif
