#ifndef TOKAMAK_BNDCOND_HPP
#define TOKAMAK_BNDCOND_HPP

#include <string>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include "../Misc/VariableConverter.hpp"

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;

namespace NESO::Solvers::tokamak
{
class TokamakBndCond;

/// A shared pointer to a boundary condition object
typedef std::shared_ptr<TokamakBndCond> TokamakBndCondSharedPtr;

/// Declaration of the boundary condition factory
typedef LU::NekFactory<
    std::string, TokamakBndCond, const LU::SessionReaderSharedPtr &,
    const Array<OneD, MR::ExpListSharedPtr> &,
    const Array<OneD, Array<OneD, NekDouble>> &, const int, const int>
    TokamakBndCondFactory;

/// Declaration of the boundary condition factory singleton
TokamakBndCondFactory &GetTokamakBndCondFactory();

/**
 * @class TokamakBndCond
 * @brief Encapsulates the user-defined boundary conditions for the
 *        tokamak solver.
 */
class TokamakBndCond
{
public:
    virtual ~TokamakBndCond()
    {
    }

    /// Apply the boundary condition
    void Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
               const NekDouble &time = 0);

    /// Apply the Weight of boundary condition
    void ApplyBwdWeight()
    {
        v_ApplyBwdWeight();
    }

protected:
    /// Session reader
    LU::SessionReaderSharedPtr m_session;
    /// Array of fields
    Array<OneD, MR::ExpListSharedPtr> m_fields;
    /// Trace normals
    Array<OneD, Array<OneD, NekDouble>> m_normals;
    /// Unit magnetic field
    Array<OneD, Array<OneD, NekDouble>> m_b;
    /// Space dimension
    int m_spacedim;
    /// Auxiliary object to convert variables
    VariableConverterSharedPtr m_varConv;

    /// Weight for average calculation of diffusion term
    NekDouble m_diffusionAveWeight;

    /// Id of the boundary region
    int m_bcRegion;

    int m_nEdgePts;
    int m_nEdgeCoeffs;

    /// Expansion of boundary adjacent elements
    MultiRegions::ExpListSharedPtr m_bndElmtExp;
    /// Expansion of boundary
    Array<OneD, MultiRegions::ExpListSharedPtr> m_bndExp;

    /// Constructor
    TokamakBndCond(const LU::SessionReaderSharedPtr &pSession,
                   const Array<OneD, MR::ExpListSharedPtr> &pFields,
                   const Array<OneD, Array<OneD, NekDouble>> &pMagneticField,
                   const int pSpaceDim, const int bcRegion);

    virtual void v_Apply(Array<OneD, Array<OneD, NekDouble>> &physarray,
                         const NekDouble &time) = 0;

    virtual void v_ApplyBwdWeight();
};
} // namespace NESO::Solvers::tokamak
#endif