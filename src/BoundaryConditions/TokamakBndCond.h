#ifndef TOKAMAK_BNDCOND
#define TOKAMAK_BNDCOND

#include <string>

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>

namespace NESO::Solvers::tokamak
{
class TokamakBndCond;

/// A shared pointer to a boundary condition object
typedef std::shared_ptr<TokamakBndCond> TokamakBndCondSharedPtr;

/// Declaration of the boundary condition factory
typedef LibUtilities::NekFactory<
    std::string, TokamakBndCond, const LibUtilities::SessionReaderSharedPtr &,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &,
    const Array<OneD, Array<OneD, NekDouble>> &,
    const Array<OneD, Array<OneD, NekDouble>> &, const int, const int,
    const int>
    TokamakBndCondFactory;

/// Declaration of the boundary condition factory singleton
TokamakBndCondFactory &GetDiffBndCondFactory();

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
    void Apply(Array<OneD, Array<OneD, NekDouble>> &magnetic,
               Array<OneD, Array<OneD, NekDouble>> &physarray,
               const NekDouble &time = 0);

    /// Apply the Weight of boundary condition
    void ApplyBwdWeight()
    {
        v_ApplyBwdWeight();
    }

protected:
    /// Session reader
    LibUtilities::SessionReaderSharedPtr m_session;
    /// Array of fields
    Array<OneD, MultiRegions::ExpListSharedPtr> m_fields;
    /// Trace normals
    Array<OneD, Array<OneD, NekDouble>> m_normals;
    /// Oblique Field
    Array<OneD, Array<OneD, NekDouble>> m_obliqueField;
    /// Space dimension
    int m_spacedim;
    /// Weight for average calculation of diffusion term
    NekDouble m_diffusionAveWeight;

    /// Id of the boundary region
    int m_bcRegion;
    /// Offset
    int m_offset;

    /// Constructor
    TokamakBndCond(const LibUtilities::SessionReaderSharedPtr &pSession,
                   const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
                   const Array<OneD, Array<OneD, NekDouble>> &pTraceNormals,
                   const Array<OneD, Array<OneD, NekDouble>> &pObliqueField,
                   const int pSpaceDim, const int bcRegion, const int cnt);

    virtual void v_Apply(Array<OneD, Array<OneD, NekDouble>> &FwdOblique,
                         Array<OneD, Array<OneD, NekDouble>> &physarray,
                         const NekDouble &time) = 0;

    virtual void v_ApplyBwdWeight();
};
} // namespace NESO::Solvers::tokamak
#endif