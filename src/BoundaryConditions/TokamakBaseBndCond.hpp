#ifndef TOKAMAK_BASEBNDCOND_HPP
#define TOKAMAK_BASEBNDCOND_HPP

#include <string>

#include "../Misc/VariableConverter.hpp"
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <nektar_interface/solver_base/neso_reader.hpp>

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace MR = Nektar::MultiRegions;
namespace SD = Nektar::SpatialDomains;

namespace NESO::Solvers::tokamak
{
class TokamakBaseBndCond;

/// A shared pointer to a boundary condition object
typedef std::shared_ptr<TokamakBaseBndCond> TokamakBaseBndCondSharedPtr;

/// Declaration of the boundary condition factory
typedef LU::NekFactory<
    std::string, TokamakBaseBndCond, const LU::SessionReaderSharedPtr &,
    const Array<OneD, MR::ExpListSharedPtr> &,
    const Array<OneD, MR::DisContFieldSharedPtr> &,
    const Array<OneD, MR::DisContFieldSharedPtr> &,
    Array<OneD, SpatialDomains::BoundaryConditionShPtr>,
    Array<OneD, MultiRegions::ExpListSharedPtr>, const int, const int>
    TokamakBaseBndCondFactory;

/// Declaration of the boundary condition factory singleton
TokamakBaseBndCondFactory &GetTokamakBaseBndCondFactory();

/**
 * @class TokamakBaseBndCond
 * @brief Encapsulates the user-defined boundary conditions for the
 *        tokamak solver.
 */
class TokamakBaseBndCond
{
public:
    virtual ~TokamakBaseBndCond()
    {
    }

    /// Apply the boundary condition
    void Apply(const Array<OneD, const Array<OneD, NekDouble>> &physarray,
               const NekDouble &time = 0);

    /// Apply the Weight of boundary condition
    void ApplyBwdWeight()
    {
        v_ApplyBwdWeight();
    }

    int pe_idx;
    int omega_idx;
    int phi_idx;

    std::vector<int> ni_idx;
    std::vector<int> vi_idx;
    std::vector<int> pi_idx;
    NESOReaderSharedPtr neso_config;

protected:
    /// Session reader
    LU::SessionReaderSharedPtr m_session;

    std::map<int, SpatialDomains::BoundaryConditionShPtr> m_bndConds;
    std::map<int, MultiRegions::ExpListSharedPtr> m_bndExp;
    /// Expansion of boundary adjacent elements
    MultiRegions::ExpListSharedPtr m_bndElmtExp;

    /// Array of fields
    Array<OneD, MR::ExpListSharedPtr> m_fields;
    /// Trace normals
    Array<OneD, Array<OneD, NekDouble>> m_normals;

    // EM fields
    Array<OneD, MR::DisContFieldSharedPtr> B;
    Array<OneD, MR::DisContFieldSharedPtr> E;

    /// Trace EM fields
    Array<OneD, Array<OneD, NekDouble>> B_bnd;
    Array<OneD, Array<OneD, NekDouble>> E_bnd;

    Array<OneD, NekDouble> mag_B;

    /// Space dimension
    int m_spacedim;
    /// Auxiliary object to convert variables
    VariableConverterSharedPtr m_varConv;

    /// Weight for average calculation of diffusion term
    NekDouble m_diffusionAveWeight;

    int m_variable;
    /// Id of the boundary region
    int m_bcRegion;

    int m_nEdgePts;
    int m_nEdgeCoeffs;

    /// Constructor
    TokamakBaseBndCond(const LU::SessionReaderSharedPtr &pSession,
                       const Array<OneD, MR::ExpListSharedPtr> &pFields,
                       const Array<OneD, MR::DisContFieldSharedPtr> &pB,
                       const Array<OneD, MR::DisContFieldSharedPtr> &pE,
                       Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
                       Array<OneD, MultiRegions::ExpListSharedPtr> exp,
                       const int pSpaceDim, const int bcRegion);

    virtual void v_Apply(const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
        const Array<OneD, const Array<OneD, NekDouble>> &physarray,
        const NekDouble &time) = 0;

    virtual void v_ApplyBwdWeight();
};
} // namespace NESO::Solvers::tokamak
#endif