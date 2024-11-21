#include "ReactionsCoupling.hpp"
using namespace std;

namespace NESO::Solvers::tokamak
{
std::string ReactionsCoupling::className =
    SolverUtils::GetForcingFactory().RegisterCreatorFunction(
        "ReactionsCoupling", ReactionsCoupling::create,
        "Coupling to Reactions");

ReactionsCoupling::ReactionsCoupling(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<SolverUtils::EquationSystem> &pEquation)
    : Forcing(pSession, pEquation), field_to_index(pSession->GetVariables())
{
}

void ReactionsCoupling::v_InitObject(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const unsigned int &pNumForcingFields, const TiXmlElement *pForce)
{
    m_NumVariable = pNumForcingFields;
}

void ReactionsCoupling::v_Apply(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray,
    [[maybe_unused]] const NekDouble &time)
{
    int npoints = pFields[0]->GetTotPoints();

    int ni_field_idx = this->field_to_index.get_idx("n");
    int ni_src_idx   = this->field_to_index.get_idx("n_src");
    Vmath::Vadd(npoints, inarray[ni_src_idx], 1, outarray[ni_field_idx], 1,
                outarray[ni_field_idx], 1);

    int p_field_idx = this->field_to_index.get_idx("p");
    int E_src_idx   = this->field_to_index.get_idx("E_src");
    Vmath::Vadd(npoints, inarray[E_src_idx], 1, outarray[p_field_idx], 1,
                outarray[p_field_idx], 1);
}
} // namespace NESO::Solvers::tokamak