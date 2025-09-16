#ifndef VORTICITY_UPWINDSOLVER
#define VORTICITY_UPWINDSOLVER

#include "TokamakSolver.hpp"

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace SU = Nektar::SolverUtils;
namespace NESO::Solvers::tokamak
{
class VorticityUpwindSolver : public TokamakSolver
{
public:
    static SU::RiemannSolverSharedPtr create(
        const LU::SessionReaderSharedPtr &pSession)
    {
        return SU::RiemannSolverSharedPtr(new VorticityUpwindSolver(pSession));
    }

    static std::string solverName;

protected:
    VorticityUpwindSolver(const LU::SessionReaderSharedPtr &pSession);

    void v_ArraySolve(const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
                      const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
                      Array<OneD, Array<OneD, NekDouble>> &flux) final;
};
} // namespace NESO::Solvers::tokamak

#endif
