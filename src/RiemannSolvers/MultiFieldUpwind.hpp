#ifndef MULTIFIELD_UPWINDSOLVER
#define MULTIFIELD_UPWINDSOLVER

#include <SolverUtils/RiemannSolvers/RiemannSolver.h>

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace SU = Nektar::SolverUtils;
namespace NESO::Solvers::tokamak
{
class MultiFieldUpwindSolver : public SU::RiemannSolver
{
public:
    static SU::RiemannSolverSharedPtr create(
        const LU::SessionReaderSharedPtr &pSession)
    {
        return SU::RiemannSolverSharedPtr(new MultiFieldUpwindSolver(pSession));
    }

    static std::string solverName;

protected:
    MultiFieldUpwindSolver(
        const LU::SessionReaderSharedPtr &pSession);

    void v_Solve(const int nDim,
                 const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
                 const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
                 Array<OneD, Array<OneD, NekDouble>> &flux) final;
};
} // namespace NESO::Solvers::tokamak

#endif
