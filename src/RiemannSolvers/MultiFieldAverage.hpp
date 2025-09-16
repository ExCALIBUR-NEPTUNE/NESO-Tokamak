#ifndef MULTIFIELD_AVERAGESOLVER
#define MULTIFIELD_AVERAGESOLVER

#include "TokamakSolver.hpp"

using namespace Nektar;
namespace LU = Nektar::LibUtilities;
namespace SU = Nektar::SolverUtils;
namespace NESO::Solvers::tokamak
{
class MultiFieldAverageSolver : public TokamakSolver
{
public:
    static SU::RiemannSolverSharedPtr create(
        const LU::SessionReaderSharedPtr &pSession)
    {
        return SU::RiemannSolverSharedPtr(
            new MultiFieldAverageSolver(pSession));
    }

    static std::string solverName;

protected:
    MultiFieldAverageSolver(const LU::SessionReaderSharedPtr &pSession);

    void v_ArraySolve(const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
                      const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
                      Array<OneD, Array<OneD, NekDouble>> &flux) final;
};
} // namespace NESO::Solvers::tokamak

#endif
