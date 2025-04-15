#include "MultiFieldUpwind.hpp"

namespace NESO::Solvers::tokamak
{
std::string MultiFieldUpwindSolver::solverName =
    SU::GetRiemannSolverFactory().RegisterCreatorFunction(
        "MultiFieldUpwind", MultiFieldUpwindSolver::create, "Upwind solver");

/**
 * @class UpwindSolver
 *
 * @brief Upwind scheme Riemann solver.
 *
 * The upwind solver determines the flux based upon an advection
 * velocity \f$\mathbf{V}\f$ and trace normal \f$\mathbf{n}\f$. In
 * particular, the flux for each component of the velocity field is
 * deterined by:
 *
 * \f[ \mathbf{f}(u^+,u^-) = \begin{cases} \mathbf{V}u^+, &
 * \mathbf{V}\cdot\mathbf{n} \geq 0,\\ \mathbf{V}u^-, &
 * \mathbf{V}\cdot\mathbf{n} < 0.\end{cases} \f]
 *
 * Here the superscript + and - denotes forwards and backwards spaces
 * respectively.
 */

/**
 * @brief Default constructor.
 */
MultiFieldUpwindSolver::MultiFieldUpwindSolver(
    const LU::SessionReaderSharedPtr &pSession)
    : RiemannSolver(pSession)
{
}

/**
 * @brief Implementation of the upwind solver.
 *
 * The upwind solver assumes that a scalar field Vn is defined, which
 * corresponds with the dot product \f$\mathbf{V}\cdot\mathbf{n}\f$,
 * where \f$\mathbf{V}\f$ is the advection velocity and \f$\mathbf{n}\f$
 * defines the normal of a vertex, edge or face at each quadrature point
 * of the trace space.
 *
 * @param Fwd   Forwards trace space.
 * @param Bwd   Backwards trace space.
 * @param flux  Resulting flux.
 */
void MultiFieldUpwindSolver::v_Solve(
    [[maybe_unused]] const int nDim,
    const Array<OneD, const Array<OneD, NekDouble>> &Fwd,
    const Array<OneD, const Array<OneD, NekDouble>> &Bwd,
    Array<OneD, Array<OneD, NekDouble>> &flux)
{
    ASSERTL1(CheckVectors("Vn"), "Vn not defined.");
    const Array<OneD, Array<OneD, NekDouble>> &traceVel = m_vectors["Vn"]();

    for (int i = 0; i < Fwd.size(); ++i)
    {
        for (int j = 0; j < traceVel.size(); ++j)
        {
            NekDouble tmp = traceVel[i][j] >= 0 ? Fwd[i][j] : Bwd[i][j];
            flux[i][j]    = traceVel[i][j] * tmp;
        }
    }
}
} // namespace NESO::Solvers::tokamak
