#include "ReactionSystem.hpp"

namespace NESO::Solvers::tokamak
{
std::string ReactionSystem::class_name =
    GetParticleSystemFactory().RegisterCreatorFunction(
        "ReactionSystem", ReactionSystem::create, "Reaction System");

ReactionSystem::ReactionSystem(NESOReaderSharedPtr session,
                               SD::MeshGraphSharedPtr graph)
    : ParticleSystem(session, graph)
{
}

} // namespace NESO::Solvers::tokamak