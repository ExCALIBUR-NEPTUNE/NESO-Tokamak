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
void ReactionSystem::set_up_boundaries()
{
    auto store = std::make_shared<ParameterStore>();
    store->set<REAL>("ReactionsBoundary/reset_distance", 1.0e-6);
    auto mesh      = std::make_shared<ParticleMeshInterface>(this->graph);
    this->boundary = std::make_shared<ReactionsBoundary>(
        Sym<REAL>("TSP"), this->sycl_target, mesh, this->config, store);
}

} // namespace NESO::Solvers::tokamak