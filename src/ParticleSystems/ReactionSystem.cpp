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
    auto mesh = std::make_shared<ParticleMeshInterface>(this->graph);

    // std::vector<int> reflection_composites;

    // for (auto &[sk, sv] : this->config->get_particle_species_boundary(0))
    // {
    //     if (sv == ParticleBoundaryConditionType::eReflective)
    //     {
    //         reflection_composites.push_back(sk);
    //     }
    // }
    this->boundary = std::make_shared<ReactionsBoundary>(
        Sym<REAL>("TSP"), this->sycl_target, mesh, /*reflection_composites,*/
        this->config, store);
}

} // namespace NESO::Solvers::tokamak