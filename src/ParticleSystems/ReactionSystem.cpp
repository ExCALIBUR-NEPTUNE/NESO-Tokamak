#include "ReactionSystem.hpp"

namespace NESO::Solvers::tokamak
{
std::string ReactionSystem::className =
    GetParticleSystemFactory().RegisterCreatorFunction(
        "ReactionSystem", ReactionSystem::create, "Reaction System");

ReactionSystem::ReactionSystem(ParticleReaderSharedPtr session,
                               SD::MeshGraphSharedPtr graph)
    : ParticleSystem(session, graph)
{
    auto remove_transform =
        std::make_shared<SimpleRemovalTransformationStrategy>();
    remove_transform_wrapper = std::make_shared<TransformationWrapper>(
        std::vector<std::shared_ptr<MarkingStrategy>>{
            make_marking_strategy<ComparisonMarkerSingle<REAL, LessThanComp>>(
                Sym<REAL>("WEIGHT"), 1e-10)},
        std::dynamic_pointer_cast<TransformationStrategy>(remove_transform));

    auto merge_transform =
        make_transformation_strategy<MergeTransformationStrategy<2>>(
            Sym<REAL>("POSITION"), Sym<REAL>("WEIGHT"), Sym<REAL>("VELOCITY"));
    merge_transform_wrapper = std::make_shared<TransformationWrapper>(
        std::vector<std::shared_ptr<MarkingStrategy>>{
            make_marking_strategy<ComparisonMarkerSingle<REAL, LessThanComp>>(
                Sym<REAL>("WEIGHT"), 0.01)},
        merge_transform);

    auto project_transform = std::make_shared<ProjectTransformation>(
        src_fields, particle_group, cell_id_translation);
    project_transform_wrapper = std::make_shared<TransformationWrapper>(
        std::dynamic_pointer_cast<TransformationStrategy>(project_transform));

    reaction_controller = std::make_shared<ReactionController>(
        std::vector<std::shared_ptr<TransformationWrapper>>{
            project_transform_wrapper, remove_transform_wrapper,
            merge_transform_wrapper},
        std::vector<std::shared_ptr<TransformationWrapper>>{
            project_transform_wrapper, remove_transform_wrapper,
            merge_transform_wrapper},
        Sym<INT>("INTERNAL_STATE"), Sym<REAL>("TOT_REACTION_RATE"));
}

ReactionSystem::~ReactionSystem() {};

} // namespace NESO::Solvers::tokamak