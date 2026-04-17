#ifndef TRANSFORMS_HPP
#define TRANSFORMS_HPP

#include <reactions/reactions.hpp>

using namespace VANTAGE::Reactions;

namespace PENKNIFE
{
inline auto removal_wrapper(REAL threshold)
{
    auto remove_transform_wrapper = std::make_shared<TransformationWrapper>(
        std::vector<std::shared_ptr<MarkingStrategy>>{
            make_direct_marking_strategy(
                "removal_marker", [=](auto w) { return w[0] < threshold; },
                Access::read(Sym<REAL>("WEIGHT")))},
        make_transformation_strategy<SimpleRemovalTransformationStrategy>());

    return remove_transform_wrapper;
};

template <size_t ndim> inline auto legacy_merging_wrapper(REAL threshold)
{
    auto merge_transform =
        make_transformation_strategy<MergeTransformationStrategy<ndim>>();

    auto merge_transform_wrapper = std::make_shared<TransformationWrapper>(
        std::vector<std::shared_ptr<MarkingStrategy>>{
            make_direct_marking_strategy(
                "merge_marker", [=](auto w) { return w[0] < threshold; },
                Access::read(Sym<REAL>("WEIGHT")))},
        merge_transform);

    return merge_transform_wrapper;
}
template <size_t vdim>
inline auto merging_wrapper(
    ParticleGroupSharedPtr particle_group, REAL threshold,
    const std::array<REAL, vdim> &global_velocity_extents,
    const std::array<INT, vdim> &n_vel_cells)
{
    // TransformationStrategy that merges N particles in a cell such that only 2
    // particles are now present with the averaged positions of all N particles,
    // the sum of all the weights of the N particles (split between the 2
    // particles) and the appropriate velocities to ensure momentum conservation

    REAL n_merging_cells = 1;
    for (auto n_cell : n_vel_cells)
    {
        n_merging_cells *= n_cell + 2;
    }
    auto merge_transform =
        make_vranic_merging_strategy<vdim>(particle_group, n_merging_cells);

    // TransformationWrapper with a marking strategy that specifies that only
    // particles with "WEIGHT" below threshold should be considered for merging.
    auto merge_transform_wrapper = std::make_shared<TransformationWrapper>(
        std::vector<std::shared_ptr<MarkingStrategy>>{
            make_direct_marking_strategy(
                "merge_marker", [=](auto w) { return w[0] < threshold; },
                Access::read(Sym<REAL>("WEIGHT")))},
        make_transformation_strategy<CompositeTransform>(
            std::vector<std::shared_ptr<TransformationStrategy>>{
                uniform_velocity_bin_transform(
                    global_velocity_extents, n_vel_cells,
                    Sym<INT>("REACTIONS_GROUPING_INDEX"),
                    Sym<REAL>("VELOCITY")),
                merge_transform}));

    return merge_transform_wrapper;
}
} // namespace PENKNIFE
#endif