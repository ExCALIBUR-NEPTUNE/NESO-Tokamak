#include "ReactionSystem.hpp"
#include "Reactions.hpp"
#include "Transforms.hpp"

namespace PENKNIFE
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
        Sym<REAL>("TSP"), this->sycl_target, mesh, this->config,
        this->species_map, store);
}

/**
 * @brief Build the reactions and add them to the reaction controller
 */
void ReactionSystem::set_up_reactions()
{
    auto prop_map = get_default_map();

    auto sycl_target = particle_group->sycl_target;
    const int rank   = sycl_target->comm_pair.rank_parent;

    std::uint64_t root_seed = 141351;
    auto rng_kernel = get_uniform_rng_kernel(sycl_target, 40, root_seed);

    for (const auto &v : this->config->get_reactions())
    {
        std::shared_ptr<AbstractReaction> reaction;
        if (std::get<0>(v) == "Ionisation")
        {
            auto target_species = Species(
                std::get<1>(v)[0], this->species_map[std::get<1>(v)[0]].mass,
                this->species_map[std::get<1>(v)[0]].charge,
                this->species_map[std::get<1>(v)[0]].id);

            auto electron_species = Species("ELECTRON");

            if (std::get<2>(v).first == "Fixed")
            {
                auto rate        = std::get<2>(v).second;
                auto energy_rate = std::get<2>(v).second;
                if (this->vdim == 2)
                {
                    reaction = ionise_reaction_fixed<2>(
                        this->sycl_target, target_species, electron_species,
                        rate, energy_rate);
                }
                else if (this->vdim == 3)
                {
                    reaction = ionise_reaction_fixed<3>(
                        this->sycl_target, target_species, electron_species,
                        rate, energy_rate);
                }
            }
            else if (std::get<2>(v).first == "AMJUEL")
            {
                if (this->vdim == 2)
                {
                    reaction = ionise_reaction_amjuel<2>(
                        this->sycl_target, Nnorm, Tnorm, 1. / omega_c,
                        mesh_length * omega_c, target_species,
                        electron_species);
                }
                else if (this->vdim == 3)
                {
                    reaction = ionise_reaction_amjuel<3>(
                        this->sycl_target, Nnorm, Tnorm, 1. / omega_c,
                        mesh_length * omega_c, target_species,
                        electron_species);
                }
            }
        }
        else if (std::get<0>(v) == "Recombination")
        {
            auto electron_species = Species("ELECTRON", 5.5e-4, -1.0);
            auto neutral_species  = Species(
                std::get<1>(v)[0], this->species_map[std::get<1>(v)[0]].mass,
                this->species_map[std::get<1>(v)[0]].charge,
                this->species_map[std::get<1>(v)[0]].id);

            auto marker_species = Species(
                std::get<1>(v)[0], this->species_map[std::get<1>(v)[0]].mass,
                this->species_map[std::get<1>(v)[0]].charge - 1,
                -1 - this->species_map[std::get<1>(v)[0]].id);

            if (std::get<2>(v).first == "Fixed")
            {
                auto rate        = std::get<2>(v).second;
                auto energy_rate = std::get<2>(v).second;

                if (this->vdim == 2)
                {
                    reaction = recombination_reaction_fixed<2>(
                        this->sycl_target, rng_kernel, mesh_length * omega_c,
                        marker_species, electron_species, neutral_species, rate,
                        energy_rate);
                }
                else if (this->vdim == 3)
                {
                    reaction = recombination_reaction_fixed<3>(
                        this->sycl_target, rng_kernel, mesh_length * omega_c,
                        marker_species, electron_species, neutral_species, rate,
                        energy_rate);
                }
            }
            else if (std::get<2>(v).first == "AMJUEL")
            {
                if (this->vdim == 2)
                {
                    reaction = recombination_reaction_amjuel<2>(
                        this->sycl_target, rng_kernel, Nnorm, Tnorm,
                        1. / omega_c, mesh_length * omega_c, marker_species,
                        electron_species, neutral_species);
                }
                else if (this->vdim == 3)
                {
                    reaction = recombination_reaction_amjuel<3>(
                        this->sycl_target, rng_kernel, Nnorm, Tnorm,
                        1. / omega_c, mesh_length * omega_c, marker_species,
                        electron_species, neutral_species);
                }
            }
        }

        else if (std::get<0>(v) == "ChargeExchange")
        {
            auto projectile_species = Species(
                std::get<1>(v)[0], this->species_map[std::get<1>(v)[0]].mass,
                this->species_map[std::get<1>(v)[0]].charge,
                this->species_map[std::get<1>(v)[0]].id);
            auto target_species = Species(
                std::get<1>(v)[1], this->species_map[std::get<1>(v)[1]].mass,
                this->species_map[std::get<1>(v)[1]].charge,
                this->species_map[std::get<1>(v)[1]].id);

            if (std::get<2>(v).first == "Fixed")
            {
                auto rate          = std::get<2>(v).second;
                auto cross_section = std::get<3>(v).second;
                if (this->vdim == 2)
                {
                    reaction = cx_reaction_fixed<2>(
                        this->sycl_target, rng_kernel, mesh_length * omega_c,
                        target_species, projectile_species, rate,
                        cross_section);
                }
                else if (this->vdim == 3)
                {
                    reaction = cx_reaction_fixed<3>(
                        this->sycl_target, rng_kernel, mesh_length * omega_c,
                        target_species, projectile_species, rate,
                        cross_section);
                }
            }
            else if (std::get<2>(v).first == "AMJUEL")
            {
                if (this->vdim == 2)
                {
                    reaction = cx_reaction_amjuel<2>(
                        this->sycl_target, rng_kernel, Nnorm, Tnorm,
                        1. / omega_c, mesh_length * omega_c, target_species,
                        projectile_species);
                }
                else if (this->vdim == 3)
                {
                    reaction = cx_reaction_amjuel<3>(
                        this->sycl_target, rng_kernel, Nnorm, Tnorm,
                        1. / omega_c, mesh_length * omega_c, target_species,
                        projectile_species);
                }
            }
        }
        this->reaction_controller->add_reaction(reaction);
    }
}

/**
 * @brief Finish setup of the class, including projection of source fields
 */
void ReactionSystem::finish_setup(
    std::vector<std::shared_ptr<DisContField>> &src_fields,
    std::vector<Sym<REAL>> &syms, std::vector<int> &components)
{
    this->src_syms       = syms;
    this->src_components = components;

    this->field_project = std::make_shared<FieldProject<DisContField>>(
        src_fields, this->particle_group, this->cell_id_translation);

    auto project_transform = std::make_shared<ProjectTransformation>(
        src_fields, this->src_syms, this->src_components, this->particle_group,
        this->cell_id_translation);
    auto project_transform_wrapper = std::make_shared<TransformationWrapper>(
        std::dynamic_pointer_cast<TransformationStrategy>(project_transform));

    double remove_threshold;
    config->get_session()->LoadParameter("RemoveThreshold", remove_threshold,
                                         1e-12);
    auto remove_transform_wrapper = removal_wrapper(remove_threshold);

    std::vector<std::shared_ptr<TransformationWrapper>> parent_transforms,
        child_transforms;
    child_transforms.push_back(remove_transform_wrapper);

    double merge_threshold, merge_velocity_extents;
    int merge_bins;
    config->get_session()->LoadParameter("MergeThreshold", merge_threshold,
                                         1e-6);
    config->get_session()->LoadParameter("MergeVelocityExtents",
                                         merge_velocity_extents, 3.0);
    config->get_session()->LoadParameter("MergeBins", merge_bins, 5);
    std::shared_ptr<TransformationWrapper> merge_transform_wrapper;

    if (this->ndim == 2 && this->vdim == 2)
    {
        std::array<REAL, 2> velocity_extents;
        std::array<INT, 2> velocity_bins;
        velocity_bins.fill(merge_bins);
        velocity_extents.fill(merge_velocity_extents);
        merge_transform_wrapper =
            merging_wrapper<2>(this->particle_group, merge_threshold,
                               velocity_extents, velocity_bins);
        child_transforms.emplace_back(merge_transform_wrapper);
    }
    else if (this->ndim == 3 && this->vdim == 3)
    {
        std::array<REAL, 3> velocity_extents;
        std::array<INT, 3> velocity_bins;
        velocity_bins.fill(merge_bins);
        velocity_extents.fill(merge_velocity_extents);
        merge_transform_wrapper =
            merging_wrapper<3>(this->particle_group, merge_threshold,
                               velocity_extents, velocity_bins);
        child_transforms.emplace_back(merge_transform_wrapper);
    }
    else if (this->ndim == 2 && this->vdim == 3)
    {
        // TODO 2D3V merging
        // merge_transform_wrapper = legacy_merging_wrapper<3>(merge_threshold);
    }

    std::vector<std::string> src_names{"ELECTRON_SOURCE_DENSITY",
                                       "ELECTRON_SOURCE_ENERGY",
                                       "ELECTRON_SOURCE_MOMENTUM"};

    for (auto &[k, v] : this->config->get_particle_species())
    {
        src_names.push_back(k + "_SOURCE_DENSITY");
        src_names.push_back(k + "_SOURCE_ENERGY");
        src_names.push_back(k + "_SOURCE_MOMENTUM");
    }

    this->zeroer_transform =
        std::make_shared<ParticleDatZeroer<REAL>>(src_names);

    parent_transforms.insert(parent_transforms.end(), child_transforms.begin(),
                             child_transforms.end());

    this->reaction_controller = std::make_shared<ReactionController>(
        parent_transforms, child_transforms);

    set_up_reactions();
    init_output("particle_trajectory.h5part", Sym<REAL>("POSITION"),
                Sym<INT>("INTERNAL_STATE"), Sym<INT>("CELL_ID"),
                Sym<REAL>("VELOCITY"), Sym<REAL>("MAGNETIC_FIELD"),
                Sym<REAL>("ELECTRON_DENSITY"), this->src_syms,
                Sym<REAL>("WEIGHT"), Sym<INT>("ID"),
                Sym<REAL>("TOT_REACTION_RATE"), Sym<REAL>("FLUID_DENSITY"),
                Sym<REAL>("FLUID_TEMPERATURE"));
}

ReactionSystem::ReactionsBoundary::ReactionsBoundary(
    Sym<REAL> time_step_prop_sym, SYCLTargetSharedPtr sycl_target,
    std::shared_ptr<ParticleMeshInterface> mesh, NESOReaderSharedPtr config,
    std::map<std::string, SpeciesInfo> &species, ParameterStoreSharedPtr store)
    : time_step_prop_sym(time_step_prop_sym), sycl_target(sycl_target),
      ndim(mesh->get_ndim()), vdim(3), config(config)
{
    double remove_threshold;
    config->get_session()->LoadParameter("RemoveThreshold", remove_threshold,
                                         1e-12);
    this->remove_wrapper = removal_wrapper(remove_threshold);
    config->read_boundary_regions();

    for (auto &v : this->config->get_surface_reactions())
    {
        std::vector<int> boundary_ids = std::get<2>(v);
        for (int b_id : boundary_ids)
        {
            if (!this->reaction_controllers[b_id])
            {
                reaction_controllers[b_id] =
                    std::make_shared<ReactionController>(
                        std::vector<std::shared_ptr<TransformationWrapper>>{
                            remove_wrapper},
                        std::vector<std::shared_ptr<TransformationWrapper>>{
                            remove_wrapper});
            }

            if (std::get<0>(v) == "Specular")
            {
                for (const auto &s : std::get<1>(v))
                {
                    auto reflected_species = Species(
                        s, species[s].mass, species[s].charge, species[s].id);
                    auto rate = std::get<3>(v).second;
                    if (this->ndim == 2)
                    {
                        if (this->vdim == 2)
                        {
                            auto reaction = specular_reflection<2, 2>(
                                this->sycl_target, reflected_species, rate);

                            this->reaction_controllers[b_id]->add_reaction(
                                reaction);
                        }
                        else if (this->vdim == 3)
                        {
                            auto reaction = specular_reflection<2, 3>(
                                this->sycl_target, reflected_species, rate);

                            this->reaction_controllers[b_id]->add_reaction(
                                reaction);
                        }
                    }
                    else if (this->ndim == 3)
                    {
                        auto reaction = specular_reflection<3, 3>(
                            this->sycl_target, reflected_species, rate);

                        this->reaction_controllers[b_id]->add_reaction(
                            reaction);
                    }
                }
            }
            else if (std::get<0>(v) == "Thermal")
            {
                for (const auto &s : std::get<1>(v))
                {
                    double T             = 1.0;
                    auto thermal_species = Species(
                        s, species[s].mass, species[s].charge, species[s].id);
                    auto rate = std::get<3>(v).second;

                    if (this->ndim == 2)
                    {
                        if (this->vdim == 2)
                        {
                            auto reaction = thermal_reflection<2, 2>(
                                this->sycl_target, thermal_species, rate, T);

                            this->reaction_controllers[b_id]->add_reaction(
                                reaction);
                        }
                        else if (this->vdim == 3)
                        {
                            auto reaction = thermal_reflection<2, 3>(
                                this->sycl_target, thermal_species, rate, T);

                            this->reaction_controllers[b_id]->add_reaction(
                                reaction);
                        }
                    }
                    else if (this->ndim == 3)
                    {
                        auto reaction = thermal_reflection<3, 3>(
                            this->sycl_target, thermal_species, rate, T);
                        this->reaction_controllers[b_id]->add_reaction(
                            reaction);
                    }
                }
            }
            else if (std::get<0>(v) == "Absorption")
            {
                for (const auto &s : std::get<1>(v))
                {
                    auto absorbed_species = Species(
                        s, species[s].mass, species[s].charge, species[s].id);
                    auto rate = std::get<3>(v).second;

                    if (this->vdim == 2)
                    {
                        auto reaction = surface_absorption<2>(
                            this->sycl_target, absorbed_species, rate);

                        this->reaction_controllers[b_id]->add_reaction(
                            reaction);
                    }
                    else if (this->vdim == 3)
                    {
                        auto reaction = surface_absorption<3>(
                            this->sycl_target, absorbed_species, rate);

                        this->reaction_controllers[b_id]->add_reaction(
                            reaction);
                    }
                }
            }
        }
    }

    this->composite_intersection =
        std::make_shared<CompositeInteraction::CompositeIntersection>(
            this->sycl_target, mesh, config->get_boundary_regions());

    this->reset_distance =
        store->get<REAL>("ReactionsBoundary/reset_distance", 1.0e-4);

    this->boundary_truncation =
        std::make_shared<BoundaryTruncation>(this->ndim, this->reset_distance);
}

} // namespace PENKNIFE