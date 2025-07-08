#ifndef REACTIONSYSTEM_HPP
#define REACTIONSYSTEM_HPP

#include "AMJUEL.hpp"
#include "ParticleSystem.hpp"
#include <neso_rng_toolkit.hpp>
#include <reactions.hpp>

using namespace Reactions;

namespace NESO::Solvers::tokamak
{

class ReactionSystem : public ParticleSystem
{

public:
    static std::string class_name;
    static ParticleSystemSharedPtr create(const NESOReaderSharedPtr session,
                                          const SD::MeshGraphSharedPtr graph)
    {
        ParticleSystemSharedPtr p =
            MemoryManager<ReactionSystem>::AllocateSharedPtr(session, graph);
        return p;
    }

    ReactionSystem(NESOReaderSharedPtr session, SD::MeshGraphSharedPtr graph);

    ~ReactionSystem() override = default;

    std::shared_ptr<TransformationWrapper> zeroer_transform_wrapper;

    inline void integrate(const double time_end, const double dt) override
    {
        this->zeroer_transform_wrapper->transform(this->particle_group);
        ParticleSystem::integrate(time_end, dt);
        this->field_project->project(this->particle_group, this->src_syms,
                                     this->src_components);
    }

    inline void set_up_reactions()
    {
        auto prop_map = get_default_map();

        auto sycl_target = particle_group->sycl_target;
        const int rank     = sycl_target->comm_pair.rank_parent;

        std::uint64_t root_seed = 141351;
        std::uint64_t seed = NESO::RNGToolkit::create_seeds(
            sycl_target->comm_pair.size_parent, rank, root_seed);

        // Create a Normal distribution with mean 4.0 and standard
        // deviation 2.0.
        auto rng_normal = NESO::RNGToolkit::create_rng<REAL>(
            NESO::RNGToolkit::Distribution::Uniform<REAL>{
                NESO::RNGToolkit::Distribution::next_value(0.0), 1.0},
            seed, sycl_target->device, sycl_target->device_index);

        // Create an interface between NESO-RNG-Toolkit and NESO-Particles
        // KernelRNG
        auto rng_interface =
            make_rng_generation_function<GenericDeviceRNGGenerationFunction,
                                         REAL>(
                [=](REAL *d_ptr, const std::size_t num_samples) -> int
                { return rng_normal->get_samples(d_ptr, num_samples); });

        auto rng_kernel =
            host_atomic_block_kernel_rng<REAL>(rng_interface, 4 * 10);

        for (const auto &[k, v] : this->config->get_reactions())
        {
            std::shared_ptr<AbstractReaction> reaction;
            if (std::get<0>(v) == "Ionisation")
            {
                auto target_species =
                    Species(this->species_map[std::get<1>(v)[0]].name,
                            this->species_map[std::get<1>(v)[0]].mass,
                            this->species_map[std::get<1>(v)[0]].charge,
                            std::get<1>(v)[0]);

                auto electron_species = Species("ELECTRON");

                if (std::get<2>(v).first == "Fixed")
                {
                    auto ionise_data        = FixedRateData(1.0);
                    auto ionise_energy_data = FixedRateData(1.0);
                    if (this->ndim == 2)
                    {
                        reaction = std::make_shared<ElectronImpactIonisation<
                            decltype(ionise_data), decltype(ionise_energy_data),
                            2>>(this->particle_group->sycl_target, ionise_data,
                                ionise_energy_data, target_species,
                                electron_species);
                    }
                    else if (this->ndim == 3)
                    {
                        reaction = std::make_shared<ElectronImpactIonisation<
                            decltype(ionise_data), decltype(ionise_energy_data),
                            3>>(this->particle_group->sycl_target, ionise_data,
                                ionise_energy_data, target_species,
                                electron_species);
                    }
                }
                else if (std::get<2>(v).first == "AMJUEL")
                {
                    auto ionise_data        = AMJUEL::ionise_data;
                    auto ionise_energy_data = AMJUEL::ionise_energy_data;

                    if (this->ndim == 2)
                    {
                        reaction = std::make_shared<ElectronImpactIonisation<
                            decltype(ionise_data), decltype(ionise_energy_data),
                            2>>(this->particle_group->sycl_target, ionise_data,
                                ionise_energy_data, target_species,
                                electron_species);
                    }
                    else if (this->ndim == 3)
                    {
                        reaction = std::make_shared<ElectronImpactIonisation<
                            decltype(ionise_data), decltype(ionise_energy_data),
                            3>>(this->particle_group->sycl_target, ionise_data,
                                ionise_energy_data, target_species,
                                electron_species);
                    }
                }
            }
            else if (std::get<0>(v) == "Recombination")
            {
                auto electron_species = Species("ELECTRON", 5.5e-4, -1.0);
                auto neutral_species =
                    Species(this->species_map[std::get<1>(v)[0]].name,
                            this->species_map[std::get<1>(v)[0]].mass,
                            this->species_map[std::get<1>(v)[0]].charge,
                            std::get<1>(v)[0]);

                auto marker_species =
                    Species(this->species_map[std::get<1>(v)[0]].name,
                            this->species_map[std::get<1>(v)[0]].mass,
                            this->species_map[std::get<1>(v)[0]].charge - 1,
                            -1 - std::get<1>(v)[0]);

                if (std::get<2>(v).first == "Fixed")
                {
                    auto recomb_data = FixedRateData(1.0);
                    auto data1       = FixedRateData(1.0);
                    auto data2       = FixedRateData(1.0);
                    auto data_calculator =
                        DataCalculator<FixedRateData, FixedRateData,
                                       FixedRateData>(data1, data1, data2);
                    if (this->ndim == 2)
                    {
                        auto recomb_reaction_kernel = RecombReactionKernels<2>(
                            marker_species, electron_species,
                            norm::potential_energy);
                        reaction = std::make_shared<
                            LinearReactionBase<1, decltype(recomb_data),
                                               decltype(recomb_reaction_kernel),
                                               decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            marker_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(neutral_species.get_id())},
                            recomb_data, recomb_reaction_kernel,
                            data_calculator);
                    }
                    else if (this->ndim == 3)
                    {
                        auto recomb_reaction_kernel = RecombReactionKernels<3>(
                            marker_species, electron_species,
                            norm::potential_energy);
                        reaction = std::make_shared<
                            LinearReactionBase<1, decltype(recomb_data),
                                               decltype(recomb_reaction_kernel),
                                               decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            marker_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(neutral_species.get_id())},
                            recomb_data, recomb_reaction_kernel,
                            data_calculator);
                    }
                }
                else if (std::get<2>(v).first == "AMJUEL")
                {
                    auto recomb_data        = AMJUEL::recomb_data;
                    auto recomb_energy_data = AMJUEL::recomb_energy_data;

                    auto constant_rate_cross_section =
                        ConstantRateCrossSection(1.0);
                    auto recomb_data_calc_sampler = FilteredMaxwellianSampler<
                        2, decltype(constant_rate_cross_section)>(
                        (norm::temp_SI * norm::kB) /
                            (marker_species.get_mass() * norm::mass_amu_SI *
                             norm::vel * norm::vel),
                        constant_rate_cross_section, rng_kernel);
                    auto data_calculator =
                        DataCalculator<decltype(recomb_energy_data),
                                       decltype(recomb_data_calc_sampler)>(
                            recomb_energy_data, recomb_data_calc_sampler);

                    if (this->ndim == 2)
                    {
                        auto recomb_reaction_kernel = RecombReactionKernels<2>(
                            marker_species, electron_species,
                            norm::potential_energy);
                        reaction = std::make_shared<
                            LinearReactionBase<1, decltype(recomb_data),
                                               decltype(recomb_reaction_kernel),
                                               decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            marker_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(neutral_species.get_id())},
                            recomb_data, recomb_reaction_kernel,
                            data_calculator);
                    }
                    else if (this->ndim == 3)
                    {
                        auto recomb_reaction_kernel = RecombReactionKernels<3>(
                            marker_species, electron_species,
                            norm::potential_energy);
                        reaction = std::make_shared<
                            LinearReactionBase<1, decltype(recomb_data),
                                               decltype(recomb_reaction_kernel),
                                               decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            marker_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(neutral_species.get_id())},
                            recomb_data, recomb_reaction_kernel,
                            data_calculator);
                    }
                }
            }

            else if (std::get<0>(v) == "ChargeExchange")
            {
                auto projectile_species =
                    Species(this->species_map[std::get<1>(v)[0]].name,
                            this->species_map[std::get<1>(v)[0]].mass,
                            this->species_map[std::get<1>(v)[0]].charge,
                            std::get<1>(v)[0]);
                auto target_species =
                    Species(this->species_map[std::get<1>(v)[1]].name,
                            this->species_map[std::get<1>(v)[1]].mass,
                            this->species_map[std::get<1>(v)[1]].charge,
                            std::get<1>(v)[1]);
                auto parent_mass = target_species.get_mass();
                auto child_mass  = projectile_species.get_mass();
                auto reduced_mass =
                    (parent_mass * child_mass) / (parent_mass + child_mass);

                if (std::get<2>(v).first == "Fixed")
                {
                    auto rate_data              = FixedRateData(1.0);
                    auto constant_cross_section = ConstantRateCrossSection(1.0);

                    auto data_calc_sampler = FilteredMaxwellianSampler<
                        2, decltype(constant_cross_section)>(
                        (norm::temp_SI * norm::kB) /
                            (child_mass * norm::mass_amu_SI * norm::vel *
                             norm::vel),
                        constant_cross_section, rng_kernel);
                    auto data_calculator =
                        DataCalculator<decltype(data_calc_sampler)>(
                            data_calc_sampler);
                    if (this->ndim == 2)
                    {
                        auto cx_kernel = CXReactionKernels<2>(
                            target_species, projectile_species, prop_map);
                        reaction = std::make_shared<LinearReactionBase<
                            1, decltype(rate_data), decltype(cx_kernel),
                            decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            projectile_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(target_species.get_id())},
                            rate_data, cx_kernel, data_calculator);
                    }
                    else if (this->ndim == 3)
                    {
                        auto cx_kernel = CXReactionKernels<3>(
                            target_species, projectile_species, prop_map);
                        reaction = std::make_shared<LinearReactionBase<
                            1, decltype(rate_data), decltype(cx_kernel),
                            decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            projectile_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(target_species.get_id())},
                            rate_data, cx_kernel, data_calculator);
                    }
                }
                else if (std::get<2>(v).first == "AMJUEL")
                {
                    auto rate_data =
                        AMJUEL::cx_rate_data(parent_mass, child_mass);
                    auto cross_section =
                        AMJUEL::amjuel_fit_cross_section(reduced_mass);

                    auto data_calc_sampler =
                        FilteredMaxwellianSampler<2, decltype(cross_section)>(
                            (norm::temp_SI * norm::kB) /
                                (child_mass * norm::mass_amu_SI * norm::vel *
                                 norm::vel),
                            cross_section, rng_kernel);

                    auto data_calculator =
                        DataCalculator<decltype(data_calc_sampler)>(
                            data_calc_sampler);

                    if (this->ndim == 2)
                    {
                        auto cx_kernel = CXReactionKernels<2>(
                            target_species, projectile_species, prop_map);
                        reaction = std::make_shared<LinearReactionBase<
                            1, decltype(rate_data), decltype(cx_kernel),
                            decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            projectile_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(target_species.get_id())},
                            rate_data, cx_kernel, data_calculator);
                    }
                    else if (this->ndim == 3)
                    {
                        auto cx_kernel = CXReactionKernels<3>(
                            target_species, projectile_species, prop_map);
                        reaction = std::make_shared<LinearReactionBase<
                            1, decltype(rate_data), decltype(cx_kernel),
                            decltype(data_calculator)>>(
                            this->particle_group->sycl_target,
                            projectile_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(target_species.get_id())},
                            rate_data, cx_kernel, data_calculator);
                    }
                }
            }
            this->reaction_controller->add_reaction(reaction);
        }
    }

    // class ReactionsBoundary
    // {

    // };

    // void pre_advection(ParticleSubGroupSharedPtr sg) override
    // {

    // };

    // void apply_boundary_conditions(ParticleSubGroupSharedPtr sg) override
    // {

    // };

    inline void integrate_inner(ParticleSubGroupSharedPtr sg,
                                const double dt_inner) override
    {
        ParticleSystem::integrate_inner(sg, dt_inner);
        reaction_controller->apply_reactions(this->particle_group, dt_inner);
    }

    inline void finish_setup(
        std::vector<std::shared_ptr<DisContField>> &src_fields,
        std::vector<Sym<REAL>> &syms, std::vector<int> &components) override
    {
        this->src_syms       = syms;
        this->src_components = components;

        this->field_project = std::make_shared<FieldProject<DisContField>>(
            src_fields, this->particle_group, this->cell_id_translation);

        auto project_transform = std::make_shared<ProjectTransformation>(
            src_fields, this->src_syms, this->src_components,
            this->particle_group, this->cell_id_translation);
        auto project_transform_wrapper =
            std::make_shared<TransformationWrapper>(
                std::dynamic_pointer_cast<TransformationStrategy>(
                    project_transform));

        auto remove_transform =
            std::make_shared<SimpleRemovalTransformationStrategy>();
        auto remove_transform_wrapper = std::make_shared<TransformationWrapper>(
            std::vector<std::shared_ptr<MarkingStrategy>>{make_marking_strategy<
                ComparisonMarkerSingle<REAL, LessThanComp>>(Sym<REAL>("WEIGHT"),
                                                            1e-10)},
            std::dynamic_pointer_cast<TransformationStrategy>(
                remove_transform));

        std::shared_ptr<TransformationStrategy> merge_transform;

        if (this->ndim == 2)
        {
            merge_transform =
                make_transformation_strategy<MergeTransformationStrategy<2>>();
        }
        else if (this->ndim == 3)
        {
            merge_transform =
                make_transformation_strategy<MergeTransformationStrategy<3>>();
        }

        auto merge_transform_wrapper = std::make_shared<TransformationWrapper>(
            std::vector<std::shared_ptr<MarkingStrategy>>{make_marking_strategy<
                ComparisonMarkerSingle<REAL, LessThanComp>>(Sym<REAL>("WEIGHT"),
                                                            0.01)},
            merge_transform);

        std::vector<std::string> src_names{"ELECTRON_SOURCE_DENSITY",
                                           "ELECTRON_SOURCE_ENERGY",
                                           "ELECTRON_SOURCE_MOMENTUM"};

        for (auto &[k, v] : this->config->get_particle_species())
        {
            std::string name = std::get<0>(v);
            src_names.push_back(name + "_SOURCE_DENSITY");
            src_names.push_back(name + "_SOURCE_ENERGY");
            src_names.push_back(name + "_SOURCE_MOMENTUM");
        }

        auto zeroer_transform =
            make_transformation_strategy<ParticleDatZeroer<REAL>>(src_names);

        this->zeroer_transform_wrapper =
            std::make_shared<TransformationWrapper>(
                std::dynamic_pointer_cast<TransformationStrategy>(
                    zeroer_transform));

        this->reaction_controller = std::make_shared<ReactionController>(
            std::vector<std::shared_ptr<TransformationWrapper>>{
                /*zeroer_transform_wrapper,*/ merge_transform_wrapper,
                remove_transform_wrapper},
            std::vector<std::shared_ptr<TransformationWrapper>>{
                /*project_transform_wrapper,*/ merge_transform_wrapper,
                remove_transform_wrapper});

        set_up_reactions();
        init_output("particle_trajectory.h5part", Sym<REAL>("POSITION"),
                    Sym<INT>("INTERNAL_STATE"), Sym<INT>("CELL_ID"),
                    Sym<REAL>("VELOCITY"), Sym<REAL>("MAGNETIC_FIELD"),
                    Sym<REAL>("ELECTRON_DENSITY"), this->src_syms,
                    Sym<REAL>("WEIGHT"), Sym<INT>("ID"),
                    Sym<REAL>("TOT_REACTION_RATE"), Sym<REAL>("FLUID_DENSITY"),
                    Sym<REAL>("FLUID_TEMPERATURE"));
    }

protected:
    class ProjectTransformation : public TransformationStrategy
    {

    public:
        ProjectTransformation(
            std::vector<std::shared_ptr<DisContField>> &src_fields,
            std::vector<Sym<REAL>> &src_syms, std::vector<int> &src_components,
            ParticleGroupSharedPtr particle_group,
            std::shared_ptr<CellIDTranslation> cell_id_translation)
            : syms(src_syms), components(src_components)
        {
            this->field_project = std::make_shared<FieldProject<DisContField>>(
                src_fields, particle_group, cell_id_translation);
        }

        void transform(ParticleSubGroupSharedPtr sub_group) override
        {
            this->field_project->project(sub_group, syms, components);
        }

    private:
        std::vector<Sym<REAL>> syms;
        std::vector<int> components;

        std::shared_ptr<FieldProject<DisContField>> field_project;
    };

    /// Reaction Controller
    std::shared_ptr<ReactionController> reaction_controller;
};

} // namespace NESO::Solvers::tokamak
#endif