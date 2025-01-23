#ifndef REACTIONSYSTEM_HPP
#define REACTIONSYSTEM_HPP

#include "ParticleSystem.hpp"
#include <reactions.hpp>

using namespace Reactions;

namespace NESO::Solvers::tokamak
{
class ReactionSystem : public ParticleSystem
{

public:
    static std::string class_name;
    static ParticleSystemSharedPtr create(const ParticleReaderSharedPtr session,
                                          const SD::MeshGraphSharedPtr graph)
    {
        ParticleSystemSharedPtr p =
            MemoryManager<ReactionSystem>::AllocateSharedPtr(session, graph);
        return p;
    }

    ReactionSystem(ParticleReaderSharedPtr session,
                   SD::MeshGraphSharedPtr graph);

    ~ReactionSystem() override = default;

    inline void init_spec() override
    {
        ParticleSystem::init_spec();
        this->particle_spec.push(ParticleProp(Sym<REAL>("TOT_REACTION_RATE"), 1));
        this->particle_spec.push(ParticleProp(Sym<REAL>("WEIGHT"), 1));
    }

    inline void set_up_reactions()
    {
        auto prop_map = default_map;

        for (const auto &[k, v] : session->get_reactions())
        {
            if (std::get<0>(v) == "Ionisation")
            {

                auto electron_species = Species("ELECTRON", 5.5e-4, -1.0);
                auto target_species   = Species(
                    species_map[std::get<1>(v)[0]].name, std::get<1>(v)[0],
                    species_map[std::get<1>(v)[0]].mass,
                    species_map[std::get<1>(v)[0]].charge);

                auto test_data = FixedRateData(1.0);
                auto reaction =
                    ElectronImpactIonisation<FixedRateData, FixedRateData>(
                        particle_group->sycl_target,
                        Sym<REAL>(
                            prop_map[default_properties.tot_reaction_rate]),
                        Sym<REAL>(prop_map[default_properties.weight]),
                        test_data, test_data, target_species, electron_species,
                        particle_spec);

                reaction_controller->add_reaction(
                    std::make_shared<decltype(reaction)>(reaction));
            }
            else if (std::get<0>(v) == "Recombination")
            {
                NESOASSERT(false, "Recombination kernels not yet implemented");
                /*
                auto electron_species = Species("ELECTRON", 5.5e-4, -1.0);
                auto species_null     = Species("", 1.0, 0.0, -1);
                auto target_species =
                Species(species_map[std::get<1>(v)[0]].name,
                            std::get<1>(v)[0],
                            species_map[std::get<1>(v)[0]].mass,
                            species_map[std::get<1>(v)[0]].charge);
                auto recomb_reaction_kernel = RecombReactionKernels<>(
                    target_species, electron_species, prop_map);
                std::array<int, 1> recomb_out_states = {0};
                auto recomb_data                     = FixedRateData(1.0);

                auto vx_beam_data = FixedRateData(1.0);
                auto vy_beam_data = FixedRateData(-1.0);

                auto data_calculator =
                    DataCalculator<FixedRateData, FixedRateData>(
                        particle_spec, vx_beam_data, vy_beam_data);
                auto reaction = LinearReactionBase<
                    1, decltype(recomb_data), decltype(recomb_reaction_kernel),
                    decltype(data_calculator)>(
                    particle_group->sycl_target,
                    Sym<REAL>(prop_map[default_properties.tot_reaction_rate]),
                    Sym<REAL>(prop_map[default_properties.weight]),
                    species_null.get_id(), recomb_out_states, recomb_data,
                    recomb_reaction_kernel, particle_spec, data_calculator);

                reaction_controller->add_reaction(
                    std::make_shared<decltype(reaction)>(reaction));
                    */
            }

            else if (std::get<0>(v) == "ChargeExchange")
            {
                auto projectile_species = Species(
                    species_map[std::get<1>(v)[0]].name, std::get<1>(v)[0],
                    species_map[std::get<1>(v)[0]].mass,
                    species_map[std::get<1>(v)[0]].charge);
                auto target_species = Species(
                    species_map[std::get<1>(v)[1]].name, std::get<1>(v)[1],
                    species_map[std::get<1>(v)[1]].mass,
                    species_map[std::get<1>(v)[1]].charge);

                auto cx_kernel = CXReactionKernels<2>(
                    target_species, projectile_species, prop_map);
                auto rate_data    = FixedRateData(1.0);
                auto vx_beam_data = FixedRateData(1.0);
                auto vy_beam_data = FixedRateData(-1.0);

                auto data_calculator =
                    DataCalculator<FixedRateData, FixedRateData>(
                        particle_spec, vx_beam_data, vy_beam_data);
                auto reaction = LinearReactionBase<
                    1, FixedRateData, CXReactionKernels<2>,
                    DataCalculator<FixedRateData, FixedRateData>>(
                    particle_group->sycl_target,
                    Sym<REAL>(prop_map[default_properties.tot_reaction_rate]),
                    Sym<REAL>(prop_map[default_properties.weight]),
                    projectile_species.get_id(),
                    std::array<int, 1>{target_species.get_id()}, rate_data,
                    cx_kernel, particle_spec, data_calculator);

                reaction_controller->add_reaction(
                    std::make_shared<decltype(reaction)>(reaction));
            }
        }
    }

    inline void integrate_inner(const double dt_inner) override
    {
        ParticleSystem::integrate_inner(dt_inner);
        reaction_controller->apply_reactions(this->particle_group, dt_inner);
    }

    void set_up_particles() override
    {
        ParticleSystem::set_up_particles();
        set_up_reactions();
    }

protected:
    class ProjectTransformation : public TransformationStrategy
    {

    public:
        ProjectTransformation(
            std::vector<std::shared_ptr<DisContField>> &src_fields,
            ParticleGroupSharedPtr particle_group,
            std::shared_ptr<CellIDTranslation> cell_id_translation)
        {
            this->field_project = std::make_shared<FieldProject<DisContField>>(
                src_fields, particle_group, cell_id_translation);
        }

        void transform(ParticleSubGroupSharedPtr sub_group)
        {
            this->field_project->project(sub_group, syms, components);
        }

    private:
        std::vector<Sym<REAL>> syms = {Sym<REAL>("SOURCE_DENSITY"),
                                       Sym<REAL>("SOURCE_ENERGY")};
        std::vector<int> components = {0, 0};

        std::shared_ptr<FieldProject<DisContField>> field_project;
    };

    std::shared_ptr<TransformationWrapper> project_transform_wrapper;
    std::shared_ptr<TransformationWrapper> merge_transform_wrapper;
    std::shared_ptr<TransformationWrapper> remove_transform_wrapper;

    /// Reaction Controller
    std::shared_ptr<ReactionController> reaction_controller;

    std::vector<DisContFieldSharedPtr> src_fields;
};

} // namespace NESO::Solvers::tokamak
#endif