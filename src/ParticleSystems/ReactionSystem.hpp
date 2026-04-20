#ifndef REACTIONSYSTEM_HPP
#define REACTIONSYSTEM_HPP

#include "ParticleSystem.hpp"
#include <reactions/reactions.hpp>

using namespace VANTAGE::Reactions;

namespace PENKNIFE
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

    std::shared_ptr<ParticleDatZeroer<REAL>> zeroer_transform;

    inline void integrate(const double time_end, const double dt) override
    {
        ParticleSystem::integrate(time_end, dt);
    }

    inline void apply_timestep(const double dt) override
    {
        this->boundary->init_funcs();
        ParticleSystem::apply_timestep(dt);
        this->boundary->finalise_funcs();
        if (this->config->get_reactions().size())
        {
            reaction_controller->apply(this->particle_group, dt);
        }
    }

    inline void zero_source_dats() override
    {
        this->zeroer_transform->transform(
            particle_sub_group(this->particle_group));
    }

    void set_up_reactions();
    void set_up_boundaries() override;

    inline void pre_advection(ParticleSubGroupSharedPtr sg) override
    {
        this->boundary->pre_advection(sg);
    };

    inline void apply_boundary_conditions(ParticleSubGroupSharedPtr sg,
                                          ParticleGroupSharedPtr cg,
                                          double dt) override
    {
        this->boundary->execute(sg, cg, dt);
    };

    inline void integrate_inner(ParticleSubGroupSharedPtr sg,
                                const double dt_inner) override
    {
        ParticleSystem::integrate_inner(sg, dt_inner);
    }

    void finish_setup(std::vector<std::shared_ptr<DisContField>> &src_fields,
                      std::vector<Sym<REAL>> &syms,
                      std::vector<int> &components) override;

    class ReactionsBoundary
    {

    public:
        ReactionsBoundary(
            Sym<REAL> time_step_prop_sym, SYCLTargetSharedPtr sycl_target,
            std::shared_ptr<ParticleMeshInterface> mesh,
            NESOReaderSharedPtr config,
            std::map<std::string, SpeciesInfo> &species,
            ParameterStoreSharedPtr store = std::make_shared<ParameterStore>());

        inline void pre_advection(ParticleSubGroupSharedPtr particle_sub_group)
        {
            this->composite_intersection->pre_integration(particle_sub_group);
        }

        inline void init_funcs()
        {

            for (int b_id : this->boundary_ids)
            {
                this->composite_intersection->function_project_initialise(
                    this->funcs_dens[b_id]);
                this->composite_intersection->function_project_initialise(
                    this->funcs_energy[b_id]);
            }
        }

        inline void finalise_funcs()
        {
            for (int b_id : this->boundary_ids)
            {
                this->composite_intersection->function_project_finalise(
                    this->funcs_dens[b_id]);
                this->composite_intersection->function_project_finalise(
                    this->funcs_energy[b_id]);
            }
        }

        inline void execute(ParticleSubGroupSharedPtr particle_sub_group,
                            ParticleGroupSharedPtr child_group, double dt)
        {
            NESOASSERT(this->ndim == 3 || this->ndim == 2,
                       "Unexpected number of dimensions.");

            auto groups = this->composite_intersection->get_intersections(
                particle_sub_group);

            for (auto &[id, sg] : groups)
            {
                copy_ephemeral_dat_to_particle_dat(
                    sg, Sym<REAL>("NESO_PARTICLES_BOUNDARY_INTERSECTION_POINT"),
                    Sym<REAL>("NESO_PARTICLES_BOUNDARY_INTERSECTION_POINT"));
                copy_ephemeral_dat_to_particle_dat(
                    sg, Sym<REAL>("NESO_PARTICLES_BOUNDARY_NORMAL"),
                    Sym<REAL>("NESO_PARTICLES_BOUNDARY_NORMAL"));
                copy_ephemeral_dat_to_particle_dat(
                    sg, Sym<INT>("NESO_PARTICLES_BOUNDARY_METADATA"),
                    Sym<INT>("NESO_PARTICLES_BOUNDARY_METADATA"));

                this->boundary_truncation->execute(
                    sg,
                    get_particle_group(particle_sub_group)->position_dat->sym,
                    this->time_step_prop_sym,
                    this->composite_intersection->previous_position_sym);
                this->reaction_controllers[id]->apply(
                    sg, dt, child_group, ControllerMode::surface_mode);

                this->composite_intersection->function_project_contribute(
                    sg, Sym<REAL>("SURFACE_DENSITY_SOURCE"), 0, false,
                    funcs_dens[id]);
                this->composite_intersection->function_project_contribute(
                    sg, Sym<REAL>("SURFACE_ENERGY_SOURCE"), 0, false,
                    funcs_energy[id]);
            }
            remove_wrapper->transform(particle_sub_group);
        }

    private:
        Sym<REAL> time_step_prop_sym;

        SYCLTargetSharedPtr sycl_target;
        std::shared_ptr<CompositeInteraction::CompositeIntersection>
            composite_intersection;
        std::shared_ptr<BoundaryTruncation> boundary_truncation;

        std::vector<int> boundary_ids;

        std::map<int, std::shared_ptr<ReactionController>> reaction_controllers;
        std::shared_ptr<TransformationWrapper> remove_wrapper;

        const int ndim;
        const int vdim;
        REAL reset_distance;

        NESOReaderSharedPtr config;

        std::map<int, CompositeInteraction::CompositeFunctionSharedPtr>
            funcs_dens;
        std::map<int, CompositeInteraction::CompositeFunctionSharedPtr>
            funcs_energy;
    };

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

    std::shared_ptr<ReactionsBoundary> boundary;
};

} // namespace PENKNIFE
#endif