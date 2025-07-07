#ifndef TOKAMAKPARTICLE_SYSTEM_H
#define TOKAMAKPARTICLE_SYSTEM_H

#include <array>

#include <nektar_interface/function_evaluation.hpp>
#include <nektar_interface/function_projection.hpp>
#include <nektar_interface/particle_boundary_conditions.hpp>
#include <nektar_interface/particle_cell_mapping/particle_cell_mapping_common.hpp>
#include <nektar_interface/solver_base/partsys_base.hpp>
#include <nektar_interface/utilities.hpp>
#include <neso_particles.hpp>

// namespace LU = Nektar::LibUtilities;
// namespace NP = NESO::Particles;

namespace NESO::Solvers::tokamak
{

/**
 * @brief
 */
class ParticleSystem : public PartSysBase
{

public:
    static std::string class_name;
    /**
     * @brief Create an instance of this class and initialise it.
     */
    static ParticleSystemSharedPtr create(const NESOReaderSharedPtr &session,
                                          const SD::MeshGraphSharedPtr &graph)
    {
        ParticleSystemSharedPtr p =
            MemoryManager<ParticleSystem>::AllocateSharedPtr(session, graph);
        return p;
    }

    /**
     *  Create a new instance.
     *
     *  @param session Particle reader to use for parameters and simulation
     * specification.
     *  @param graph Nektar++ MeshGraph on which particles exist.
     *  @param comm (optional) MPI communicator to use - default MPI_COMM_WORLD.
     *
     */
    ParticleSystem(NESOReaderSharedPtr session, SD::MeshGraphSharedPtr graph,
                   MPI_Comm comm = MPI_COMM_WORLD);

    virtual ~ParticleSystem() override = default;

    /// Disable (implicit) copies.
    ParticleSystem(const ParticleSystem &st) = delete;
    /// Disable (implicit) copies.
    ParticleSystem &operator=(ParticleSystem const &a) = delete;

    inline virtual void init_spec() override
    {
        this->particle_spec = {
            ParticleProp(Sym<REAL>("POSITION"), this->ndim, true),
            ParticleProp(Sym<REAL>("VELOCITY"), this->ndim),
            ParticleProp(Sym<INT>("CELL_ID"), 1, true),
            ParticleProp(Sym<INT>("ID"), 1),
            ParticleProp(Sym<INT>("INTERNAL_STATE"), 1),
            ParticleProp(Sym<REAL>("M"), 1),
            ParticleProp(Sym<REAL>("Q"), 1),
            ParticleProp(Sym<REAL>("ELECTRON_DENSITY"), 1),
            ParticleProp(Sym<REAL>("ELECTRON_TEMPERATURE"), 1),
            ParticleProp(Sym<REAL>("ELECTRON_SOURCE_ENERGY"), 1),
            ParticleProp(Sym<REAL>("ELECTRON_SOURCE_MOMENTUM"), this->ndim),
            ParticleProp(Sym<REAL>("ELECTRON_SOURCE_DENSITY"), 1),
            ParticleProp(Sym<REAL>("ELECTRIC_FIELD"), 3),
            ParticleProp(Sym<REAL>("MAGNETIC_FIELD"), 3),
            ParticleProp(Sym<REAL>("TSP"), 2)};

        for (auto &[k, v] : this->config->get_particle_species())
        {
            std::string name = std::get<0>(v);
            this->particle_spec.push(
                ParticleProp(Sym<REAL>(name + "_SOURCE_DENSITY"), 1));
            this->particle_spec.push(
                ParticleProp(Sym<REAL>(name + "_SOURCE_ENERGY"), 1));
            this->particle_spec.push(
                ParticleProp(Sym<REAL>(name + "_SOURCE_MOMENTUM"), this->ndim));
        }
        this->particle_spec.push(ParticleProp(Sym<REAL>("WEIGHT"), 1));
        this->particle_spec.push(
            ParticleProp(Sym<REAL>("TOT_REACTION_RATE"), 1));
        this->particle_spec.push(
            ParticleProp(Sym<INT>("REACTIONS_PANIC_FLAG"), 1));
        this->particle_spec.push(
            ParticleProp(Sym<INT>("PARTICLE_REACTED_FLAG"), 1));
        this->particle_spec.push(ParticleProp(Sym<REAL>("FLUID_DENSITY"), 1));
        this->particle_spec.push(
            ParticleProp(Sym<REAL>("FLUID_TEMPERATURE"), 1));
        this->particle_spec.push(
            ParticleProp(Sym<REAL>("FLUID_FLOW_SPEED"), this->ndim));
    }

    virtual void init_object() override
    {
        PartSysBase::init_object();
        this->particle_remover =
            std::make_shared<ParticleRemover>(this->sycl_target);

        this->transfer_particles();
        // for (auto &[k, v] : this->get_species())
        //{
        pre_advection(particle_sub_group(this->particle_group));
        apply_boundary_conditions(particle_sub_group(this->particle_group));
        //}
    }

    virtual void set_up_species() override;

    struct SpeciesInfo
    {
        std::string name;
        double mass;
        double charge;

        std::shared_ptr<ParticleSubGroup> sub_group;

        /// Reflective Boundary Conditions
        // std::shared_ptr<NektarCompositeTruncatedReflection> reflection;
        // std::vector<int> reflection_composites;
    };
    virtual std::map<int, SpeciesInfo> &get_species()
    {
        return species_map;
    }

    void set_up_boundaries();

    /**
     *  Integrate the particle system forward to the requested time using
     *  (at most) the requested time step.
     *
     *  @param time_end Target time to integrate to.
     *  @param dt Time step size.
     */
    inline virtual void integrate(const double time_end, const double dt)
    {

        // Get the current simulation time.
        NESOASSERT(time_end >= this->simulation_time,
                   "Cannot integrate backwards in time.");
        if (time_end == this->simulation_time)
        {
            return;
        }

        double time_tmp = this->simulation_time;
        while (time_tmp < time_end)
        {
            const double dt_inner = std::min(dt, time_end - time_tmp);
            this->add_sources(time_tmp, dt_inner);
            this->add_sinks(time_tmp, dt_inner);
            this->apply_timestep(particle_sub_group(this->particle_group),
                                 dt_inner);
            this->transfer_particles();

            time_tmp += dt_inner;
        }

        this->simulation_time = time_end;
    }

    /**
     * Setup the projection object to use the src fields from the eqnsys.
     *
     * @param src_fields Nektar++ fields to project ionised particle data onto.
     * @param syms Corresponding Particle Syms
     * @param syms Corresponding components
     */
    inline virtual void finish_setup(
        std::vector<std::shared_ptr<DisContField>> &src_fields,
        std::vector<Sym<REAL>> &syms, std::vector<int> &components)
    {
        this->src_syms       = syms;
        this->src_components = components;
        this->field_project  = std::make_shared<FieldProject<DisContField>>(
            src_fields, this->particle_group, this->cell_id_translation);
        init_output("particle_trajectory.h5part", Sym<REAL>("POSITION"),
                    Sym<INT>("INTERNAL_STATE"), Sym<INT>("CELL_ID"),
                    Sym<REAL>("VELOCITY"), Sym<REAL>("MAGNETIC_FIELD"),
                    Sym<REAL>("ELECTRON_DENSITY"), this->src_syms,
                    Sym<INT>("ID"), Sym<REAL>("TOT_REACTION_RATE"));
    }

    void add_sources(double time, double dt);
    void add_sinks(double time, double dt);

    /**
     *  Project the plasma source and momentum contributions from particle data
     *  onto field data.
     */
    inline virtual void project_source_terms()
    {
        NESOASSERT(this->field_project != nullptr,
                   "Field project object is null. Was setup_project called?");

        this->field_project->project(this->particle_group, this->src_syms,
                                     this->src_components);

        // remove fully ionised particles from the simulation
        remove_marked_particles();
    }

    inline virtual void setup_evaluate_fields(
        Array<OneD, std::shared_ptr<DisContField>> &E,
        Array<OneD, std::shared_ptr<DisContField>> &B,
        std::shared_ptr<DisContField> ne, std::shared_ptr<DisContField> Te,
        Array<OneD, std::shared_ptr<DisContField>> &ve)
    {
        auto mesh = std::dynamic_pointer_cast<ParticleMeshInterface>(
            particle_group->domain->mesh);
        this->field_evaluate_ne =
            std::make_shared<FunctionEvaluateBasis<DisContField>>(
                ne, mesh, this->cell_id_translation);
        if (Te)
        {
            this->field_evaluate_Te =
                std::make_shared<FunctionEvaluateBasis<DisContField>>(
                    Te, mesh, this->cell_id_translation);
        }
        this->field_evaluate_ve =
            std::vector<std::shared_ptr<FunctionEvaluateBasis<DisContField>>>(
                this->ndim);
        for (int d = 0; d < this->ndim; ++d)
        {
            if (ve[d])
            {
                this->field_evaluate_ve[d] =
                    std::make_shared<FunctionEvaluateBasis<DisContField>>(
                        ve[d], mesh, this->cell_id_translation);
            }

            this->field_evaluate_E.emplace_back(
                std::make_shared<FunctionEvaluateBasis<DisContField>>(
                    E[d], mesh, this->cell_id_translation));
            this->field_evaluate_B.emplace_back(
                std::make_shared<FunctionEvaluateBasis<DisContField>>(
                    B[d], mesh, this->cell_id_translation));
        }
    }

    /**
     * Evaluate E and B at the particle locations.
     */
    inline virtual void evaluate_fields(
        Array<OneD, std::shared_ptr<DisContField>> &E,
        Array<OneD, std::shared_ptr<DisContField>> &B,
        std::shared_ptr<DisContField> ne, std::shared_ptr<DisContField> Te,
        Array<OneD, std::shared_ptr<DisContField>> &ve)
    {

        for (int d = 0; d < this->ndim; ++d)
        {
            NESOASSERT(this->field_evaluate_E[d] != nullptr,
                       "FieldEvaluate not setup.");
            this->field_evaluate_E[d]->evaluate(this->particle_group,
                                                Sym<REAL>("ELECTRIC_FIELD"), d,
                                                E[d]->GetCoeffs());
            NESOASSERT(this->field_evaluate_B[d] != nullptr,
                       "FieldEvaluate not setup.");
            this->field_evaluate_B[d]->evaluate(this->particle_group,
                                                Sym<REAL>("MAGNETIC_FIELD"), d,
                                                E[d]->GetCoeffs());
        }

        NESOASSERT(this->field_evaluate_ne != nullptr,
                   "FieldEvaluate not setup.");
        this->field_evaluate_ne->evaluate(this->particle_group,
                                          Sym<REAL>("FLUID_DENSITY"), 0,
                                          ne->GetCoeffs());
        if (field_evaluate_Te)                
        {
            this->field_evaluate_Te->evaluate(this->particle_group,
                                              Sym<REAL>("FLUID_TEMPERATURE"), 0,
                                              Te->GetCoeffs());
        }
        for (int d = 0; d < this->ndim; ++d)
        {
            if (this->field_evaluate_ve[d])
            {
                this->field_evaluate_ve[d]->evaluate(
                    this->particle_group, Sym<REAL>("FLUID_FLOW_SPEED"), d,
                    ve[d]->GetCoeffs());
            }
        }
    }

    inline void remove_marked_particles()
    {
        this->particle_remover->remove(this->particle_group,
                                       (*this->particle_group)[Sym<INT>("ID")],
                                       this->particle_remove_key);
    }

    /**
     * Adds a velocity drift of \alpha * V_drift where V_drift = E X B
     * evaluation to the velocity of each particle. The coefficient \alpha is
     * the read from the session file key `particle_v_drift_scaling`.
     */
    inline void initialise_particles_from_fields(
        Array<OneD, std::shared_ptr<DisContField>> &E,
        Array<OneD, std::shared_ptr<DisContField>> &B,
        std::shared_ptr<DisContField> ne, std::shared_ptr<DisContField> Te,
        Array<OneD, std::shared_ptr<DisContField>> &ve)
    {
        double h_alpha;
        this->config->load_parameter("particle_v_drift_scaling", h_alpha, 1.0);
        const double k_alpha = h_alpha;

        this->evaluate_fields(E, B, ne, Te, ve);
        if (this->ndim == 3)
        {
            particle_loop(
                "ParticleSystem:initialise_particles_from_fields",
                this->particle_group,
                [=](auto VELOCITY, auto Ek, auto Bk)
                {
                    // Ei contains d(phi)/dx_i.
                    const auto mE0 = Ek.at(0);
                    const auto mE1 = Ek.at(1);
                    const auto mE2 = Ek.at(2);
                    REAL exb0, exb1, exb2;
                    MAPPING_CROSS_PRODUCT_3D(mE0, mE1, mE2, Bk.at(0), Bk.at(1),
                                             Bk.at(2), exb0, exb1, exb2);
                    VELOCITY.at(0) += k_alpha * exb0;
                    VELOCITY.at(1) += k_alpha * exb1;
                    VELOCITY.at(2) += k_alpha * exb2;
                },
                Access::write(Sym<REAL>("VELOCITY")),
                Access::read(Sym<REAL>("ELECTRIC_FIELD")),
                Access::read(Sym<REAL>("MAGNETIC_FIELD")))
                ->execute();
        }
        else if (this->ndim == 2)
        {
            particle_loop(
                "ParticleSystem:initialise_particles_from_fields",
                this->particle_group,
                [=](auto VELOCITY, auto Ek, auto Bk)
                {
                    // Ei contains d(phi)/dx_i.
                    const auto mE0 = Ek.at(0);
                    const auto mE1 = Ek.at(1);

                    VELOCITY.at(0) += k_alpha * mE0;
                    VELOCITY.at(1) += k_alpha * mE1;
                },
                Access::write(Sym<REAL>("VELOCITY")),
                Access::read(Sym<REAL>("ELECTRIC_FIELD")),
                Access::read(Sym<REAL>("MAGNETIC_FIELD")))
                ->execute();
        }
    }

protected:
    virtual inline void integrate_inner(ParticleSubGroupSharedPtr sg,
                                        const double dt_inner)
    {
        const auto k_dt = dt_inner;
        if (this->ndim == 3)
        {
            particle_loop(
                "ParticleSystem:boris", sg,
                [=](auto E, auto B, auto Q, auto M, auto P, auto V, auto TSP)
                {
                    const REAL dt_left  = k_dt - TSP.at(0);
                    const REAL hdt_left = dt_left * 0.5;
                    if (dt_left > 0.0)
                    {
                        const REAL QoM = Q.at(0) / M.at(0);

                        const REAL scaling_t = QoM * hdt_left;
                        const REAL t_0       = B.at(0) * scaling_t;
                        const REAL t_1       = B.at(1) * scaling_t;
                        const REAL t_2       = B.at(2) * scaling_t;

                        const REAL tmagsq = t_0 * t_0 + t_1 * t_1 + t_2 * t_2;
                        const REAL scaling_s = 2.0 / (1.0 + tmagsq);

                        const REAL s_0 = scaling_s * t_0;
                        const REAL s_1 = scaling_s * t_1;
                        const REAL s_2 = scaling_s * t_2;

                        const REAL V_0 = V.at(0);
                        const REAL V_1 = V.at(1);
                        const REAL V_2 = V.at(2);

                        const REAL v_minus_0 = V_0 + (E.at(0)) * scaling_t;
                        const REAL v_minus_1 = V_1 + (E.at(1)) * scaling_t;
                        const REAL v_minus_2 = V_2 + (E.at(2)) * scaling_t;

                        REAL v_prime_0, v_prime_1, v_prime_2;
                        MAPPING_CROSS_PRODUCT_3D(
                            v_minus_0, v_minus_1, v_minus_2, t_0, t_1, t_2,
                            v_prime_0, v_prime_1, v_prime_2)

                        v_prime_0 += v_minus_0;
                        v_prime_1 += v_minus_1;
                        v_prime_2 += v_minus_2;

                        REAL v_plus_0, v_plus_1, v_plus_2;
                        MAPPING_CROSS_PRODUCT_3D(v_prime_0, v_prime_1,
                                                 v_prime_2, s_0, s_1, s_2,
                                                 v_plus_0, v_plus_1, v_plus_2)

                        v_plus_0 += v_minus_0;
                        v_plus_1 += v_minus_1;
                        v_plus_2 += v_minus_2;

                        V.at(0) = v_plus_0 + scaling_t * (E.at(0));
                        V.at(1) = v_plus_1 + scaling_t * (E.at(1));
                        V.at(2) = v_plus_2 + scaling_t * (E.at(2));

                        // update of position to next time step
                        P.at(0) += dt_left * V.at(0);
                        P.at(1) += dt_left * V.at(1);
                        P.at(2) += dt_left * V.at(2);
                        TSP.at(0) = k_dt;
                        TSP.at(1) = dt_left;
                    }
                },
                Access::read(Sym<REAL>("ELECTRIC_FIELD")),
                Access::read(Sym<REAL>("MAGNETIC_FIELD")),
                Access::read(Sym<REAL>("Q")), Access::read(Sym<REAL>("M")),
                Access::write(Sym<REAL>("POSITION")),
                Access::write(Sym<REAL>("VELOCITY")),
                Access::write(Sym<REAL>("TSP")))
                ->execute();
        }
        else if (ndim == 2)
        {
            particle_loop(
                "euler_advection", sg,
                [=](auto V, auto P, auto TSP)
                {
                    const REAL dt_left = k_dt - TSP.at(0);
                    if (dt_left > 0.0)
                    {

                        P.at(0) += dt_left * V.at(0);
                        P.at(1) += dt_left * V.at(1);

                        TSP.at(0) = k_dt;
                        TSP.at(1) = dt_left;
                    }
                },
                Access::read(Sym<REAL>("VELOCITY")),
                Access::write(Sym<REAL>("POSITION")),
                Access::write(Sym<REAL>("TSP")))
                ->execute();
        }
    };

    uint64_t total_num_particles_added = 0;

    const int particle_remove_key = -1;
    std::shared_ptr<ParticleRemover> particle_remover;

    std::map<int, SpeciesInfo> species_map;

    std::vector<Sym<REAL>> src_syms;
    std::vector<int> src_components;

    std::shared_ptr<FieldProject<DisContField>> field_project;

    std::shared_ptr<FunctionEvaluateBasis<DisContField>> field_evaluate_ne;
    std::shared_ptr<FunctionEvaluateBasis<DisContField>> field_evaluate_Te;
    std::vector<std::shared_ptr<FunctionEvaluateBasis<DisContField>>>
        field_evaluate_ve;

    std::vector<std::shared_ptr<FunctionEvaluateBasis<DisContField>>>
        field_evaluate_E;

    std::vector<std::shared_ptr<FunctionEvaluateBasis<DisContField>>>
        field_evaluate_B;

    std::shared_ptr<NektarCompositeTruncatedReflection> reflection;

    /// Simulation time
    double simulation_time;

    inline void apply_timestep_reset(ParticleSubGroupSharedPtr sg)
    {
        particle_loop(
            sg,
            [=](auto TSP)
            {
                TSP.at(0) = 0.0;
                TSP.at(1) = 0.0;
            },
            Access::write(Sym<REAL>("TSP")))
            ->execute();
    }

    void pre_advection(ParticleSubGroupSharedPtr sg)
    {
        // if (auto r = this->species_map[k].reflection)
        //{

        reflection->pre_advection(sg);
        //}
    };

    void apply_boundary_conditions(ParticleSubGroupSharedPtr sg)
    {
        // if (auto r = this->species_map[k].reflection)
        //{
        reflection->execute(sg);
        //}
    };

    auto find_partial_moves(ParticleSubGroupSharedPtr sg, const double dt)
    {
        return particle_sub_group(
            sg, [=](auto TSP) { return TSP.at(0) < dt; },
            Access::read(Sym<REAL>("TSP")));
    };

    bool partial_moves_remaining(ParticleSubGroupSharedPtr sg)
    {
        const int size = sg->get_npart_local();
        int size_global;
        MPICHK(MPI_Allreduce(&size, &size_global, 1, MPI_INT, MPI_SUM,
                             sycl_target->comm_pair.comm_parent));
        return size_global > 0;
    };

    inline void apply_timestep(ParticleSubGroupSharedPtr sg, const double dt)
    {
        apply_timestep_reset(sg);
        pre_advection(sg);
        integrate_inner(sg, dt);
        apply_boundary_conditions(sg);
        sg = find_partial_moves(sg, dt);
        while (partial_moves_remaining(sg))
        {
            pre_advection(sg);
            integrate_inner(sg, dt);
            apply_boundary_conditions(sg);
            sg = find_partial_moves(sg, dt);
        }
    }

    /**
     *  Apply boundary conditions and transfer particles between MPI ranks.
     * // Move some of this to PartSysBase / make it a pure-virtual func?
     */
    inline void transfer_particles()
    {
        auto t0 = profile_timestamp();

        this->particle_group->hybrid_move();
        this->cell_id_translation->execute();
        this->particle_group->cell_move();
        this->sycl_target->profile_map.inc(
            "ParticleSystem", "transfer_particles", 1,
            profile_elapsed(t0, profile_timestamp()));
    }
};
} // namespace NESO::Solvers::tokamak
#endif
