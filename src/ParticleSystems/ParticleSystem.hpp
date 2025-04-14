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
    static ParticleSystemSharedPtr create(
        const ParticleReaderSharedPtr &session,
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
    ParticleSystem(ParticleReaderSharedPtr session,
                   SD::MeshGraphSharedPtr graph,
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
            ParticleProp(Sym<REAL>("ELECTRON_SOURCE_ENERGY"), 1),
            ParticleProp(Sym<REAL>("ELECTRON_SOURCE_MOMENTUM"), this->ndim),
            ParticleProp(Sym<REAL>("ELECTRON_SOURCE_DENSITY"), 1),
            ParticleProp(Sym<REAL>("E0"), 1),
            ParticleProp(Sym<REAL>("E1"), 1),
            ParticleProp(Sym<REAL>("E2"), 1),
            ParticleProp(Sym<REAL>("B0"), 1),
            ParticleProp(Sym<REAL>("B1"), 1),
            ParticleProp(Sym<REAL>("B2"), 1),
            ParticleProp(Sym<REAL>("TSP"), 2)};

        for (auto &[k, v] : this->config->get_species())
        {
            std::string name = std::get<0>(v);
            this->particle_spec.push(
                ParticleProp(Sym<REAL>(name + "_SOURCE_DENSITY"), 1));
            this->particle_spec.push(
                ParticleProp(Sym<REAL>(name + "_SOURCE_ENERGY"), 1));
            this->particle_spec.push(
                ParticleProp(Sym<REAL>(name + "_SOURCE_MOMENTUM"), this->ndim));
        }
    }

    virtual void init_object() override
    {
        this->config->read_particles();
        this->init_spec();
        this->read_params();
        this->particle_group = std::make_shared<ParticleGroup>(
            this->domain, this->particle_spec, this->sycl_target);
        this->cell_id_translation = std::make_shared<CellIDTranslation>(
            this->sycl_target, this->particle_group->cell_id_dat,
            this->particle_mesh_interface);
        this->particle_remover =
            std::make_shared<ParticleRemover>(this->sycl_target);
        this->set_up_species();
        this->set_up_boundaries();
        // PartSysBase::init_object();

        // parallel_advection_initialisation(this->particle_group);
        // parallel_advection_store(this->particle_group);
        // const int num_steps = 20;
        // for (int stepx = 0; stepx < num_steps; stepx++)
        // {

        //     parallel_advection_step(this->particle_group, num_steps, stepx);
        //     this->transfer_particles();
        // }
        // parallel_advection_restore(this->particle_group);
        //  Move particles to the owning ranks and correct cells.
        this->transfer_particles();
        pre_advection(particle_sub_group(this->particle_group));
        apply_boundary_conditions(particle_sub_group(this->particle_group));
    }
    virtual void set_up_particles() override
    {
        PartSysBase::set_up_particles();
    }

    virtual void set_up_species() override;

    struct SpeciesInfo
    {
        std::string name;
        double mass;
        double charge;

        std::shared_ptr<ParticleSubGroup> sub_group;
    };
    virtual std::map<int, SpeciesInfo> &get_species()
    {
        return species_map;
    }

    virtual void set_up_boundaries() override;

    /**
     *  Integrate the particle system forward to the requested time using
     *  (at most) the requested time step.
     *
     *  @param time_end Target time to integrate to.
     *  @param dt Time step size.
     */
    inline void integrate(const double time_end, const double dt)
    {

        // Get the current simulation time.
        NESOASSERT(time_end >= this->simulation_time,
                   "Cannot integrate backwards in time.");
        if (time_end == this->simulation_time || this->num_parts_tot == 0)
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
                    Sym<REAL>("VELOCITY"), Sym<REAL>("B0"), Sym<REAL>("B1"),
                    Sym<REAL>("B2"), Sym<REAL>("ELECTRON_DENSITY"),
                    this->src_syms, Sym<INT>("ID"));
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

    inline void setup_evaluate_ne(std::shared_ptr<DisContField> n)
    {
        this->field_evaluate_ne = std::make_shared<FieldEvaluate<DisContField>>(
            n, this->particle_group, this->cell_id_translation);
        this->fields["ne"] = n;
    }

    /**
     * Set up the evaluation of E.
     *
     * @param E0 Nektar++ field storing E in direction 0.
     * @param E1 Nektar++ field storing E in direction 1.
     * @param E2 Nektar++ field storing E in direction 2.
     */
    inline void setup_evaluate_E(std::shared_ptr<DisContField> E0,
                                 std::shared_ptr<DisContField> E1,
                                 std::shared_ptr<DisContField> E2)
    {
        // TODO redo with redone Bary Evaluate.
        this->field_evaluate_E0 = std::make_shared<FieldEvaluate<DisContField>>(
            E0, this->particle_group, this->cell_id_translation);
        this->field_evaluate_E1 = std::make_shared<FieldEvaluate<DisContField>>(
            E1, this->particle_group, this->cell_id_translation);
        this->field_evaluate_E2 = std::make_shared<FieldEvaluate<DisContField>>(
            E2, this->particle_group, this->cell_id_translation);

        this->fields["E0"] = E0;
        this->fields["E1"] = E1;
        this->fields["E2"] = E2;
    }

    /**
     * Set up the evaluation of B.
     *
     * @param B0 Nektar++ field storing B in direction 0.
     * @param B1 Nektar++ field storing B in direction 1.
     * @param B2 Nektar++ field storing B in direction 2.
     */
    inline void setup_evaluate_B(std::shared_ptr<DisContField> B0,
                                 std::shared_ptr<DisContField> B1,
                                 std::shared_ptr<DisContField> B2)
    {
        // TODO redo with redone Bary Evaluate.
        this->field_evaluate_B0 = std::make_shared<FieldEvaluate<DisContField>>(
            B0, this->particle_group, this->cell_id_translation);
        this->field_evaluate_B1 = std::make_shared<FieldEvaluate<DisContField>>(
            B1, this->particle_group, this->cell_id_translation);
        this->field_evaluate_B2 = std::make_shared<FieldEvaluate<DisContField>>(
            B2, this->particle_group, this->cell_id_translation);

        this->fields["B0"] = B0;
        this->fields["B1"] = B1;
        this->fields["B2"] = B2;
    }

    /**
     * Evaluate E and B at the particle locations.
     */
    inline void evaluate_fields()
    {
        NESOASSERT(this->field_evaluate_E0 != nullptr,
                   "FieldEvaluate not setup.");
        NESOASSERT(this->field_evaluate_E1 != nullptr,
                   "FieldEvaluate not setup.");
        NESOASSERT(this->field_evaluate_E2 != nullptr,
                   "FieldEvaluate not setup.");
        this->field_evaluate_E0->evaluate(Sym<REAL>("E0"));
        this->field_evaluate_E1->evaluate(Sym<REAL>("E1"));
        this->field_evaluate_E2->evaluate(Sym<REAL>("E2"));

        NESOASSERT(this->field_evaluate_B0 != nullptr,
                   "FieldEvaluate not setup.");
        NESOASSERT(this->field_evaluate_B1 != nullptr,
                   "FieldEvaluate not setup.");
        NESOASSERT(this->field_evaluate_B2 != nullptr,
                   "FieldEvaluate not setup.");
        this->field_evaluate_B0->evaluate(Sym<REAL>("B0"));
        this->field_evaluate_B1->evaluate(Sym<REAL>("B1"));
        this->field_evaluate_B2->evaluate(Sym<REAL>("B2"));

        NESOASSERT(this->field_evaluate_ne != nullptr,
                   "FieldEvaluate not setup.");
        this->field_evaluate_ne->evaluate(Sym<REAL>("ELECTRON_DENSITY"));
    }

    inline void remove_marked_particles()
    {
        this->particle_remover->remove(this->particle_group,
                                       (*this->particle_group)[Sym<INT>("ID")],
                                       this->particle_remove_key);
    }

    inline void pre_integration()
    {
        reflection->pre_advection(particle_sub_group(this->particle_group));
    }
    /**
     * Apply boundary conditions.
     */
    inline void boundary_conditions()
    {
        this->reflection->execute(particle_sub_group(this->particle_group));
    }

    /**
     * Adds a velocity drift of \alpha * V_drift where V_drift = E X B
     * evaluation to the velocity of each particle. The coefficient \alpha is
     * the read from the session file key `particle_v_drift_scaling`.
     */
    inline void initialise_particles_from_fields()
    {
        double h_alpha;
        this->config->load_parameter("particle_v_drift_scaling", h_alpha, 1.0);
        const double k_alpha = h_alpha;

        this->evaluate_fields();
        if (this->ndim == 3)
        {
            particle_loop(
                "ParticleSystem:initialise_particles_from_fields",
                this->particle_group,
                [=](auto VELOCITY, auto E0, auto E1, auto E2, auto B0, auto B1,
                    auto B2)
                {
                    // Ei contains d(phi)/dx_i.
                    const auto mE0 = E0.at(0);
                    const auto mE1 = E1.at(0);
                    const auto mE2 = E2.at(0);
                    REAL exb0, exb1, exb2;
                    MAPPING_CROSS_PRODUCT_3D(mE0, mE1, mE2, B0.at(0), B1.at(0),
                                             B2.at(0), exb0, exb1, exb2);
                    VELOCITY.at(0) += k_alpha * exb0;
                    VELOCITY.at(1) += k_alpha * exb1;
                    VELOCITY.at(2) += k_alpha * exb2;
                },
                Access::write(Sym<REAL>("VELOCITY")),
                Access::read(Sym<REAL>("E0")), Access::read(Sym<REAL>("E1")),
                Access::read(Sym<REAL>("E2")), Access::read(Sym<REAL>("B0")),
                Access::read(Sym<REAL>("B1")), Access::read(Sym<REAL>("B2")))
                ->execute();
        }
        else if (this->ndim == 2)
        {
            particle_loop(
                "ParticleSystem:initialise_particles_from_fields",
                this->particle_group,
                [=](auto VELOCITY, auto E0, auto E1, auto B0, auto B1)
                {
                    // Ei contains d(phi)/dx_i.
                    const auto mE0 = E0.at(0);
                    const auto mE1 = E1.at(0);

                    VELOCITY.at(0) += k_alpha * mE0;
                    VELOCITY.at(1) += k_alpha * mE1;
                },
                Access::write(Sym<REAL>("VELOCITY")),
                Access::read(Sym<REAL>("E0")), Access::read(Sym<REAL>("E1")),
                Access::read(Sym<REAL>("B0")), Access::read(Sym<REAL>("B1")))
                ->execute();
        }
    }

protected:
    inline void ionise(const double dt)
    {

        // Evaluate the density and temperature fields at the particle locations
        this->evaluate_fields();

        const double k_dt      = dt;
        const REAL rate        = 100;
        const INT k_remove_key = particle_remove_key;

        auto t0 = profile_timestamp();
        for (auto &[k, v] : this->species_map)
        {
            particle_loop(
                "NeutralParticleSystem::ionise", this->particle_group,
                [=](auto k_ID, auto k_n, auto k_SD, auto k_W)
                {
                    const REAL n      = k_n.at(0);
                    const REAL weight = k_W.at(0);
                    // note that the rate will be a positive number, so minus
                    // sign here
                    REAL deltaweight = -rate * weight * k_dt * n;

                    /* Check whether weight is about to drop below zero.
                       If so, flag particle for removal and adjust deltaweight.
                       These particles are removed after the project call.
                    */
                    if ((weight + deltaweight) <= 0)
                    {
                        k_ID.at(0)  = k_remove_key;
                        deltaweight = -weight;
                    }

                    // Mutate the weight on the particle
                    k_W.at(0) += deltaweight;
                    // Set value for fluid density source (num / Nektar unit
                    // time)
                    k_SD.at(0) = -deltaweight / k_dt;
                },
                Access::write(Sym<INT>("ID")),
                Access::read(Sym<REAL>("ELECTRON_DENSITY")),
                Access::write(Sym<REAL>(v.name + "_SOURCE_DENSITY")),
                Access::write(Sym<REAL>("WEIGHT")))
                ->execute();
        }
        this->sycl_target->profile_map.inc(
            "NeutralParticleSystem", "Ionisation_Execute", 1,
            profile_elapsed(t0, profile_timestamp()));
    }

    virtual inline void integrate_inner(ParticleSubGroupSharedPtr sg,
                                        const double dt_inner)
    {
        const auto k_dt = dt_inner;
        if (this->ndim == 3)
        {
            particle_loop(
                "ParticleSystem:boris", sg,
                [=](auto E0, auto E1, auto E2, auto B0, auto B1, auto B2,
                    auto Q, auto M, auto P, auto V, auto TSP)
                {
                    const REAL dt_left  = k_dt - TSP.at(0);
                    const REAL hdt_left = dt_left * 0.5;
                    if (dt_left > 0.0)
                    {
                        const REAL QoM = Q.at(0) / M.at(0);

                        const REAL scaling_t = QoM * hdt_left;
                        const REAL t_0       = B0.at(0) * scaling_t;
                        const REAL t_1       = B1.at(0) * scaling_t;
                        const REAL t_2       = B2.at(0) * scaling_t;

                        const REAL tmagsq = t_0 * t_0 + t_1 * t_1 + t_2 * t_2;
                        const REAL scaling_s = 2.0 / (1.0 + tmagsq);

                        const REAL s_0 = scaling_s * t_0;
                        const REAL s_1 = scaling_s * t_1;
                        const REAL s_2 = scaling_s * t_2;

                        const REAL V_0 = V.at(0);
                        const REAL V_1 = V.at(1);
                        const REAL V_2 = V.at(2);

                        const REAL v_minus_0 = V_0 + (E0.at(0)) * scaling_t;
                        const REAL v_minus_1 = V_1 + (E1.at(0)) * scaling_t;
                        const REAL v_minus_2 = V_2 + (E2.at(0)) * scaling_t;

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

                        V.at(0) = v_plus_0 + scaling_t * (E0.at(0));
                        V.at(1) = v_plus_1 + scaling_t * (E1.at(0));
                        V.at(2) = v_plus_2 + scaling_t * (E2.at(0));

                        // update of position to next time step
                        P.at(0) += dt_left * V.at(0);
                        P.at(1) += dt_left * V.at(1);
                        P.at(2) += dt_left * V.at(2);
                        TSP.at(0) = k_dt;
                        TSP.at(1) = dt_left;
                    }
                },
                Access::read(Sym<REAL>("E0")), Access::read(Sym<REAL>("E1")),
                Access::read(Sym<REAL>("E2")), Access::read(Sym<REAL>("B0")),
                Access::read(Sym<REAL>("B1")), Access::read(Sym<REAL>("B2")),
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
        // ionise(dt_inner);
    };

    uint64_t total_num_particles_added = 0;
    ;
    const int particle_remove_key = -1;
    std::shared_ptr<ParticleRemover> particle_remover;

    std::map<int, SpeciesInfo> species_map;

    std::vector<Sym<REAL>> src_syms;
    std::vector<int> src_components;

    std::shared_ptr<FieldProject<DisContField>> field_project;

    /// Object used to evaluate Nektar electric field
    std::shared_ptr<FieldEvaluate<DisContField>> field_evaluate_E0;
    std::shared_ptr<FieldEvaluate<DisContField>> field_evaluate_E1;
    std::shared_ptr<FieldEvaluate<DisContField>> field_evaluate_E2;
    /// Object used to evaluate Nektar magnetic field
    std::shared_ptr<FieldEvaluate<DisContField>> field_evaluate_B0;
    std::shared_ptr<FieldEvaluate<DisContField>> field_evaluate_B1;
    std::shared_ptr<FieldEvaluate<DisContField>> field_evaluate_B2;
    /// Object used to project onto Nektar number density field
    std::shared_ptr<FieldEvaluate<DisContField>> field_evaluate_ne;
    std::map<std::string, std::shared_ptr<DisContField>> fields;
    /// Simulation time
    double simulation_time;

    /// Reflective Boundary Conditions
    std::shared_ptr<NektarCompositeTruncatedReflection> reflection;

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
        reflection->pre_advection(sg);
    };

    void apply_boundary_conditions(ParticleSubGroupSharedPtr sg)
    {
        reflection->execute(sg);
    };

    auto find_partial_moves(ParticleSubGroupSharedPtr sg, const double dt)
    {
        return particle_sub_group(
            this->particle_group, [=](auto TSP) { return TSP.at(0) < dt; },
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
