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

#include <LibUtilities/BasicUtils/SessionReader.h>

// namespace LU = Nektar::LibUtilities;
// namespace NP = NESO::Particles;

namespace NESO::Solvers::tokamak
{

/**
 * @brief
 */
class ParticleSystem : public PartSysBase
{

    static inline ParticleSpec particle_spec{
        ParticleProp(Sym<INT>("CELL_ID"), 1, true),
        ParticleProp(Sym<INT>("PARTICLE_ID"), 1),
        ParticleProp(Sym<REAL>("POSITION"), 2, true),
        ParticleProp(Sym<REAL>("VELOCITY"), 2),
        ParticleProp(Sym<REAL>("M"), 1),
        ParticleProp(Sym<REAL>("Q"), 1),
        ParticleProp(Sym<REAL>("E0"), 1),
        ParticleProp(Sym<REAL>("E1"), 1),
        ParticleProp(Sym<REAL>("E2"), 1),
        ParticleProp(Sym<REAL>("B0"), 1),
        ParticleProp(Sym<REAL>("B1"), 1),
        ParticleProp(Sym<REAL>("B2"), 1)};

public:
    /**
     *  Create a new instance.
     *
     *  @param session Nektar++ session to use for parameters and simulation
     * specification.
     *  @param graph Nektar++ MeshGraph on which particles exist.
     *  @param comm (optional) MPI communicator to use - default MPI_COMM_WORLD.
     *
     */
    ParticleSystem(LU::SessionReaderSharedPtr session,
                   SD::MeshGraphSharedPtr graph, MPI_Comm comm = MPI_COMM_WORLD)
        : PartSysBase(session, graph, particle_spec, comm), simulation_time(0.0)
    {
        // this->pbc = std::make_shared<NektarCartesianPeriodic>(
        //     this->sycl_target, this->graph,
        //     this->particle_group->position_dat);
        auto config = std::make_shared<ParameterStore>();
        config->set<REAL>("MapParticlesNewton/newton_tol", 1.0e-10);
        config->set<REAL>("MapParticlesNewton/contained_tol", 1.0e-6);
        config->set<REAL>("CompositeIntersection/newton_tol", 1.0e-10);
        config->set<REAL>("CompositeIntersection/line_intersection_tol",
                          1.0e-10);
        config->set<REAL>("NektarCompositeTruncatedReflection/reset_distance",
                          1.0e-6);
        auto mesh = std::make_shared<ParticleMeshInterface>(this->graph);
        std::vector<int> reflection_composites = {100};
        this->reflection = std::make_shared<NektarCompositeTruncatedReflection>(
            Sym<REAL>("VELOCITY"), Sym<REAL>("TSP"), this->sycl_target, mesh,
            reflection_composites, config);

        // V initial velocity
        // P uniform sample?

        long rstart, rend;
        const long size = this->sycl_target->comm_pair.size_parent;
        const long rank = this->sycl_target->comm_pair.rank_parent;
        get_decomp_1d(size, (long)this->num_parts_tot, rank, &rstart, &rend);
        const long N = rend - rstart;

        // get seed from file
        std::srand(std::time(nullptr));
        int seed;
        this->session->LoadParameter("particle_position_seed", seed,
                                     std::rand());
        std::mt19937 rng_phasespace(seed + rank);

        // get particle mass from file
        double particle_mass;
        this->session->LoadParameter("particle_mass", particle_mass, 1.0);
        // get particle charge from file
        double particle_charge;
        this->session->LoadParameter("particle_charge", particle_charge, 1.0);

        double particle_B_scaling;
        this->session->LoadParameter("particle_B_scaling", particle_B_scaling,
                                     1.0);

        double particle_thermal_velocity;
        this->session->LoadParameter("particle_thermal_velocity",
                                     particle_thermal_velocity, 0.0);
        double particle_velocity_B_scaling;
        this->session->LoadParameter("particle_velocity_B_scaling",
                                     particle_velocity_B_scaling, 0.0);

        this->particle_remover =
            std::make_shared<ParticleRemover>(this->sycl_target);
        // Report param values (later in the initialisation)
        report_param("Random seed", seed);
        report_param("Mass", particle_mass);
        report_param("Charge", particle_charge);
        report_param("Thermal velocity", particle_thermal_velocity);
        report_param("B scaling", particle_B_scaling);
        report_param("Velocity scaling", particle_velocity_B_scaling);

        std::normal_distribution d{0.0, particle_thermal_velocity};

        if (N > 0)
        {
            ParticleSet initial_distribution(
                N, this->particle_group->get_particle_spec());
            // auto positions = uniform_within_extents(
            //     N, ndim, this->pbc->global_extent, rng_phasespace);
            const int npart_per_cell = 10;
            std::mt19937 rng(534234 + rank);
            std::vector<std::vector<double>> positions;
            std::vector<int> cells;
            rng = uniform_within_elements(graph, npart_per_cell, positions,
                                          cells, 1.0e-10, rng);

            for (int px = 0; px < N; px++)
            {
                for (int dimx = 0; dimx < this->graph->GetMeshDimension();
                     dimx++)
                {
                    // Create the position in dimx
                    // const REAL pos_shift = this->pbc->global_origin[dimx];
                    const REAL pos_tmp = /*pos_shift + */ positions[dimx][px];
                    initial_distribution[Sym<REAL>("POSITION")][px][dimx] =
                        pos_tmp;
                    // Create the velocity in dimx
                    REAL vtmp = 0.0;

                    if (particle_thermal_velocity > 0.0)
                    {
                        vtmp += d(rng_phasespace);
                    }
                    initial_distribution[Sym<REAL>("VELOCITY")][px][dimx] =
                        vtmp;
                }
                initial_distribution[Sym<REAL>("Q")][px][0] = particle_charge;
                initial_distribution[Sym<REAL>("M")][px][0] = particle_mass;
                initial_distribution[Sym<INT>("PARTICLE_ID")][px][0] =
                    px + rstart;
            }

            this->particle_group->add_particles_local(initial_distribution);
        }

        parallel_advection_initialisation(this->particle_group);
        parallel_advection_store(this->particle_group);
        const int num_steps = 20;
        for (int stepx = 0; stepx < num_steps; stepx++)
        {
            parallel_advection_step(this->particle_group, num_steps, stepx);
            this->transfer_particles();
        }
        parallel_advection_restore(this->particle_group);
        // Move particles to the owning ranks and correct cells.
        this->transfer_particles();

        init_output("particle_trajectory.h5part", Sym<INT>("CELL_ID"),
                    Sym<REAL>("VELOCITY"), Sym<REAL>("E0"), Sym<REAL>("E1"),
                    Sym<REAL>("E2"), Sym<INT>("PARTICLE_ID"));
    };

    /// Disable (implicit) copies.
    ParticleSystem(const ParticleSystem &st) = delete;
    /// Disable (implicit) copies.
    ParticleSystem &operator=(ParticleSystem const &a) = delete;

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
            this->integrate_inner(dt_inner);
            time_tmp += dt_inner;
        }

        this->simulation_time = time_end;
        this->transfer_particles();
    }

    /**
     * Setup the projection object to use the following fields.
     *
     * @param rho_src Nektar++ fields to project ionised particle data onto.
     */
    inline void setup_project(std::shared_ptr<DisContField> ni_src,
                              std::shared_ptr<DisContField> E_src)
    {
        std::vector<std::shared_ptr<DisContField>> fields = {ni_src, E_src};
        this->field_project = std::make_shared<FieldProject<DisContField>>(
            fields, this->particle_group, this->cell_id_translation);

        // Setup debugging output for each field
        this->fields["ni_src"] = ni_src;
        this->fields["E_src"]  = E_src;
    }

    /**
     *  Project the plasma source and momentum contributions from particle data
     *  onto field data.
     */
    inline void project_source_terms()
    {
        NESOASSERT(this->field_project != nullptr,
                   "Field project object is null. Was setup_project called?");

        std::vector<Sym<REAL>> syms = {
            Sym<REAL>("SOURCE_DENSITY"), Sym<REAL>("SOURCE_MOMENTUM"),
            Sym<REAL>("SOURCE_MOMENTUM"), Sym<REAL>("SOURCE_ENERGY")};
        std::vector<int> components = {0, 0, 1, 0};
        this->field_project->project(syms, components);

        // remove fully ionised particles from the simulation
        remove_marked_particles();
    }

    /**
     * Set up the evaluation of grad(phi).
     *
     * @param gradphi0 Nektar++ field storing grad(phi) in direction 0.
     * @param gradphi1 Nektar++ field storing grad(phi) in direction 1.
     * @param gradphi2 Nektar++ field storing grad(phi) in direction 2.
     */
    inline void setup_evaluate_grad_phi(std::shared_ptr<DisContField> gradphi0,
                                        std::shared_ptr<DisContField> gradphi1,
                                        std::shared_ptr<DisContField> gradphi2)
    {
        // TODO redo with redone Bary Evaluate.
        this->field_evaluate_gradphi0 =
            std::make_shared<FieldEvaluate<DisContField>>(
                gradphi0, this->particle_group, this->cell_id_translation);
        this->field_evaluate_gradphi1 =
            std::make_shared<FieldEvaluate<DisContField>>(
                gradphi1, this->particle_group, this->cell_id_translation);
        this->field_evaluate_gradphi2 =
            std::make_shared<FieldEvaluate<DisContField>>(
                gradphi2, this->particle_group, this->cell_id_translation);

        this->fields["gradphi0"] = gradphi0;
        this->fields["gradphi1"] = gradphi1;
        this->fields["gradphi2"] = gradphi2;
    }

    /**
     * Set up the evaluation of B.
     *
     * @param gradphi0 Nektar++ field storing B in direction 0.
     * @param gradphi1 Nektar++ field storing B in direction 1.
     * @param gradphi2 Nektar++ field storing B in direction 2.
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
     * Evaluate grad(phi) and B at the particle locations.
     */
    inline void evaluate_fields()
    {
        NESOASSERT(this->field_evaluate_gradphi0 != nullptr,
                   "FieldEvaluate not setup.");
        NESOASSERT(this->field_evaluate_gradphi1 != nullptr,
                   "FieldEvaluate not setup.");
        NESOASSERT(this->field_evaluate_gradphi2 != nullptr,
                   "FieldEvaluate not setup.");
        this->field_evaluate_gradphi0->evaluate(Sym<REAL>("E0"));
        this->field_evaluate_gradphi1->evaluate(Sym<REAL>("E1"));
        this->field_evaluate_gradphi2->evaluate(Sym<REAL>("E2"));

        NESOASSERT(this->field_evaluate_B0 != nullptr,
                   "FieldEvaluate not setup.");
        NESOASSERT(this->field_evaluate_B1 != nullptr,
                   "FieldEvaluate not setup.");
        NESOASSERT(this->field_evaluate_B2 != nullptr,
                   "FieldEvaluate not setup.");
        this->field_evaluate_B0->evaluate(Sym<REAL>("B0"));
        this->field_evaluate_B1->evaluate(Sym<REAL>("B1"));
        this->field_evaluate_B2->evaluate(Sym<REAL>("B2"));
    }

    inline void remove_marked_particles()
    {
        this->particle_remover->remove(
            this->particle_group,
            (*this->particle_group)[Sym<INT>("PARTICLE_ID")],
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
        // this->pbc->execute();
        this->reflection->execute(particle_sub_group(this->particle_group));
    }

    /**
     * Adds a velocity drift of \alpha * V_drift where V_drift = -grad(phi) X B
     * evaluation to the velocity of each particle. The coefficient \alpha is
     * the read from the session file key `particle_v_drift_scaling`.
     */
    inline void initialise_particles_from_fields()
    {
        double h_alpha;
        this->session->LoadParameter("particle_v_drift_scaling", h_alpha, 1.0);
        const double k_alpha = h_alpha;

        this->evaluate_fields();
        particle_loop(
            "ParticleSystem:initialise_particles_from_fields",
            this->particle_group,
            [=](auto VELOCITY, auto E0, auto E1, auto E2, auto B0, auto B1,
                auto B2)
            {
                // Ei contains d(phi)/dx_i.
                const auto mE0 = -1.0 * E0.at(0);
                const auto mE1 = -1.0 * E1.at(0);
                const auto mE2 = -1.0 * E2.at(0);
                REAL exb0, exb1, exb2;
                MAPPING_CROSS_PRODUCT_3D(mE0, mE1, mE2, B0.at(0), B1.at(0),
                                         B2.at(0), exb0, exb1, exb2);
                VELOCITY.at(0) += k_alpha * exb0;
                VELOCITY.at(1) += k_alpha * exb1;
                VELOCITY.at(2) += k_alpha * exb2;
            },
            Access::write(Sym<REAL>("VELOCITY")), Access::read(Sym<REAL>("E0")),
            Access::read(Sym<REAL>("E1")), Access::read(Sym<REAL>("E2")),
            Access::read(Sym<REAL>("B0")), Access::read(Sym<REAL>("B1")),
            Access::read(Sym<REAL>("B2")))
            ->execute();
    }

protected:
    inline void integrate_inner(const double dt_inner)
    {
        const auto k_dt  = dt_inner;
        const auto k_dht = dt_inner * 0.5;

        particle_loop(
            "ParticleSystem:boris", this->particle_group,
            [=](auto E0, auto E1, auto E2, auto B0, auto B1, auto B2, auto Q,
                auto M, auto P, auto V)
            {
                const REAL QoM = Q.at(0) / M.at(0);

                const REAL scaling_t = QoM * k_dht;
                const REAL t_0       = B0.at(0) * scaling_t;
                const REAL t_1       = B1.at(0) * scaling_t;
                const REAL t_2       = B2.at(0) * scaling_t;

                const REAL tmagsq    = t_0 * t_0 + t_1 * t_1 + t_2 * t_2;
                const REAL scaling_s = 2.0 / (1.0 + tmagsq);

                const REAL s_0 = scaling_s * t_0;
                const REAL s_1 = scaling_s * t_1;
                const REAL s_2 = scaling_s * t_2;

                const REAL V_0 = V.at(0);
                const REAL V_1 = V.at(1);
                const REAL V_2 = V.at(2);

                // The E dat contains d(phi)/dx not E -> multiply by -1.
                const REAL v_minus_0 = V_0 + (-1.0 * E0.at(0)) * scaling_t;
                const REAL v_minus_1 = V_1 + (-1.0 * E1.at(0)) * scaling_t;
                const REAL v_minus_2 = V_2 + (-1.0 * E2.at(0)) * scaling_t;

                REAL v_prime_0, v_prime_1, v_prime_2;
                MAPPING_CROSS_PRODUCT_3D(v_minus_0, v_minus_1, v_minus_2, t_0,
                                         t_1, t_2, v_prime_0, v_prime_1,
                                         v_prime_2)

                v_prime_0 += v_minus_0;
                v_prime_1 += v_minus_1;
                v_prime_2 += v_minus_2;

                REAL v_plus_0, v_plus_1, v_plus_2;
                MAPPING_CROSS_PRODUCT_3D(v_prime_0, v_prime_1, v_prime_2, s_0,
                                         s_1, s_2, v_plus_0, v_plus_1, v_plus_2)

                v_plus_0 += v_minus_0;
                v_plus_1 += v_minus_1;
                v_plus_2 += v_minus_2;

                // The E dat contains d(phi)/dx not E -> multiply by -1.
                V.at(0) = v_plus_0 + scaling_t * (-1.0 * E0.at(0));
                V.at(1) = v_plus_1 + scaling_t * (-1.0 * E1.at(0));
                V.at(2) = v_plus_2 + scaling_t * (-1.0 * E2.at(0));

                // update of position to next time step
                P.at(0) += k_dt * V.at(0);
                P.at(1) += k_dt * V.at(1);
                P.at(2) += k_dt * V.at(2);
            },
            Access::read(Sym<REAL>("E0")), Access::read(Sym<REAL>("E1")),
            Access::read(Sym<REAL>("E2")), Access::read(Sym<REAL>("B0")),
            Access::read(Sym<REAL>("B1")), Access::read(Sym<REAL>("B2")),
            Access::read(Sym<REAL>("Q")), Access::read(Sym<REAL>("M")),
            Access::write(Sym<REAL>("POSITION")),
            Access::write(Sym<REAL>("VELOCITY")))
            ->execute();
    };
    const int particle_remove_key = -1;
    std::shared_ptr<ParticleRemover> particle_remover;

    std::shared_ptr<FieldProject<DisContField>> field_project;
    /// Object used to evaluate Nektar electric field
    std::shared_ptr<FieldEvaluate<DisContField>> field_evaluate_gradphi0;
    std::shared_ptr<FieldEvaluate<DisContField>> field_evaluate_gradphi1;
    std::shared_ptr<FieldEvaluate<DisContField>> field_evaluate_gradphi2;
    /// Object used to evaluate Nektar magnetic field
    std::shared_ptr<FieldEvaluate<DisContField>> field_evaluate_B0;
    std::shared_ptr<FieldEvaluate<DisContField>> field_evaluate_B1;
    std::shared_ptr<FieldEvaluate<DisContField>> field_evaluate_B2;
    /// Object used to project onto Nektar number density field
    std::map<std::string, std::shared_ptr<DisContField>> fields;
    /// Simulation time
    double simulation_time;

    /// Periodic Boundary Conditions
    // std::shared_ptr<NektarCartesianPeriodic> pbc;
    /// Reflective Boundary Conditions
    std::shared_ptr<NektarCompositeTruncatedReflection> reflection;

    /**
     *  Apply boundary conditions and transfer particles between MPI ranks.
     * // Move some of this to PartSysBase / make it a pure-virtual func?
     */
    inline void transfer_particles()
    {
        auto t0 = profile_timestamp();
        this->pre_integration();
        this->boundary_conditions();
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
