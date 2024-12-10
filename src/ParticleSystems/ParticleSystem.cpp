#include "ParticleSystem.hpp"

namespace NESO::Solvers::tokamak
{
ParticleSystemFactory &GetParticleSystemFactory()
{
    static ParticleSystemFactory instance;
    return instance;
}

void ParticleSystem::InitSpec()
{
    particle_spec = {ParticleProp(Sym<INT>("CELL_ID"), 1, true),
                     ParticleProp(Sym<INT>("PARTICLE_ID"), 1),
                     ParticleProp(Sym<REAL>("POSITION"), 3, true),
                     ParticleProp(Sym<REAL>("VELOCITY"), 3),
                     ParticleProp(Sym<REAL>("M"), 1),
                     ParticleProp(Sym<REAL>("Q"), 1),
                     ParticleProp(Sym<REAL>("SOURCE_DENSITY"), 1),
                     ParticleProp(Sym<REAL>("SOURCE_ENERGY"), 1),
                     ParticleProp(Sym<REAL>("E0"), 1),
                     ParticleProp(Sym<REAL>("E1"), 1),
                     ParticleProp(Sym<REAL>("E2"), 1),
                     ParticleProp(Sym<REAL>("B0"), 1),
                     ParticleProp(Sym<REAL>("B1"), 1),
                     ParticleProp(Sym<REAL>("B2"), 1)};
}

ParticleSystem::ParticleSystem(LU::SessionReaderSharedPtr session,
                               SD::MeshGraphSharedPtr graph, MPI_Comm comm)
    : PartSysBase(session, graph, particle_spec, comm), simulation_time(0.0)
{
    InitSpec();
    auto config = std::make_shared<ParameterStore>();
    config->set<REAL>("MapParticlesNewton/newton_tol", 1.0e-10);
    config->set<REAL>("MapParticlesNewton/contained_tol", 1.0e-6);
    config->set<REAL>("CompositeIntersection/newton_tol", 1.0e-10);
    config->set<REAL>("CompositeIntersection/line_intersection_tol", 1.0e-10);
    config->set<REAL>("NektarCompositeTruncatedReflection/reset_distance",
                      1.0e-6);
    auto mesh = std::make_shared<ParticleMeshInterface>(this->graph);
    std::vector<int> reflection_composites = {1, 2, 3, 4};
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
    this->session->LoadParameter("particle_position_seed", seed, std::rand());
    std::mt19937 rng_phasespace(seed + rank);

    // get particle mass from file
    double particle_mass;
    this->session->LoadParameter("particle_mass", particle_mass, 1.0);
    // get particle charge from file
    double particle_charge;
    this->session->LoadParameter("particle_charge", particle_charge, 1.0);

    double particle_B_scaling;
    this->session->LoadParameter("particle_B_scaling", particle_B_scaling, 1.0);

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

        const int npart_per_cell = 10;
        std::mt19937 rng(534234 + rank);
        std::vector<std::vector<double>> positions;
        std::vector<int> cells;
        rng = uniform_within_elements(graph, npart_per_cell, positions, cells,
                                      1.0e-10, rng);
        for (int px = 0; px < N; px++)
        {
            for (int dimx = 0; dimx < this->graph->GetMeshDimension(); dimx++)
            {
                // Create the position in dimx
                // const REAL pos_shift = this->pbc->global_origin[dimx];
                const REAL pos_tmp = /*pos_shift + */ positions[dimx][px];
                initial_distribution[Sym<REAL>("POSITION")][px][dimx] = pos_tmp;
                // Create the velocity in dimx
                REAL vtmp = 0.0;

                if (particle_thermal_velocity > 0.0)
                {
                    vtmp += d(rng_phasespace);
                }
                initial_distribution[Sym<REAL>("VELOCITY")][px][dimx] = vtmp;
            }
            initial_distribution[Sym<REAL>("Q")][px][0] = particle_charge;
            initial_distribution[Sym<REAL>("M")][px][0] = particle_mass;
            initial_distribution[Sym<INT>("PARTICLE_ID")][px][0] = px + rstart;
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

inline void ParticleSystem::integrate(const double time_end, const double dt)
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

inline void ParticleSystem::setup_evaluate_E(std::shared_ptr<DisContField> E0,
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

inline void ParticleSystem::setup_evaluate_B(std::shared_ptr<DisContField> B0,
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

inline void ParticleSystem::evaluate_fields()
{
    NESOASSERT(this->field_evaluate_E0 != nullptr, "FieldEvaluate not setup.");
    NESOASSERT(this->field_evaluate_E1 != nullptr, "FieldEvaluate not setup.");
    NESOASSERT(this->field_evaluate_E2 != nullptr, "FieldEvaluate not setup.");
    this->field_evaluate_E0->evaluate(Sym<REAL>("E0"));
    this->field_evaluate_E1->evaluate(Sym<REAL>("E1"));
    this->field_evaluate_E2->evaluate(Sym<REAL>("E2"));

    NESOASSERT(this->field_evaluate_B0 != nullptr, "FieldEvaluate not setup.");
    NESOASSERT(this->field_evaluate_B1 != nullptr, "FieldEvaluate not setup.");
    NESOASSERT(this->field_evaluate_B2 != nullptr, "FieldEvaluate not setup.");
    this->field_evaluate_B0->evaluate(Sym<REAL>("B0"));
    this->field_evaluate_B1->evaluate(Sym<REAL>("B1"));
    this->field_evaluate_B2->evaluate(Sym<REAL>("B2"));
}
} // namespace NESO::Solvers::tokamak