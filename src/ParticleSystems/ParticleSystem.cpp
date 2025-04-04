#include "ParticleSystem.hpp"

namespace NESO::Solvers::tokamak
{

std::string ParticleSystem::class_name =
    GetParticleSystemFactory().RegisterCreatorFunction(
        "ParticleSystem", ParticleSystem::create, "Particle System");

ParticleSystem::ParticleSystem(ParticleReaderSharedPtr session,
                               SD::MeshGraphSharedPtr graph, MPI_Comm comm)
    : PartSysBase(session, graph, comm), simulation_time(0.0) {};

void ParticleSystem::set_up_species()
{
    // get seed from file
    std::srand(std::time(nullptr));
    int seed;

    this->config->load_parameter("particle_position_seed", seed, std::rand());
    double particle_B_scaling;
    this->config->load_parameter("particle_B_scaling", particle_B_scaling, 1.0);

    double particle_thermal_velocity;
    this->config->load_parameter("particle_thermal_velocity",
                                 particle_thermal_velocity, 0.0);
    double particle_velocity_B_scaling;
    this->config->load_parameter("particle_velocity_B_scaling",
                                 particle_velocity_B_scaling, 0.0);

    for (const auto &[k, v] : this->config->get_species())
    {
        std::string name = std::get<0>(v);
        double particle_mass, particle_charge;
        this->config->load_species_parameter(k, "Mass", particle_mass);
        this->config->load_species_parameter(k, "Charge", particle_charge);
        int particle_number = this->config->get_species_initial_N(k);
        long rstart, rend;
        const long size = this->sycl_target->comm_pair.size_parent;
        const long rank = this->sycl_target->comm_pair.rank_parent;
        get_decomp_1d(size, (long)particle_number, rank, &rstart, &rend);
        // const long N = rend - rstart;
        std::mt19937 rng_phasespace(seed + rank);

        // Report param values (later in the initialisation)
        report_param("Random seed", seed);
        report_param("Mass", particle_mass);
        report_param("Charge", particle_charge);
        report_param("Thermal velocity", particle_thermal_velocity);
        report_param("B scaling", particle_B_scaling);
        report_param("Velocity scaling", particle_velocity_B_scaling);
        std::vector<std::vector<double>> positions, velocities;
        std::normal_distribution d{0.0, particle_thermal_velocity};
        std::vector<int> cells;

        auto neqn = this->config->get_species_initial(k, "n");
        rng_phasespace =
            dist_within_extents(graph, neqn, particle_number, positions, cells,
                                1.0e-10, rng_phasespace);
        const long N = cells.size();

        auto veqn  = this->config->get_species_initial(k, "v");
        velocities = NESO::Particles::normal_distribution(
            particle_number, this->ndim, 0.0, particle_thermal_velocity,
            rng_phasespace);

        int id_offset = 0;
        MPICHK(MPI_Exscan(&N, &id_offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD));
        if (N > 0)
        {
            ParticleSet initial_distribution(
                particle_number, this->particle_group->get_particle_spec());

            for (int px = 0; px < particle_number; px++)
            {
                for (int dimx = 0; dimx < this->graph->GetMeshDimension();
                     dimx++)
                {

                    initial_distribution[Sym<REAL>("POSITION")][px][dimx] =
                        positions[dimx][px];
                    initial_distribution[Sym<REAL>("VELOCITY")][px][dimx] =
                        velocities[dimx][px];

                    initial_distribution[Sym<REAL>("ELECTRON_SOURCE_MOMENTUM")]
                                        [px][dimx] = 0.0;
                }
                initial_distribution[Sym<REAL>("Q")][px][0] = particle_charge;
                initial_distribution[Sym<REAL>("M")][px][0] = particle_mass;
                initial_distribution[Sym<INT>("ID")][px][0] = px + id_offset;
                initial_distribution[Sym<INT>("CELL_ID")][px][0] = cells.at(px);
                initial_distribution[Sym<INT>("INTERNAL_STATE")][px][0] = k;
                initial_distribution[Sym<REAL>("WEIGHT")][px][0]        = 1.0;
                initial_distribution[Sym<REAL>("TOT_REACTION_RATE")][px][0] =
                    0.0;
                initial_distribution[Sym<REAL>("ELECTRON_DENSITY")][px][0] =
                    0.0;
                initial_distribution[Sym<REAL>("ELECTRON_SOURCE_ENERGY")][px]
                                    [0] = 0.0;
                initial_distribution[Sym<REAL>("ELECTRON_SOURCE_DENSITY")][px]
                                    [0] = 0.0;
            }

            this->particle_group->add_particles_local(initial_distribution);
        }
        int s = k;
        ParticleSubGroupSharedPtr sub_group =
            std::make_shared<ParticleSubGroup>(
                this->particle_group, [s](auto sid) { return sid[0] == s; },
                Access::read(Sym<INT>("INTERNAL_STATE")));
        species_map[k] =
            SpeciesInfo{name, particle_mass, particle_charge, sub_group};
        this->total_num_particles_added += particle_number;
    }
}

void ParticleSystem::add_sources(double time, double dt)
{
    // get seed from file
    std::srand(std::time(nullptr));
    int seed;

    this->config->load_parameter("particle_position_seed", seed, std::rand());
    double particle_B_scaling;
    this->config->load_parameter("particle_B_scaling", particle_B_scaling, 1.0);

    double particle_thermal_velocity;
    this->config->load_parameter("particle_thermal_velocity",
                                 particle_thermal_velocity, 0.0);
    double particle_velocity_B_scaling;
    this->config->load_parameter("particle_velocity_B_scaling",
                                 particle_velocity_B_scaling, 0.0);

    for (const auto &[k, v] : this->config->get_species())
    {
        std::string name = std::get<0>(v);

        double particle_mass   = species_map[k].mass;
        double particle_charge = species_map[k].charge;

        long rstart, rend;
        const long size = this->sycl_target->comm_pair.size_parent;
        const long rank = this->sycl_target->comm_pair.rank_parent;
        // get_decomp_1d(size, (long)particle_number, rank, &rstart, &rend);
        //  const long N = rend - rstart;
        std::mt19937 rng_phasespace(seed + rank);

        std::vector<std::vector<double>> positions, velocities;
        std::normal_distribution d{0.0, particle_thermal_velocity};

        for (auto &source : this->config->get_species_sources(k))
        {
            int particle_number = source.first;
            std::vector<int> cells;
            LU::EquationSharedPtr neqn = source.second.at(std::pair("n", 0)).m_expression;

            rng_phasespace =
                dist_within_extents(graph, neqn, particle_number, positions,
                                    cells, 1.0e-10, rng_phasespace);
            const long N = cells.size();

            int id_offset = 0;
            MPICHK(MPI_Exscan(&N, &id_offset, 1, MPI_INT, MPI_SUM,
                              MPI_COMM_WORLD));
            if (N > 0)
            {
                ParticleSet src_distribution(
                    particle_number, this->particle_group->get_particle_spec());

                for (int px = 0; px < particle_number; px++)
                {
                    for (int dimx = 0; dimx < this->graph->GetMeshDimension();
                         dimx++)
                    {

                        src_distribution[Sym<REAL>("POSITION")][px][dimx] =
                            positions[dimx][px];
                        src_distribution[Sym<REAL>("VELOCITY")][px][dimx] =
                            velocities[dimx][px];

                        src_distribution[Sym<REAL>("ELECTRON_SOURCE_MOMENTUM")]
                                        [px][dimx] = 0.0;
                    }
                    src_distribution[Sym<REAL>("Q")][px][0] = particle_charge;
                    src_distribution[Sym<REAL>("M")][px][0] = particle_mass;
                    src_distribution[Sym<INT>("ID")][px][0] =
                        px + id_offset + this->total_num_particles_added;
                    src_distribution[Sym<INT>("CELL_ID")][px][0] = cells.at(px);
                    src_distribution[Sym<INT>("INTERNAL_STATE")][px][0] = k;
                    src_distribution[Sym<REAL>("WEIGHT")][px][0]        = 1.0;
                    src_distribution[Sym<REAL>("TOT_REACTION_RATE")][px][0] =
                        0.0;
                    src_distribution[Sym<REAL>("ELECTRON_DENSITY")][px][0] =
                        0.0;
                    src_distribution[Sym<REAL>("ELECTRON_SOURCE_ENERGY")][px]
                                    [0] = 0.0;
                    src_distribution[Sym<REAL>("ELECTRON_SOURCE_DENSITY")][px]
                                    [0] = 0.0;
                }

                this->particle_group->add_particles_local(src_distribution);
            }
            this->total_num_particles_added += particle_number;
        }
        int s = k;
        ParticleSubGroupSharedPtr sub_group =
            std::make_shared<ParticleSubGroup>(
                this->particle_group, [s](auto sid) { return sid[0] == s; },
                Access::read(Sym<INT>("INTERNAL_STATE")));
        species_map[k] =
            SpeciesInfo{name, particle_mass, particle_charge, sub_group};
    }
}

void ParticleSystem::set_up_boundaries()
{
    auto config = std::make_shared<ParameterStore>();
    config->set<REAL>("MapParticlesNewton/newton_tol", 1.0e-10);
    config->set<REAL>("MapParticlesNewton/contained_tol", 1.0e-6);
    config->set<REAL>("CompositeIntersection/newton_tol", 1.0e-10);
    config->set<REAL>("CompositeIntersection/line_intersection_tol", 1.0e-10);
    config->set<REAL>("NektarCompositeTruncatedReflection/reset_distance",
                      1.0e-6);
    auto mesh = std::make_shared<ParticleMeshInterface>(this->graph);
    std::vector<int> reflection_composites;

    for (auto &[k, v] : this->config->get_boundaries())
    {
        for (auto &[sk, sv] : v)
        {
            if (sv == ParticleBoundaryConditionType::eReflective)
            {
                reflection_composites.push_back(k);
            }
        }
    }
    this->reflection = std::make_shared<NektarCompositeTruncatedReflection>(
        Sym<REAL>("VELOCITY"), Sym<REAL>("TSP"), this->sycl_target, mesh,
        reflection_composites, config);
}

} // namespace NESO::Solvers::tokamak