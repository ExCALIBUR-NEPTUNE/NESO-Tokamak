#include "ParticleSystem.hpp"

namespace NESO::Solvers::tokamak
{

std::string ParticleSystem::class_name =
    GetParticleSystemFactory().RegisterCreatorFunction(
        "ParticleSystem", ParticleSystem::create, "Particle System");

ParticleSystem::ParticleSystem(NESOReaderSharedPtr session,
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
    const long rank = this->sycl_target->comm_pair.rank_parent;
    std::mt19937 rng_phasespace(seed + rank);
    for (const auto &[k, v] : this->config->get_particle_species())
    {
        std::string name = std::get<0>(v);
        double particle_mass, particle_charge;
        this->config->load_particle_species_parameter(k, "Mass", particle_mass);
        this->config->load_particle_species_parameter(k, "Charge",
                                                      particle_charge);
        int particle_number = this->config->get_particle_species_initial_N(k);

        if (particle_number > 0)
        {
            std::vector<std::vector<double>> positions, velocities;
            std::vector<int> cells;
            auto neqn = this->config->get_particle_species_initial(k, "n");
            rng_phasespace =
                dist_within_extents(this->graph, neqn, 0, particle_number,
                                    positions, cells, 1.0e-10, rng_phasespace);
            // auto veqn  = this->config->get_species_initial(k, "v");
            int N         = cells.size();
            int id_offset = 0;
            MPICHK(MPI_Exscan(&N, &id_offset, 1, MPI_INT, MPI_SUM,
                              this->sycl_target->comm));
            if (N > 0)
            {

                velocities = NESO::Particles::normal_distribution(
                    N, this->ndim, 0.0, particle_thermal_velocity,
                    rng_phasespace);
                ParticleSet initial_distribution(
                    N, this->particle_group->get_particle_spec());

                for (int px = 0; px < N; px++)
                {
                    for (int dimx = 0; dimx < this->graph->GetMeshDimension();
                         dimx++)
                    {
                        initial_distribution[Sym<REAL>("POSITION")][px][dimx] =
                            positions[dimx][px];
                        initial_distribution[Sym<REAL>("VELOCITY")][px][dimx] =
                            velocities[dimx][px];

                        initial_distribution[Sym<REAL>(
                            "ELECTRON_SOURCE_MOMENTUM")][px][dimx] = 0.0;
                        initial_distribution[Sym<REAL>("FLUID_FLOW_SPEED")][px]
                                            [dimx] = 0;
                    }
                    initial_distribution[Sym<REAL>("Q")][px][0] =
                        particle_charge;
                    initial_distribution[Sym<REAL>("M")][px][0] = particle_mass;
                    initial_distribution[Sym<INT>("ID")][px][0] =
                        px + id_offset + this->total_num_particles_added;
                    initial_distribution[Sym<INT>("CELL_ID")][px][0] =
                        cells.at(px);
                    initial_distribution[Sym<INT>("INTERNAL_STATE")][px][0] = k;
                    initial_distribution[Sym<REAL>("WEIGHT")][px][0] = 0.02;
                    initial_distribution[Sym<REAL>("TOT_REACTION_RATE")][px]
                                        [0] = 0.0;
                    initial_distribution[Sym<REAL>("ELECTRON_DENSITY")][px][0] =
                        2.0;
                    initial_distribution[Sym<REAL>("ELECTRON_TEMPERATURE")][px]
                                        [0] = 2.0;
                    initial_distribution[Sym<REAL>("ELECTRON_SOURCE_ENERGY")]
                                        [px][0] = 0.0;
                    initial_distribution[Sym<REAL>("ELECTRON_SOURCE_DENSITY")]
                                        [px][0] = 0.0;
                    initial_distribution[Sym<REAL>("FLUID_DENSITY")][px][0] =
                        2.0; // 1e18 m^-3
                    initial_distribution[Sym<REAL>("FLUID_TEMPERATURE")][px]
                                        [0] = 2.0; // eV
                }

                this->particle_group->add_particles_local(initial_distribution);
                this->total_num_particles_added += N;
            }
        }

        int s = k;
        ParticleSubGroupSharedPtr sub_group =
            std::make_shared<ParticleSubGroup>(
                this->particle_group, [s](auto sid) { return sid[0] == s; },
                Access::read(Sym<INT>("INTERNAL_STATE")));
        species_map[k] =
            SpeciesInfo{name, particle_mass, particle_charge, sub_group};
    }
    set_up_boundaries();
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
    const long rank = this->sycl_target->comm_pair.rank_parent;
    std::mt19937 rng_phasespace(seed + rank);

    for (const auto &[k, v] : this->config->get_particle_species())
    {
        std::string name = std::get<0>(v);

        double particle_mass   = species_map[k].mass;
        double particle_charge = species_map[k].charge;

        for (auto &source : this->config->get_particle_species_sources(k))
        {
            long particle_number = source.first;
            if (particle_number > 0)
            {
                std::vector<std::vector<double>> positions, velocities;
                std::vector<int> cells;
                if (auto v = source.second.find(std::pair("n", 0));
                    v != source.second.end()) // Diffuse source
                {
                    rng_phasespace =
                        dist_within_extents(this->graph, v->second.m_expression,
                                            time, particle_number, positions,
                                            cells, 1.0e-10, rng_phasespace);
                }
                else // Point source
                {
                    long rstart, rend;
                    const long size = this->sycl_target->comm_pair.size_parent;
                    const long rank = this->sycl_target->comm_pair.rank_parent;
                    get_decomp_1d(size, particle_number, rank, &rstart, &rend);
                    const long local_particle_number = rend - rstart;

                    double x = source.second.at(std::pair("X", 0))
                                   .m_expression->Evaluate();

                    double y = source.second.at(std::pair("Y", 0))
                                   .m_expression->Evaluate();

                    positions.push_back(
                        std::vector<double>(local_particle_number, x));

                    positions.push_back(
                        std::vector<double>(local_particle_number, y));
                    cells = std::vector<int>(local_particle_number, 0);

                    if (ndim == 3)
                    {
                        double z = source.second.at(std::pair("Z", 0))
                                       .m_expression->Evaluate();
                        positions.push_back(
                            std::vector<double>(local_particle_number, z));
                    }
                }
                // auto veqn  = this->config->get_species_initial(k, "v");
                int N         = cells.size();
                int id_offset = 0;
                MPICHK(MPI_Exscan(&N, &id_offset, 1, MPI_INT, MPI_SUM,
                                  this->sycl_target->comm));
                if (N > 0)
                {

                    velocities = NESO::Particles::normal_distribution(
                        N, this->ndim, 0.0, particle_thermal_velocity,
                        rng_phasespace);
                    ParticleSet src_distribution(
                        N, this->particle_group->get_particle_spec());

                    for (int px = 0; px < N; px++)
                    {
                        for (int dimx = 0;
                             dimx < this->graph->GetMeshDimension(); dimx++)
                        {
                            src_distribution[Sym<REAL>("POSITION")][px][dimx] =
                                positions[dimx][px];
                            src_distribution[Sym<REAL>("VELOCITY")][px][dimx] =
                                velocities[dimx][px];

                            src_distribution[Sym<REAL>(
                                "ELECTRON_SOURCE_MOMENTUM")][px][dimx] = 0.0;
                            src_distribution[Sym<REAL>("FLUID_FLOW_SPEED")][px]
                                            [dimx] = 0;
                        }
                        src_distribution[Sym<REAL>("Q")][px][0] =
                            particle_charge;
                        src_distribution[Sym<REAL>("M")][px][0] = particle_mass;
                        src_distribution[Sym<INT>("ID")][px][0] =
                            px + id_offset + this->total_num_particles_added;
                        src_distribution[Sym<INT>("CELL_ID")][px][0] =
                            cells.at(px);
                        src_distribution[Sym<INT>("INTERNAL_STATE")][px][0] = k;
                        src_distribution[Sym<REAL>("WEIGHT")][px][0] = 0.02;
                        src_distribution[Sym<REAL>("TOT_REACTION_RATE")][px]
                                        [0] = 0.0;
                        src_distribution[Sym<REAL>("ELECTRON_DENSITY")][px][0] =
                            2.0;
                        src_distribution[Sym<REAL>("ELECTRON_TEMPERATURE")][px]
                                        [0] = 2.0;
                        src_distribution[Sym<REAL>("ELECTRON_SOURCE_ENERGY")]
                                        [px][0] = 0.0;
                        src_distribution[Sym<REAL>("ELECTRON_SOURCE_DENSITY")]
                                        [px][0] = 0.0;
                        src_distribution[Sym<REAL>("FLUID_DENSITY")][px][0] =
                            2.0; // 1e18 m^-3
                        src_distribution[Sym<REAL>("FLUID_TEMPERATURE")][px]
                                        [0] = 2.0; // eV
                    }

                    this->particle_group->add_particles_local(src_distribution);
                    this->total_num_particles_added += N;
                }
            }
        }
        int s = k;
        ParticleSubGroupSharedPtr sub_group =
            std::make_shared<ParticleSubGroup>(
                this->particle_group, [s](auto sid) { return sid[0] == s; },
                Access::read(Sym<INT>("INTERNAL_STATE")));
        species_map[k] =
            SpeciesInfo{name, particle_mass, particle_charge, sub_group};
    }
    transfer_particles();
}

void ParticleSystem::add_sinks(double time, double dt)
{
    // get seed from file
    std::srand(std::time(nullptr));
    int seed;

    const long rank = this->sycl_target->comm_pair.rank_parent;
    std::mt19937 rng_phasespace(seed + rank);

    std::uniform_real_distribution<> rng_dist(0, 1);

    Array<OneD, Array<OneD, NekDouble>> posarr(3);
    for (int dimx = 0; dimx < 3; dimx++)
    {
        posarr[dimx] = Array<OneD, NekDouble>(
            this->particle_group->get_npart_local(), 0.0);
    }

    // for each cell
    const int cell_count = this->domain->mesh->get_cell_count();
    for (int particle_index = 0, cellx = 0; cellx < cell_count; cellx++)
    {
        // for each particle in the cell

        auto cells =
            (*particle_group)[Sym<INT>("CELL_ID")]->cell_dat.get_cell(cellx);
        const int nrow = cells->nrow;
        // many cells will be empty so we check before issuing more
        // copies
        if (nrow > 0)
        {
            auto phys_positions =
                (*particle_group)[Sym<REAL>("POSITION")]->cell_dat.get_cell(
                    cellx);

            for (int rowx = 0; rowx < nrow; rowx++)
            {
                // copy the particle data into the store of points
                for (int dimx = 0; dimx < ndim; dimx++)
                {
                    posarr[dimx][particle_index] =
                        (*phys_positions)[dimx][rowx];
                }
                particle_index++;
            }
        }
    }
    for (const auto &[k, v] : this->config->get_particle_species())
    {
        for (auto &sink : this->config->get_particle_species_sinks(k))
        {
            LU::EquationSharedPtr neqn =
                sink.at(std::pair("n", 0)).m_expression;
            Array<OneD, NekDouble> eqnarr(
                this->particle_group->get_npart_local());

            neqn->Evaluate(posarr[0], posarr[1], posarr[2], time, eqnarr);

            for (int particle_index = 0, cellx = 0; cellx < cell_count; cellx++)
            {
                auto cells =
                    (*particle_group)[Sym<INT>("CELL_ID")]->cell_dat.get_cell(
                        cellx);
                const int nrow = cells->nrow;
                if (nrow > 0)
                {
                    auto id =
                        (*particle_group)[Sym<INT>("ID")]->cell_dat.get_cell(
                            cellx);
                    auto is = (*particle_group)[Sym<INT>("INTERNAL_STATE")]
                                  ->cell_dat.get_cell(cellx);
                    for (int rowx = 0; rowx < nrow; rowx++)
                    {
                        if (eqnarr[particle_index] > rng_dist(rng_phasespace) &&
                            k == (*is)[0][rowx])
                        {
                            (*id)[0][rowx] = particle_remove_key;
                        }
                        particle_index++;
                    }
                    (*particle_group)[Sym<INT>("ID")]->cell_dat.set_cell(cellx,
                                                                         id);
                }
            }
        }
    }
    remove_marked_particles();
}

void ParticleSystem::set_up_boundaries()
{
    auto store = std::make_shared<ParameterStore>();
    store->set<REAL>("MapParticlesNewton/newton_tol", 1.0e-10);
    store->set<REAL>("MapParticlesNewton/contained_tol", 1.0e-6);
    store->set<REAL>("CompositeIntersection/newton_tol", 1.0e-10);
    store->set<REAL>("CompositeIntersection/line_intersection_tol", 1.0e-10);
    store->set<REAL>("NektarCompositeTruncatedReflection/reset_distance",
                     1.0e-6);
    auto mesh = std::make_shared<ParticleMeshInterface>(this->graph);

    std::vector<int> reflection_composites;

    for (auto &[sk, sv] : this->config->get_particle_species_boundary(0))
    {
        if (sv == ParticleBoundaryConditionType::eReflective)
        {
            reflection_composites.push_back(sk);
        }
    }
    this->reflection = std::make_shared<NektarCompositeTruncatedReflection>(
        Sym<REAL>("VELOCITY"), Sym<REAL>("TSP"), this->sycl_target, mesh,
        reflection_composites, store);

    // for (auto &[k, v] : this->get_species())
    // {
    //     std::cout << "Setup\n";

    //     for (auto &[sk, sv] : this->config->get_particle_species_boundary(k))
    //     {
    //         if (sv == ParticleBoundaryConditionType::eReflective)
    //         {
    //             std::cout << "Composites\n";

    //             v.reflection_composites.push_back(sk);
    //         }
    //     }
    //     if (!v.reflection_composites.empty())
    //     {
    //         std::cout << "Reflection\n";

    //         v.reflection =
    //         std::make_shared<NektarCompositeTruncatedReflection>(
    //             Sym<REAL>("VELOCITY"), Sym<REAL>("TSP"), this->sycl_target,
    //             mesh, v.reflection_composites, store);
    //     }
    // }
}

} // namespace NESO::Solvers::tokamak