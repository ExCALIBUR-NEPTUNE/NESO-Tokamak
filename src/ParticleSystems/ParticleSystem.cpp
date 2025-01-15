#include "ParticleSystem.hpp"

namespace NESO::Solvers::tokamak
{

std::string ParticleSystem::className =
    GetParticleSystemFactory().RegisterCreatorFunction(
        "ParticleSystem", ParticleSystem::create, "Particle System");

ParticleSystem::ParticleSystem(ParticleReaderSharedPtr session,
                               SD::MeshGraphSharedPtr graph, MPI_Comm comm)
    : PartSysBase(session, graph, comm), simulation_time(0.0) {

          // V initial velocity
          // P uniform sample?
      };

void ParticleSystem::SetUpSpecies()
{
    // get seed from file
    std::srand(std::time(nullptr));
    int seed;
    this->session->LoadParameter("particle_position_seed", seed, std::rand());
    double particle_B_scaling;
    this->session->LoadParameter("particle_B_scaling", particle_B_scaling, 1.0);

    double particle_thermal_velocity;
    this->session->LoadParameter("particle_thermal_velocity",
                                 particle_thermal_velocity, 0.0);
    double particle_velocity_B_scaling;
    this->session->LoadParameter("particle_velocity_B_scaling",
                                 particle_velocity_B_scaling, 0.0);

    for (const auto &[k, v] : this->session->GetSpecies())
    {
        double particle_mass, particle_charge;

        this->session->LoadSpeciesParameter(k, "Mass", particle_mass);
        this->session->LoadSpeciesParameter(k, "Charge", particle_charge);

        long rstart, rend;
        const long size = this->sycl_target->comm_pair.size_parent;
        const long rank = this->sycl_target->comm_pair.rank_parent;
        get_decomp_1d(size, (long)this->num_parts_tot, rank, &rstart, &rend);
        const long N = rend - rstart;

        std::mt19937 rng_phasespace(seed + rank);

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
        ParticleSubGroupSharedPtr sub_group =
            std::make_shared<ParticleSubGroup>(
                this->particle_group, [k](auto sid) { return sid[0] == k; },
                Access::read(Sym<INT>("SPECIES")));
        species_map[k] = SpeciesInfo{std::get<0>(v), particle_mass,
                                     particle_charge, sub_group};

        particle_spec.push(ParticleProp(Sym<REAL>(k + "_src"), 1));
    }
}

void ParticleSystem::SetUpBoundaries()
{
    auto config = std::make_shared<ParameterStore>();
    config->set<REAL>("MapParticlesNewton/newton_tol", 1.0e-10);
    config->set<REAL>("MapParticlesNewton/contained_tol", 1.0e-6);
    config->set<REAL>("CompositeIntersection/newton_tol", 1.0e-10);
    config->set<REAL>("CompositeIntersection/line_intersection_tol", 1.0e-10);
    config->set<REAL>("NektarCompositeTruncatedReflection/reset_distance",
                      1.0e-6);
    auto mesh = std::make_shared<ParticleMeshInterface>(this->graph);
    std::vector<int> reflection_composites = {1, 2, 3, 4};

    for (auto& [k,v] : this->session->GetBoundaries())
    {
        for(auto& [sk, sv] : v)
        {
            if (sv == ParticleBoundaryConditionType::eReflective)
            {
                reflection_composites.push_back
            }
        }
    }
    this->reflection = std::make_shared<NektarCompositeTruncatedReflection>(
        Sym<REAL>("VELOCITY"), Sym<REAL>("TSP"), this->sycl_target, mesh,
        reflection_composites, config);
}

} // namespace NESO::Solvers::tokamak