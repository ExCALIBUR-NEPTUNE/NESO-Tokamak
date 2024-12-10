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

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>

// namespace LU = Nektar::LibUtilities;
// namespace NP = NESO::Particles;

namespace NESO::Solvers::tokamak
{
class ParticleSystem;

typedef ParticleSystemSharedPtr std::shared_ptr<ParticleSystem>;
typedef ParticleSystemFactory LU::NekFactory<std::string, ParticleSystem,
                                             const LU::SessionReaderSharedPtr,
                                             const SD::MeshGraphSharedPtr>;
ParticleSystemFactory &GetParticleSystemFactory();
/**
 * @brief
 */
class ParticleSystem : public PartSysBase
{

public:
    /**
     * @brief Create an instance of this class and initialise it.
     */
    static ParticleSystemSharedPtr create(
        const LU::SessionReaderSharedPtr &session,
        const SD::MeshGraphSharedPtr &graph)
    {
        ParticleSystemSharedPtr p =
            MemoryManager<ParticleSystem>::AllocateSharedPtr(session, graph);
        return p;
    }

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
                   SD::MeshGraphSharedPtr graph,
                   MPI_Comm comm = MPI_COMM_WORLD);

    virtual ~ParticleSystem() = 0;

    /// Disable (implicit) copies.
    ParticleSystem(const ParticleSystem &st) = delete;
    /// Disable (implicit) copies.
    ParticleSystem &operator=(ParticleSystem const &a) = delete;

    virtual void InitSpec();

    /**
     *  Integrate the particle system forward to the requested time using
     *  (at most) the requested time step.
     *
     *  @param time_end Target time to integrate to.
     *  @param dt Time step size.
     */
    inline void integrate(const double time_end, const double dt);

    /**
     * Setup the projection object to use the following fields.
     *
     * @param rho_src Nektar++ fields to project ionised particle data onto.
     */
    inline void setup_project(
        std::vector<std::shared_ptr<DisContField>> &src_fields)
    {

        this->field_project = std::make_shared<FieldProject<DisContField>>(
            src_fields, this->particle_group, this->cell_id_translation);
    }

    /**
     *  Project the plasma source and momentum contributions from particle data
     *  onto field data.
     */
    inline void project_source_terms()
    {
        NESOASSERT(this->field_project != nullptr,
                   "Field project object is null. Was setup_project called?");

        std::vector<Sym<REAL>> syms = {Sym<REAL>("SOURCE_DENSITY"),
                                       Sym<REAL>("SOURCE_ENERGY")};
        std::vector<int> components = {0, 0};

        for (Species &s : species_list)
        {
            // Uncomment when sub_group projection enabled
            // this->field_project->project(s.sub_group, vsyms, components);
        }
        // remove fully ionised particles from the simulation
        remove_marked_particles();
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
                                 std::shared_ptr<DisContField> E2);

    /**
     * Set up the evaluation of B.
     *
     * @param B0 Nektar++ field storing B in direction 0.
     * @param B1 Nektar++ field storing B in direction 1.
     * @param B2 Nektar++ field storing B in direction 2.
     */
    inline void setup_evaluate_B(std::shared_ptr<DisContField> B0,
                                 std::shared_ptr<DisContField> B1,
                                 std::shared_ptr<DisContField> B2);

    /**
     * Evaluate E and B at the particle locations.
     */
    inline void evaluate_fields();

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
            Access::write(Sym<REAL>("VELOCITY")), Access::read(Sym<REAL>("E0")),
            Access::read(Sym<REAL>("E1")), Access::read(Sym<REAL>("E2")),
            Access::read(Sym<REAL>("B0")), Access::read(Sym<REAL>("B1")),
            Access::read(Sym<REAL>("B2")))
            ->execute();
    }

protected:
    virtual inline void integrate_inner(const double dt_inner)
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

                const REAL v_minus_0 = V_0 + (E0.at(0)) * scaling_t;
                const REAL v_minus_1 = V_1 + (E1.at(0)) * scaling_t;
                const REAL v_minus_2 = V_2 + (E2.at(0)) * scaling_t;

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

                V.at(0) = v_plus_0 + scaling_t * (E0.at(0));
                V.at(1) = v_plus_1 + scaling_t * (E1.at(0));
                V.at(2) = v_plus_2 + scaling_t * (E2.at(0));

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
    ParticleSpec particle_spec;

    const int particle_remove_key = -1;
    std::shared_ptr<ParticleRemover> particle_remover;

    struct Species
    {
        std::string name;
        double mass;
        double charge;

        std::shared_ptr<ParticleSubGroup> sub_group;
    };

    std::vector<Species> species_list;

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
    std::map<std::string, std::shared_ptr<DisContField>> fields;
    /// Simulation time
    double simulation_time;

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

    // To be moved into NESO

    virtual void ReadParticles()
    {
        // Check we actually have a document loaded.
        ASSERTL0(&session->GetDocument(), "No XML document loaded.");

        TiXmlHandle docHandle(&session->GetDocument());
        TiXmlElement *particles;

        // Look for all data in PARTICLES block.
        particles = docHandle.FirstChildElement("NEKTAR")
                        .FirstChildElement("PARTICLES")
                        .Element();

        if (!particles)
        {
            return;
        }

        TiXmlElement *species = particles->FirstChildElement("SPECIES");
        TiXmlElement *specie  = species->FirstChildElement("S");

        while (specie)
        {
            nSpecies++;
            std::stringstream tagcontent;
            tagcontent << *specie;

            ASSERTL0(specie->Attribute("ID"),
                     "Missing ID attribute in Species XML "
                     "element: \n\t'" +
                         tagcontent.str() + "'");
            std::string name = species->Attribute("NAME");
            ASSERTL0(!name.empty(),
                     "NAME attribute must be non-empty in XML element:\n\t'" +
                         tagcontent.str() + "'");

            // generate a list of species.
            std::vector<std::string> species;
            bool valid = ParseUtils::GenerateVector(species, varStrings);

            ASSERTL0(valid, "Unable to process list of variable in XML "
                            "element \n\t'" +
                                tagcontent.str() + "'");

            if (varStrings.size())
            {
                TiXmlElement *info = specie->FirstChildElement("P");

                while (info)
                {
                    tagcontent.clear();
                    tagcontent << *info;
                    // read the property name
                    ASSERTL0(info->Attribute("PROPERTY"),
                             "Missing PROPERTY attribute in "
                             "Species  '" +
                                 name + "' in XML element: \n\t'" +
                                 tagcontent.str() + "'");
                    std::string property = info->Attribute("PROPERTY");
                    ASSERTL0(!property.empty(),
                             "Species properties must have a "
                             "non-empty name for Species '" +
                                 name + "' in XML element: \n\t'" +
                                 tagcontent.str() + "'");

                    // make sure that solver property is capitalised
                    std::string propertyUpper = boost::to_upper_copy(property);

                    // read the value
                    ASSERTL0(info->Attribute("VALUE"),
                             "Missing VALUE attribute in Species '" + name +
                                 "' in XML element: \n\t" + tagcontent.str() +
                                 "'");
                    std::string value = info->Attribute("VALUE");
                    ASSERTL0(!value.empty(), "Species properties must have a "
                                             "non-empty value for Species '" +
                                                 name +
                                                 "' in XML element: \n\t'" +
                                                 tagcontent.str() + "'");

                    // Store values under variable map.
                    for (int i = 0; i < varStrings.size(); ++i)
                    {
                        auto x = GetGloSysSolnList().find(varStrings[i]);
                        if (x == GetGloSysSolnList().end())
                        {
                            (GetGloSysSolnList()[varStrings[i]])
                                [propertyUpper] = value;
                        }
                        else
                        {
                            x->second[propertyUpper] = value;
                        }
                    }
                    info = info->NextSiblingElement("P");
                }
                specie = specie->NextSiblingElement("S");
            }
        }
    }
};
} // namespace NESO::Solvers::tokamak
#endif
