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
    static ParticleSystemSharedPtr create(
        const LU::SessionReaderSharedPtr session,
        const SD::MeshGraphSharedPtr graph)
    {
        return MemoryManager<ReactionSystem>(session, graph);
    }

    static std::string className;

    ReactionSystem(LU::SessionReaderSharedPtr session,
                   SD::MeshGraphSharedPtr graph);

    ~ReactionSystem() override;

    inline void setup_reaction()
    {
        for (auto r : reactions)
        {
            auto reaction =
                LinearReactionBase<num_products_per_parent,
                                   decltype(r.data),
                                   decltype(r.reaction_kernel),
                                   decltype(r.data_calc_obj)>(
                    particle_group->sycl_target, Sym<REAL>("TOT_REACTION_RATE"),
                    Sym<REAL>("WEIGHT"), species_null.get_id(),
                    r.out_states, r.data, r.reaction_kernel,
                    particle_spec, r.data_calc_obj);
            reaction_controller.add_reaction(std::make_shared<decltype(reaction)>(reaction));
        }
    }

    inline void integrate_inner(const double dt_inner) override
    {
        ParticleSystem::integrate_inner(dt_inner);
        reaction_controller->apply_reactions(this->particle_group, dt_inner);
    }

    void ReadParticles() override
    {
        ParticleSystem::ReadParticles();
        ReadReactions();
    }

protected:
    class ProjectTransformation : TransformationStrategy
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
            field_project->project(sub_group, syms, components);
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

    void ReadReactions()
    {
        // Check we actually have a document loaded.
        ASSERTL0(&session->GetDocument(), "No XML document loaded.");

        TiXmlHandle docHandle(&session->GetDocument());
        TiXmlElement *particles;

        // Look for all data in PARTICLES block.
        particles = docHandle.FirstChildElement("NEKTAR")
                        .FirstChildElement("PARTICLES")
                        .Element();

        TiXmlElement *reactions = particles->FirstChildElement("REACTION");
        if (reactions)
        {
            TiXmlElement *reaction = reactions->FirstChildElement("R");

            while (reaction)
            {
                ReadReaction(reaction);
                reaction = reaction->NextSiblingElement("R");
            }
        }
    }

    void ReadReaction(TiXmlElement *reaction)
    {
        nReactions++;
        std::stringstream tagcontent;
        tagcontent << *reaction;

        ASSERTL0(reaction->Attribute("ID"),
                 "Missing ID attribute in Reaction XML "
                 "element: \n\t'" +
                     tagcontent.str() + "'");
        std::string name = reaction->Attribute("NAME");
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
            TiXmlElement *info = reaction->FirstChildElement("P");

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
                ASSERTL0(!property.empty(), "Species properties must have a "
                                            "non-empty name for Species '" +
                                                name +
                                                "' in XML element: \n\t'" +
                                                tagcontent.str() + "'");

                // make sure that solver property is capitalised
                std::string propertyUpper = boost::to_upper_copy(property);

                // read the value
                ASSERTL0(info->Attribute("VALUE"),
                         "Missing VALUE attribute in Species '" + name +
                             "' in XML element: \n\t" + tagcontent.str() + "'");
                std::string value = info->Attribute("VALUE");
                ASSERTL0(!value.empty(), "Species properties must have a "
                                         "non-empty value for Species '" +
                                             name + "' in XML element: \n\t'" +
                                             tagcontent.str() + "'");

                // Store values under variable map.
                for (int i = 0; i < varStrings.size(); ++i)
                {
                    auto x = GetGloSysSolnList().find(varStrings[i]);
                    if (x == GetGloSysSolnList().end())
                    {
                        (GetGloSysSolnList()[varStrings[i]])[propertyUpper] =
                            value;
                    }
                    else
                    {
                        x->second[propertyUpper] = value;
                    }
                }
                info = info->NextSiblingElement("P");
            }
        }
    }
};

} // namespace NESO::Solvers::tokamak
#endif