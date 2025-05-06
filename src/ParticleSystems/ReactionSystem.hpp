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

    inline void init_spec() override
    {
        ParticleSystem::init_spec();
        this->particle_spec.push(
            ParticleProp(Sym<REAL>("TOT_REACTION_RATE"), 1));
        this->particle_spec.push(ParticleProp(Sym<REAL>("WEIGHT"), 1));
        this->particle_spec.push(
            ParticleProp(Sym<INT>("REACTIONS_PANIC_FLAG"), 1));
        this->particle_spec.push(
            ParticleProp(Sym<INT>("PARTICLE_REACTED_FLAG"), 1));
    }
    template <int n1, int n2>
    inline auto fetch_amjuel_coeffs(const std::string &filename)
    {
        std::ifstream input{filename};

        std::vector<std::vector<float>> csv_rows;
        for (std::string line; std::getline(input, line);)
        {
            std::istringstream ss(std::move(line));
            std::vector<float> row;
            for (std::string value; std::getline(ss, value, ',');)
            {
                row.push_back(std::move(strtod(value.c_str(), 0)));
            }
            csv_rows.push_back(std::move(row));
        }

        std::array<std::array<REAL, n2>, n1> arr;
        for (int i = 0; i < n2; ++i)
        {
            for (int j = 0; j < n1; ++j)
            {
                arr[i][j] = csv_rows[i][j];
            }
        }
        return arr;
    }

    inline void get_normalisation_constants(const std::string &filename)
    {
        std::ifstream input{filename};

        for (std::string line; std::getline(input, line);)
        {
            std::istringstream ss(std::move(line));
            std::vector<std::string> key_value_pair;
            for (std::string value; std::getline(ss, value, ',');)
            {
                key_value_pair.push_back(std::move(value));
            }
            normalisations[key_value_pair[0]] =
                strtod(key_value_pair[1].c_str(), 0);
        }

        bool time_norm_default   = (fabs(normalisations["Time norm"] - -1.0) <
                                  std::numeric_limits<double>::epsilon());
        bool length_norm_default = (fabs(normalisations["Length norm"] - -1.0) <
                                    std::numeric_limits<double>::epsilon());

        REAL vel_norm_default =
            std::sqrt((2 * 1.60217663e-19) / normalisations["Mass norm"]);
        normalisations["Vel norm"] = vel_norm_default;

        if ((time_norm_default) && !(length_norm_default))
        {
            normalisations["Time norm"] =
                normalisations["Length norm"] / vel_norm_default;
        }
        else if (!(time_norm_default) && (length_norm_default))
        {
            normalisations["Length norm"] =
                normalisations["Time norm"] * vel_norm_default;
        }
        else if ((time_norm_default) && (length_norm_default))
        {
            normalisations["Time norm"]   = 1.0;
            normalisations["Length norm"] = 1.0;
        }
        else if (!(time_norm_default) && !(length_norm_default))
        {
            normalisations["Vel norm"] =
                normalisations["Length norm"] / normalisations["Time norm"];
        }
    }

    inline void set_up_reactions()
    {
        auto prop_map = default_map;
        get_normalisation_constants("data/Normalisation_constants.csv");
        REAL mass_amu       = 1.0;
        REAL time_norm      = this->normalisations["Time norm"];
        REAL vel_norm       = this->normalisations["Vel norm"];
        REAL en_norm        = mass_amu * vel_norm * vel_norm;
        REAL length_norm    = this->normalisations["Length norm"];
        REAL temp_norm      = 1.0;
        REAL temp_norm_SI   = this->normalisations["Temperature norm"];
        REAL dens_norm      = this->normalisations["Density norm"];
        REAL amjuel_cs_norm = this->normalisations["AMJUEL Cross-section norm"];
        REAL mass_norm      = this->normalisations["Mass norm"];
        REAL kB             = 1.380649e-23;
        static const int num_coeffs_E  = 9;
        static const int num_coeffs_n  = 9;
        static const int num_coeffs_np = 9;
        static const int num_coeffs_T  = 9;

        const int rank   = particle_group->sycl_target->comm_pair.rank_parent;
        std::mt19937 rng = std::mt19937(52234126 + rank);
        // std::mt19937 rng = std::mt19937(std::random_device{}());
        std::uniform_real_distribution<REAL> uniform_dist(0.0, 1.0);
        auto rng_lambda = [&]() -> REAL
        {
            REAL rng_sample;
            do
            {
                rng_sample = uniform_dist(rng);
            } while (rng_sample == 0.0);
            return rng_sample;
        };
        auto rng_kernel =
            host_atomic_block_kernel_rng<REAL>(rng_lambda, 4 * 10);

        REAL mass_amu_SI = mass_amu * mass_norm;

        for (const auto &[k, v] : this->config->get_reactions())
        {
            if (std::get<0>(v) == "Ionisation")
            {
                auto electron_species = Species("ELECTRON", 5.5e-4, -1.0);
                auto target_species =
                    Species(this->species_map[std::get<1>(v)[0]].name,
                            this->species_map[std::get<1>(v)[0]].mass,
                            this->species_map[std::get<1>(v)[0]].charge,
                            std::get<1>(v)[0]);

                auto test_data = FixedRateData(1.0);

                auto h4_2_1_5_coeffs =
                    fetch_amjuel_coeffs<num_coeffs_n, num_coeffs_T>(
                        "data/H.4_2.1.5.csv");
                auto ionise_data = AMJUEL2DData<num_coeffs_T, num_coeffs_n>(
                    1.0, dens_norm, temp_norm, time_norm, h4_2_1_5_coeffs);

                auto h10_2_1_5_coeffs =
                    fetch_amjuel_coeffs<num_coeffs_np, num_coeffs_T>(
                        "data/H.10_2.1.8.csv");
                auto ionise_energy_data =
                    AMJUEL2DData<num_coeffs_T, num_coeffs_np>(
                        mass_amu * vel_norm * vel_norm, dens_norm, temp_norm,
                        time_norm, h10_2_1_5_coeffs);
                if (this->ndim == 2)
                {
                    auto reaction = ElectronImpactIonisation<
                        decltype(ionise_data), decltype(ionise_energy_data), 2>(
                        this->particle_group->sycl_target, ionise_data,
                        ionise_energy_data, target_species, electron_species,
                        this->particle_spec);

                    this->reaction_controller->add_reaction(
                        std::make_shared<decltype(reaction)>(reaction));
                }
                else if (this->ndim == 3)
                {
                    auto reaction = ElectronImpactIonisation<
                        decltype(ionise_data), decltype(ionise_energy_data), 3>(
                        this->particle_group->sycl_target, ionise_data,
                        ionise_energy_data, target_species, electron_species,
                        this->particle_spec);
                    this->reaction_controller->add_reaction(
                        std::make_shared<decltype(reaction)>(reaction));
                }
            }
            else if (std::get<0>(v) == "Recombination")
            {
                auto normalised_potential_energy =
                    13.6 * 1.60217663e-19 / (mass_amu_SI * vel_norm * vel_norm);
                auto electron_species = Species("ELECTRON", 5.5e-4, -1.0);
                auto neutral_species =
                    Species(this->species_map[std::get<1>(v)[0]].name,
                            this->species_map[std::get<1>(v)[0]].mass,
                            this->species_map[std::get<1>(v)[0]].charge,
                            std::get<1>(v)[0]);

                auto marker_species =
                    Species(this->species_map[std::get<1>(v)[0]].name,
                            this->species_map[std::get<1>(v)[0]].mass,
                            this->species_map[std::get<1>(v)[0]].charge - 1,
                            -1 - std::get<1>(v)[0]);
                // auto recomb_data = FixedRateData(1.0);

                auto h4_2_1_8_coeffs =
                    fetch_amjuel_coeffs<num_coeffs_T, num_coeffs_n>(
                        "data/H.4_2.1.8.csv");
                auto recomb_data = AMJUEL2DData<num_coeffs_T, num_coeffs_n>(
                    1.0, dens_norm, temp_norm, time_norm, h4_2_1_8_coeffs);

                auto h10_2_1_8_coeffs =
                    fetch_amjuel_coeffs<num_coeffs_T, num_coeffs_np>(
                        "data/H.10_2.1.8.csv");
                auto recomb_energy_data =
                    AMJUEL2DData<num_coeffs_T, num_coeffs_np>(
                        mass_amu * vel_norm * vel_norm, dens_norm, temp_norm,
                        time_norm, h10_2_1_8_coeffs);

                auto constant_rate_cross_section =
                    ConstantRateCrossSection(1.0);
                auto recomb_data_calc_sampler = FilteredMaxwellianSampler<
                    2, decltype(constant_rate_cross_section)>(
                    (temp_norm_SI * kB) / (marker_species.get_mass() *
                                           mass_amu_SI * vel_norm * vel_norm),
                    constant_rate_cross_section,
                    rng_kernel); // norm_ratio = temp_norm /
                // (recomb_species.get_mass()
                // * mass_amu * vel_norm *vel_norm)
                auto data_calculator =
                    DataCalculator<decltype(recomb_energy_data),
                                   decltype(recomb_data_calc_sampler)>(
                        particle_spec, recomb_energy_data,
                        recomb_data_calc_sampler);

                // auto data1 = FixedRateData(1.0);
                // auto data2 = FixedRateData(-1.0);

                // auto data_calculator =
                //     DataCalculator<FixedRateData, FixedRateData,
                //     FixedRateData>(
                //         particle_spec, data1, data1, data2);

                if (this->ndim == 2)
                {
                    auto recomb_reaction_kernel = RecombReactionKernels<2>(
                        marker_species, electron_species,
                        normalised_potential_energy);
                    auto reaction =
                        LinearReactionBase<1, decltype(recomb_data),
                                           decltype(recomb_reaction_kernel),
                                           decltype(data_calculator)>(
                            particle_group->sycl_target,
                            marker_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(neutral_species.get_id())},
                            recomb_data, recomb_reaction_kernel, particle_spec,
                            data_calculator);
                    // auto reaction = Recombination<decltype(recomb_data),
                    //                               decltype(data_calculator)>(
                    //     this->particle_group->sycl_target,
                    //     Sym<REAL>(prop_map[default_properties.tot_reaction_rate]),
                    //     Sym<REAL>(prop_map[default_properties.weight]),
                    //     recomb_data, data_calculator, marker_species,
                    //     electron_species, neutral_species,
                    //     this->particle_spec,
                    //     normalised_potential_energy);
                    this->reaction_controller->add_reaction(
                        std::make_shared<decltype(reaction)>(reaction));
                }
                else if (this->ndim == 3)
                {
                    auto recomb_reaction_kernel = RecombReactionKernels<3>(
                        marker_species, electron_species,
                        normalised_potential_energy);
                    auto reaction =
                        LinearReactionBase<1, decltype(recomb_data),
                                           decltype(recomb_reaction_kernel),
                                           decltype(data_calculator)>(
                            particle_group->sycl_target,
                            marker_species.get_id(),
                            std::array<int, 1>{
                                static_cast<int>(neutral_species.get_id())},
                            recomb_data, recomb_reaction_kernel, particle_spec,
                            data_calculator);
                    this->reaction_controller->add_reaction(
                        std::make_shared<decltype(reaction)>(reaction));
                }
            }

            else if (std::get<0>(v) == "ChargeExchange")
            {
                auto projectile_species =
                    Species(this->species_map[std::get<1>(v)[0]].name,
                            this->species_map[std::get<1>(v)[0]].mass,
                            this->species_map[std::get<1>(v)[0]].charge,
                            std::get<1>(v)[0]);
                auto target_species =
                    Species(this->species_map[std::get<1>(v)[1]].name,
                            this->species_map[std::get<1>(v)[1]].mass,
                            this->species_map[std::get<1>(v)[1]].charge,
                            std::get<1>(v)[1]);

                auto rate_data    = FixedRateData(1.0);
                auto vx_beam_data = FixedRateData(1.0);
                auto vy_beam_data = FixedRateData(-1.0);
                auto data_calculator =
                    DataCalculator<FixedRateData, FixedRateData>(
                        particle_spec, vx_beam_data, vy_beam_data);

                if (this->ndim == 2)
                {
                    auto cx_kernel = CXReactionKernels<2>(
                        target_species, projectile_species, prop_map);
                    auto reaction = LinearReactionBase<
                        1, FixedRateData, decltype(cx_kernel),
                        DataCalculator<FixedRateData, FixedRateData>>(
                        this->particle_group->sycl_target,
                        projectile_species.get_id(),
                        std::array<int, 1>{
                            static_cast<int>(target_species.get_id())},
                        rate_data, cx_kernel, this->particle_spec,
                        data_calculator);
                    this->reaction_controller->add_reaction(
                        std::make_shared<decltype(reaction)>(reaction));
                }
                else if (this->ndim == 3)
                {
                    auto cx_kernel = CXReactionKernels<3>(
                        target_species, projectile_species, prop_map);
                    auto reaction = LinearReactionBase<
                        1, FixedRateData, decltype(cx_kernel),
                        DataCalculator<FixedRateData, FixedRateData>>(
                        this->particle_group->sycl_target,
                        projectile_species.get_id(),
                        std::array<int, 1>{
                            static_cast<int>(target_species.get_id())},
                        rate_data, cx_kernel, this->particle_spec,
                        data_calculator);
                    this->reaction_controller->add_reaction(
                        std::make_shared<decltype(reaction)>(reaction));
                }
            }
        }
    }

    inline void integrate_inner(ParticleSubGroupSharedPtr sg,
                                const double dt_inner) override
    {
        ParticleSystem::integrate_inner(sg, dt_inner);
        reaction_controller->apply_reactions(this->particle_group, dt_inner);
    }

    inline void finish_setup(
        std::vector<std::shared_ptr<DisContField>> &src_fields,
        std::vector<Sym<REAL>> &syms, std::vector<int> &components) override
    {
        this->src_syms       = syms;
        this->src_components = components;

        this->field_project = std::make_shared<FieldProject<DisContField>>(
            src_fields, this->particle_group, this->cell_id_translation);

        auto project_transform = std::make_shared<ProjectTransformation>(
            src_fields, this->src_syms, this->src_components,
            this->particle_group, this->cell_id_translation);
        auto project_transform_wrapper =
            std::make_shared<TransformationWrapper>(
                std::dynamic_pointer_cast<TransformationStrategy>(
                    project_transform));

        auto remove_transform =
            std::make_shared<SimpleRemovalTransformationStrategy>();
        auto remove_transform_wrapper = std::make_shared<TransformationWrapper>(
            std::vector<std::shared_ptr<MarkingStrategy>>{make_marking_strategy<
                ComparisonMarkerSingle<REAL, LessThanComp>>(Sym<REAL>("WEIGHT"),
                                                            1e-10)},
            std::dynamic_pointer_cast<TransformationStrategy>(
                remove_transform));

        std::shared_ptr<TransformationStrategy> merge_transform;

        if (this->ndim == 2)
        {
            merge_transform =
                make_transformation_strategy<MergeTransformationStrategy<2>>();
        }
        else if (this->ndim == 3)
        {
            merge_transform =
                make_transformation_strategy<MergeTransformationStrategy<3>>();
        }

        auto merge_transform_wrapper = std::make_shared<TransformationWrapper>(
            std::vector<std::shared_ptr<MarkingStrategy>>{make_marking_strategy<
                ComparisonMarkerSingle<REAL, LessThanComp>>(Sym<REAL>("WEIGHT"),
                                                            0.01)},
            merge_transform);

        std::vector<std::string> src_names;
        for (auto sym : this->src_syms)
        {
            src_names.push_back(sym.name);
        }
        auto zeroer_transform =
            make_transformation_strategy<ParticleDatZeroer<REAL>>(src_names);

        auto zeroer_transform_wrapper = std::make_shared<TransformationWrapper>(
            std::dynamic_pointer_cast<TransformationStrategy>(
                zeroer_transform));

        this->reaction_controller = std::make_shared<ReactionController>(
            std::vector<std::shared_ptr<TransformationWrapper>>{
                zeroer_transform_wrapper, merge_transform_wrapper,
                remove_transform_wrapper},
            std::vector<std::shared_ptr<TransformationWrapper>>{
                project_transform_wrapper, merge_transform_wrapper,
                remove_transform_wrapper});

        set_up_reactions();

        init_output("particle_trajectory.h5part", Sym<REAL>("POSITION"),
                    Sym<INT>("INTERNAL_STATE"), Sym<INT>("CELL_ID"),
                    Sym<REAL>("VELOCITY"), Sym<REAL>("B0"), Sym<REAL>("B1"),
                    Sym<REAL>("B2"), Sym<REAL>("ELECTRON_DENSITY"),
                    this->src_syms, Sym<REAL>("WEIGHT"), Sym<INT>("ID"));
    }

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

    // normalisations
    std::map<std::string, double> normalisations;
};

} // namespace NESO::Solvers::tokamak
#endif