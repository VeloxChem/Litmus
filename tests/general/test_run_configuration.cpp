// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.
// E-mail: rinkevic@kth.se
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <gtest/gtest.h>

#include <string>
#include <utility>
#include <vector>

#include "config.hpp"
#include "run_configuration.hpp"

using cfg::ConfigError;
using cfg::Hardware;
using cfg::IntegralType;
using cfg::Language;
using cfg::OperatorType;
using cfg::Signature;
using cfg::StorageForm;

TEST(RunConfigurationTest, AppliesDefaults)
{
    // a minimal config needs only the two required keys
    const auto config = cfg::parse_string(R"(
        integral_type = "four_center"
        max_ang_mom   = 2
    )");

    const auto run_config = cfg::make_run_configuration(config);

    EXPECT_EQ(run_config.integral_type, IntegralType::four_center);
    EXPECT_EQ(run_config.operator_type, OperatorType::overlap);
    EXPECT_EQ(run_config.max_ang_mom, 2);
    EXPECT_EQ(run_config.min_ang_mom, 0);
    EXPECT_EQ(run_config.hardware, Hardware::cpu);
    EXPECT_EQ(run_config.language, Language::cpp);
    EXPECT_EQ(run_config.storage_form, StorageForm::veloxchem_sparse);
    EXPECT_EQ(run_config.signature, Signature::veloxchem_screened);
}

TEST(RunConfigurationTest, ReadsAllExplicitFields)
{
    const auto config = cfg::parse_string(R"(
        integral_type = "three_center"
        operator_type = "electron_repulsion"
        hardware      = "cpu"
        language      = "C++"
        min_ang_mom   = 1
        max_ang_mom   = 4
        storage_form  = "VeloxChemSparse"
        signature     = "VeloxChemScreened"
    )");

    const auto run_config = cfg::make_run_configuration(config);

    EXPECT_EQ(run_config.integral_type, IntegralType::three_center);
    EXPECT_EQ(run_config.operator_type, OperatorType::electron_repulsion);
    EXPECT_EQ(run_config.hardware, Hardware::cpu);
    EXPECT_EQ(run_config.language, Language::cpp);
    EXPECT_EQ(run_config.min_ang_mom, 1);
    EXPECT_EQ(run_config.max_ang_mom, 4);
    EXPECT_EQ(run_config.storage_form, StorageForm::veloxchem_sparse);
    EXPECT_EQ(run_config.signature, Signature::veloxchem_screened);
}

TEST(RunConfigurationTest, NormalizesSpellings)
{
    // case, separators, language alias, and short integral-type forms all match
    const auto config = cfg::parse_string(R"(
        integral_type = "TWO-CENTER"
        language      = "cpp"
        max_ang_mom   = 1
        storage_form  = "veloxchem_sparse"
    )");

    const auto run_config = cfg::make_run_configuration(config);

    EXPECT_EQ(run_config.integral_type, IntegralType::two_center);
    EXPECT_EQ(run_config.language, Language::cpp);
    EXPECT_EQ(run_config.storage_form, StorageForm::veloxchem_sparse);
}

TEST(RunConfigurationTest, ShortIntegralTypeAliases)
{
    for (const auto& [text, expected] : std::vector<std::pair<std::string, IntegralType>>{
             {"2c", IntegralType::two_center},
             {"3c", IntegralType::three_center},
             {"4c", IntegralType::four_center}})
    {
        const auto config = cfg::parse_string("integral_type = \"" + text + "\"\nmax_ang_mom = 1");

        EXPECT_EQ(cfg::make_run_configuration(config).integral_type, expected);
    }
}

TEST(RunConfigurationTest, ParsesEveryOperatorType)
{
    // each accepted spelling (including the space/underscore variants the
    // generators use) maps to the right enum value
    const std::vector<std::pair<std::string, OperatorType>> cases{
        {"overlap", OperatorType::overlap},
        {"kinetic energy", OperatorType::kinetic_energy},
        {"nuclear_potential", OperatorType::nuclear_potential},
        {"electron repulsion", OperatorType::electron_repulsion},
        {"dipole_momentum", OperatorType::dipole_momentum},
        {"linear momentum", OperatorType::linear_momentum},
        {"local", OperatorType::local_ecp},
        {"local_ecp", OperatorType::local_ecp},
        {"projected", OperatorType::projected_ecp},
        {"three center overlap", OperatorType::three_center_overlap},
        {"three_center_r2", OperatorType::three_center_r2},
        {"three center r.r2", OperatorType::three_center_r_dot_r2},
        {"three_center_r_dot_r2", OperatorType::three_center_r_dot_r2}};

    for (const auto& [text, expected] : cases)
    {
        const auto config =
            cfg::parse_string("integral_type = \"two_center\"\nmax_ang_mom = 1\n"
                              "operator_type = \"" + text + "\"");

        EXPECT_EQ(cfg::make_run_configuration(config).operator_type, expected) << text;
    }
}

TEST(RunConfigurationTest, OperatorTypeToStringMatchesGeneratorLabels)
{
    EXPECT_EQ(cfg::to_string(OperatorType::overlap), "overlap");
    EXPECT_EQ(cfg::to_string(OperatorType::kinetic_energy), "kinetic energy");
    EXPECT_EQ(cfg::to_string(OperatorType::nuclear_potential), "nuclear potential");
    EXPECT_EQ(cfg::to_string(OperatorType::electron_repulsion), "electron repulsion");
    EXPECT_EQ(cfg::to_string(OperatorType::dipole_momentum), "dipole momentum");
    EXPECT_EQ(cfg::to_string(OperatorType::linear_momentum), "linear momentum");
    EXPECT_EQ(cfg::to_string(OperatorType::local_ecp), "local");
    EXPECT_EQ(cfg::to_string(OperatorType::projected_ecp), "projected");
    EXPECT_EQ(cfg::to_string(OperatorType::three_center_overlap), "three center overlap");
    EXPECT_EQ(cfg::to_string(OperatorType::three_center_r2), "three center r2");
    EXPECT_EQ(cfg::to_string(OperatorType::three_center_r_dot_r2), "three center r.r2");
}

TEST(RunConfigurationTest, MissingRequiredKeysThrow)
{
    EXPECT_THROW(cfg::make_run_configuration(cfg::parse_string("max_ang_mom = 2")), ConfigError);
    EXPECT_THROW(cfg::make_run_configuration(cfg::parse_string("integral_type = \"two_center\"")),
                 ConfigError);
}

TEST(RunConfigurationTest, UnknownEnumValuesThrow)
{
    const auto bad = [](const std::string& text) {
        return cfg::make_run_configuration(cfg::parse_string(text + "\nmax_ang_mom = 1"));
    };

    EXPECT_THROW(bad("integral_type = \"five_center\""), ConfigError);
    EXPECT_THROW(bad("integral_type = \"two_center\"\noperator_type = \"magnetic\""), ConfigError);
    EXPECT_THROW(bad("integral_type = \"two_center\"\nhardware = \"gpu\""), ConfigError);
    EXPECT_THROW(bad("integral_type = \"two_center\"\nlanguage = \"rust\""), ConfigError);
    EXPECT_THROW(bad("integral_type = \"two_center\"\nstorage_form = \"dense\""), ConfigError);
    EXPECT_THROW(bad("integral_type = \"two_center\"\nsignature = \"plain\""), ConfigError);
}

TEST(RunConfigurationTest, InconsistentAngularMomentumThrows)
{
    EXPECT_THROW(cfg::make_run_configuration(cfg::parse_string(R"(
                     integral_type = "two_center"
                     min_ang_mom   = 3
                     max_ang_mom   = 2
                 )")),
                 ConfigError);

    EXPECT_THROW(cfg::make_run_configuration(cfg::parse_string(R"(
                     integral_type = "two_center"
                     min_ang_mom   = -1
                     max_ang_mom   = 2
                 )")),
                 ConfigError);
}

TEST(RunConfigurationTest, ToStringRoundTrips)
{
    EXPECT_EQ(cfg::to_string(Hardware::cpu), "cpu");
    EXPECT_EQ(cfg::to_string(Language::cpp), "C++");
    EXPECT_EQ(cfg::to_string(IntegralType::two_center), "two_center");
    EXPECT_EQ(cfg::to_string(IntegralType::three_center), "three_center");
    EXPECT_EQ(cfg::to_string(IntegralType::four_center), "four_center");
    EXPECT_EQ(cfg::to_string(StorageForm::veloxchem_sparse), "VeloxChemSparse");
    EXPECT_EQ(cfg::to_string(Signature::veloxchem_screened), "VeloxChemScreened");
}
