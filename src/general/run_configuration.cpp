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

#include "run_configuration.hpp"

#include <cctype>

namespace cfg {  // cfg namespace

namespace {  // unnamed namespace for parsing helpers

/// Normalizes an enum spelling for matching: lowercased, with separators ('_',
/// '-', spaces, '.') removed (so "two_center", "two-center", and "TwoCenter"
/// coincide, as do "kinetic energy" and "kinetic_energy").
std::string
normalize(const std::string& text)
{
    std::string out;

    for (const auto ch : text)
    {
        if ((ch == '_') || (ch == '-') || (ch == '.') || std::isspace(static_cast<unsigned char>(ch)))
        {
            continue;
        }

        out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
    }

    return out;
}

Hardware
parse_hardware(const std::string& value)
{
    const auto key = normalize(value);

    if (key == "cpu") return Hardware::cpu;

    throw ConfigError("config: unknown hardware '" + value + "'; valid: cpu");
}

Language
parse_language(const std::string& value)
{
    const auto key = normalize(value);

    if ((key == "c++") || (key == "cpp")) return Language::cpp;

    throw ConfigError("config: unknown language '" + value + "'; valid: C++");
}

IntegralType
parse_integral_type(const std::string& value)
{
    const auto key = normalize(value);

    if ((key == "twocenter") || (key == "2c")) return IntegralType::two_center;

    if ((key == "threecenter") || (key == "3c")) return IntegralType::three_center;

    if ((key == "fourcenter") || (key == "4c")) return IntegralType::four_center;

    throw ConfigError("config: unknown integral_type '" + value +
                      "'; valid: two_center, three_center, four_center");
}

RecursionType
parse_recursion_type(const std::string& value)
{
    const auto key = normalize(value);

    if (key == "hrrbraket") return RecursionType::hrr_bra_ket;

    if (key == "hrrbra") return RecursionType::hrr_bra;

    if (key == "hrrket") return RecursionType::hrr_ket;

    throw ConfigError("config: unknown recursion_type '" + value +
                      "'; valid: hrr_bra_ket, hrr_bra, hrr_ket");
}

OperatorType
parse_operator_type(const std::string& value)
{
    const auto key = normalize(value);

    if (key == "overlap")           return OperatorType::overlap;
    if (key == "kineticenergy")     return OperatorType::kinetic_energy;
    if (key == "nuclearpotential")  return OperatorType::nuclear_potential;
    if (key == "electronrepulsion") return OperatorType::electron_repulsion;
    if (key == "dipolemomentum")    return OperatorType::dipole_momentum;
    if (key == "linearmomentum")    return OperatorType::linear_momentum;
    if ((key == "local") || (key == "localecp"))         return OperatorType::local_ecp;
    if ((key == "projected") || (key == "projectedecp")) return OperatorType::projected_ecp;
    if (key == "threecenteroverlap") return OperatorType::three_center_overlap;
    if (key == "threecenterr2")      return OperatorType::three_center_r2;
    // "three center r.r2" and "three_center_r_dot_r2" both reach here
    if ((key == "threecenterrr2") || (key == "threecenterrdotr2"))
    {
        return OperatorType::three_center_r_dot_r2;
    }

    throw ConfigError("config: unknown operator_type '" + value +
                      "'; valid: overlap, kinetic_energy, nuclear_potential, "
                      "electron_repulsion, dipole_momentum, linear_momentum, local_ecp, "
                      "projected_ecp, three_center_overlap, three_center_r2, three_center_r_dot_r2");
}

StorageForm
parse_storage_form(const std::string& value)
{
    const auto key = normalize(value);

    if (key == "veloxchemsparse") return StorageForm::veloxchem_sparse;

    throw ConfigError("config: unknown storage_form '" + value + "'; valid: VeloxChemSparse");
}

Signature
parse_signature(const std::string& value)
{
    const auto key = normalize(value);

    if (key == "veloxchemscreened") return Signature::veloxchem_screened;

    throw ConfigError("config: unknown signature '" + value + "'; valid: VeloxChemScreened");
}

}  // namespace

RunConfiguration
make_run_configuration(const Config& config)
{
    RunConfiguration run_config;

    // required: exactly one of integral_type / recursion_type, and max_ang_mom

    const auto has_integral  = config.has("integral_type");
    const auto has_recursion = config.has("recursion_type");

    if (has_integral && has_recursion)
    {
        throw ConfigError("config: 'integral_type' and 'recursion_type' are mutually "
                          "exclusive; specify exactly one");
    }

    if (!has_integral && !has_recursion)
    {
        throw ConfigError("config: exactly one of 'integral_type' or 'recursion_type' "
                          "is required");
    }

    if (has_integral)
    {
        run_config.integral_type = parse_integral_type(config.get_string("integral_type"));
    }
    else
    {
        run_config.recursion_type = parse_recursion_type(config.get_string("recursion_type"));
    }

    run_config.max_ang_mom = config.get_int("max_ang_mom");

    // optional with documented defaults

    run_config.min_ang_mom = config.get_int("min_ang_mom", 0);

    if (config.has("operator_type"))
    {
        run_config.operator_type = parse_operator_type(config.get_string("operator_type"));
    }

    if (config.has("hardware"))
    {
        run_config.hardware = parse_hardware(config.get_string("hardware"));
    }

    if (config.has("language"))
    {
        run_config.language = parse_language(config.get_string("language"));
    }

    if (config.has("storage_form"))
    {
        run_config.storage_form = parse_storage_form(config.get_string("storage_form"));
    }

    if (config.has("signature"))
    {
        run_config.signature = parse_signature(config.get_string("signature"));
    }

    // validate the angular momentum range

    if (run_config.min_ang_mom < 0)
    {
        throw ConfigError("config: 'min_ang_mom' must be non-negative, got " +
                          std::to_string(run_config.min_ang_mom));
    }

    if (run_config.min_ang_mom > run_config.max_ang_mom)
    {
        throw ConfigError("config: 'min_ang_mom' (" + std::to_string(run_config.min_ang_mom) +
                          ") exceeds 'max_ang_mom' (" + std::to_string(run_config.max_ang_mom) + ")");
    }

    return run_config;
}

std::string
to_string(Hardware value)
{
    switch (value)
    {
        case Hardware::cpu: return "cpu";
    }

    return "cpu";
}

std::string
to_string(Language value)
{
    switch (value)
    {
        case Language::cpp: return "C++";
    }

    return "C++";
}

std::string
to_string(IntegralType value)
{
    switch (value)
    {
        case IntegralType::two_center:   return "two_center";
        case IntegralType::three_center: return "three_center";
        case IntegralType::four_center:  return "four_center";
    }

    return "two_center";
}

std::string
to_string(RecursionType value)
{
    switch (value)
    {
        case RecursionType::hrr_bra_ket: return "hrr_bra_ket";
        case RecursionType::hrr_bra:     return "hrr_bra";
        case RecursionType::hrr_ket:     return "hrr_ket";
    }

    return "hrr_bra_ket";
}

std::string
to_string(OperatorType value)
{
    switch (value)
    {
        case OperatorType::overlap:               return "overlap";
        case OperatorType::kinetic_energy:        return "kinetic energy";
        case OperatorType::nuclear_potential:     return "nuclear potential";
        case OperatorType::electron_repulsion:    return "electron repulsion";
        case OperatorType::dipole_momentum:       return "dipole momentum";
        case OperatorType::linear_momentum:       return "linear momentum";
        case OperatorType::local_ecp:             return "local";
        case OperatorType::projected_ecp:         return "projected";
        case OperatorType::three_center_overlap:  return "three center overlap";
        case OperatorType::three_center_r2:       return "three center r2";
        case OperatorType::three_center_r_dot_r2: return "three center r.r2";
    }

    return "overlap";
}

std::string
to_string(StorageForm value)
{
    switch (value)
    {
        case StorageForm::veloxchem_sparse: return "VeloxChemSparse";
    }

    return "VeloxChemSparse";
}

std::string
to_string(Signature value)
{
    switch (value)
    {
        case Signature::veloxchem_screened: return "VeloxChemScreened";
    }

    return "VeloxChemScreened";
}

}  // namespace cfg
