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

#ifndef run_configuration_hpp
#define run_configuration_hpp

#include <optional>
#include <string>

#include "config.hpp"

namespace cfg {  // cfg namespace

/// The target hardware a code generator emits for.
enum class Hardware
{
    cpu
};

/// The source language a code generator emits.
enum class Language
{
    cpp
};

/// The number of centers an integral spans.
enum class IntegralType
{
    two_center,
    three_center,
    four_center
};

/// The horizontal-recurrence transfer a recursion-type run generates: the
/// momentum transfer to both centers, to the bra only, or to the ket only.
enum class RecursionType
{
    hrr_bra_ket,
    hrr_bra,
    hrr_ket,
    vrr_cartesian,
    vrr_spherical
};

/// The integrand operator of an integral. The spellings mirror the labels the
/// generators recognize (see to_string), so they round-trip into generated code.
enum class OperatorType
{
    overlap,
    kinetic_energy,
    nuclear_potential,
    electron_repulsion,
    dipole_momentum,
    linear_momentum,
    local_ecp,
    projected_ecp,
    three_center_overlap,
    three_center_r2,
    three_center_r_dot_r2
};

/// The storage layout of the generated integral buffers.
enum class StorageForm
{
    veloxchem_sparse
};

/// The call signature/screening convention of the generated kernels.
enum class Signature
{
    veloxchem_screened
};

/// A validated code-generation run configuration.
///
/// Built from a parsed Config by make_run_configuration(), which applies the
/// documented defaults and rejects unknown values or an inconsistent angular
/// momentum range. Once constructed every field is trusted.
struct RunConfiguration
{
    /// The number of integral centers. Exactly one of integral_type and
    /// recursion_type is set (they are mutually exclusive alternatives).
    std::optional<IntegralType> integral_type;

    /// The horizontal-recurrence transfer to generate. Exactly one of
    /// integral_type and recursion_type is set.
    std::optional<RecursionType> recursion_type;

    /// The integrand operator (default: overlap).
    OperatorType operator_type = OperatorType::overlap;

    /// The target hardware (default: cpu).
    Hardware hardware = Hardware::cpu;

    /// The emitted source language (default: C++).
    Language language = Language::cpp;

    /// The minimum angular momentum (default: 0).
    int min_ang_mom = 0;

    /// The maximum angular momentum (required; no default).
    int max_ang_mom = 0;

    /// The generated buffer storage layout (default: VeloxChemSparse).
    StorageForm storage_form = StorageForm::veloxchem_sparse;

    /// The generated kernel signature convention (default: VeloxChemScreened).
    Signature signature = Signature::veloxchem_screened;
};

/// Builds a validated run configuration from a parsed config.
/// @param config The parsed key/value configuration.
/// @return The validated run configuration (throws ConfigError on a missing
///         required key, an unknown enumerated value, or min > max angular
///         momentum).
RunConfiguration make_run_configuration(const Config& config);

/// @param value The hardware value.
/// @return The canonical string spelling of a hardware value.
std::string to_string(Hardware value);

/// @param value The language value.
/// @return The canonical string spelling of a language value.
std::string to_string(Language value);

/// @param value The integral-type value.
/// @return The canonical string spelling of an integral-type value.
std::string to_string(IntegralType value);

/// @param value The recursion-type value.
/// @return The canonical string spelling of a recursion-type value.
std::string to_string(RecursionType value);

/// @param value The operator-type value.
/// @return The canonical string spelling of an operator-type value (the label
///         the generators recognize, e.g. "electron repulsion").
std::string to_string(OperatorType value);

/// @param value The storage-form value.
/// @return The canonical string spelling of a storage-form value.
std::string to_string(StorageForm value);

/// @param value The signature value.
/// @return The canonical string spelling of a signature value.
std::string to_string(Signature value);

}  // namespace cfg

#endif /* run_configuration_hpp */
