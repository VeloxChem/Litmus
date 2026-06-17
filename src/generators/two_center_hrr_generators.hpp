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

#ifndef two_center_hrr_generators_hpp
#define two_center_hrr_generators_hpp

#include "run_configuration.hpp"

/// Generates the os2c::hrr two-center horizontal-recurrence kernels that fold the
/// Cartesian->spherical transform into the bra/ket momentum transfer. The kernels
/// consume the contracted Cartesian base integrals (reduced to one side at s) and
/// the AB distances, and write the spherical target. The recursion_type selects
/// which transfers to emit: hrr_bra (bra side, la <= lb), hrr_ket (ket side,
/// la > lb), or hrr_bra_ket (both).
class TwoCenterHrrGenerator
{
public:
    /// Creates a two-center HRR kernel generator.
    TwoCenterHrrGenerator() = default;

    /// Writes the kernel .hpp/.cpp pair for every (la|lb) target selected by the
    /// configured recursion_type and angular-momentum range into the current
    /// working directory.
    /// @param run_config The new-style run configuration (recursion_type set).
    void generate(const cfg::RunConfiguration& run_config) const;
};

#endif /* two_center_hrr_generators_hpp */
