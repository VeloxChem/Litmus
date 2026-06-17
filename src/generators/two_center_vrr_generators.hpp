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

#ifndef two_center_vrr_generators_hpp
#define two_center_vrr_generators_hpp

#include "run_configuration.hpp"

/// Generates the two-center overlap vertical-recurrence kernels that build the
/// (s|lb) overlap from the (s|s) seed. The recursion_type selects the flavor:
/// vrr_cartesian emits the single-step Cartesian building blocks (os2c::vrr::ovl,
/// fe retained, per-primitive), vrr_spherical emits the full Cartesian->spherical
/// kernels (os2c::ovl, fe cancelled, contracted).
class TwoCenterVrrGenerator
{
public:
    /// Creates a two-center VRR kernel generator.
    TwoCenterVrrGenerator() = default;

    /// Writes the kernel .hpp/.cpp pair for every ket angular momentum in the
    /// configured range (lb >= 1) into the current working directory.
    /// @param run_config The new-style run configuration (recursion_type set to a
    /// vrr_* value).
    void generate(const cfg::RunConfiguration& run_config) const;
};

#endif /* two_center_vrr_generators_hpp */
