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

#ifndef two_center_generators_hpp
#define two_center_generators_hpp

#include <array>

#include "run_configuration.hpp"
#include "t2c_defs.hpp"

/// Two-center integrals code generator driven by a typed run configuration.
/// This is the new-style counterpart to T2CCPUGenerator: instead of a loose
/// (label, max_ang_mom, geom, rec_form, use_rs) argument list it consumes a
/// validated cfg::RunConfiguration carrying the integral type, integrand
/// operator, angular-momentum range, hardware/language target, and the storage
/// form and signature of the emitted code.
class TwoCenterGenerator
{
    /// Builds the base two-center integral for the configured operator and the
    /// given bra/ket angular momenta.
    /// @param run_config The run configuration (selects the integrand operator).
    /// @param ang_moms The angular momenta of the bra (A) and ket (B) centers.
    /// @return The two-center integral.
    I2CIntegral _get_integral(const cfg::RunConfiguration& run_config,
                              const std::array<int, 2>&    ang_moms) const;

public:
    /// Creates a two-center integrals code generator.
    TwoCenterGenerator() = default;

    /// Generates the selected two-center integrals over the configured
    /// angular-momentum range [min_ang_mom, max_ang_mom] on the A and B centers.
    /// @param run_config The validated run configuration.
    void generate(const cfg::RunConfiguration& run_config) const;
};

#endif /* two_center_generators_hpp */
