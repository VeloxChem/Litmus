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

    /// Builds the HRR integrals for a target two-center integral (a|o|b): the
    /// momentum-transfer triangle that grows the smaller side from the VRR seeds
    /// up to the target (the seed row itself belongs to the VRR base group).
    /// @param run_config The run configuration (selects the integrand operator).
    /// @param integral The target two-center integral (a|o|b).
    /// @return The set of HRR transfer integrals.
    SI2CIntegrals _generate_hrr_integral_group(const cfg::RunConfiguration& run_config,
                                               const I2CIntegral&           integral) const;

    /// Builds the VRR base integrals consumed by HRR for a target (a|o|b): the
    /// seed ladder that keeps the smaller side at zero. (0|o|b)...(0|o|a+b) when
    /// a <= b, and (a|o|0)...(a+b|o|0) otherwise. These are the integrals output
    /// from HRR and fed into the VRR recursion.
    /// @param run_config The run configuration (selects the integrand operator).
    /// @param integral The target two-center integral (a|o|b).
    /// @return The set of VRR base integrals needed by HRR.
    SI2CIntegrals _generate_vrr_base_integral_group(const cfg::RunConfiguration& run_config,
                                                    const I2CIntegral&           integral) const;

    /// Runs the vertical recursion down from the given VRR base integrals,
    /// enumerating every VRR integral required to evaluate them. The returned set
    /// is the union of the base group and the remaining VRR integrals generated to
    /// produce it.
    /// @param run_config The run configuration (selects the recursion drivers).
    /// @param base The VRR base integrals (output from HRR) to recurse down from.
    /// @return The complete set of VRR integrals.
    SI2CIntegrals _generate_vrr_integral_group(const cfg::RunConfiguration& run_config,
                                               const SI2CIntegrals&         base) const;

public:
    /// Creates a two-center integrals code generator.
    TwoCenterGenerator() = default;

    /// Generates the selected two-center integrals over the configured
    /// angular-momentum range [min_ang_mom, max_ang_mom] on the A and B centers.
    /// @param run_config The validated run configuration.
    void generate(const cfg::RunConfiguration& run_config) const;
};

#endif /* two_center_generators_hpp */
