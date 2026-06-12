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

#ifndef two_center_emitters_hpp
#define two_center_emitters_hpp

#include <memory>

#include "run_configuration.hpp"
#include "t2c_defs.hpp"

/// Emits the C++ source (a header/definition pair) that computes one two-center
/// integral. An emitter owns the language-, hardware-, and signature-specific
/// formatting choices; TwoCenterGenerator hands it the recursion groups it has
/// already derived and the emitter turns them into source files.
///
/// The concrete emitter for a run is selected by make_two_center_emitter(),
/// which dispatches on the run configuration's hardware and language so an
/// unsupported target is rejected at the seam rather than silently mis-emitted.
class TwoCenterEmitter
{
public:
    virtual ~TwoCenterEmitter() = default;

    /// Emits the source files that compute the target two-center integral.
    /// @param run_config The validated run configuration (selects the signature
    ///        and storage form the emitted kernel branches on).
    /// @param integral The target two-center integral (a|o|b).
    /// @param hrr_ints The HRR transfer integrals that grow the smaller side up
    ///        to the target.
    /// @param vrr_base_ints The VRR base integrals (seeds) HRR consumes.
    /// @param vrr_rest_ints The remaining VRR integrals generated to evaluate the
    ///        base (the full VRR group minus the base).
    virtual void emit(const cfg::RunConfiguration& run_config,
                      const I2CIntegral&           integral,
                      const SI2CIntegrals&         hrr_ints,
                      const SI2CIntegrals&         vrr_base_ints,
                      const SI2CIntegrals&         vrr_rest_ints) const = 0;
};

/// Selects the two-center emitter for a run configuration.
/// @param run_config The validated run configuration.
/// @return The emitter matching the configured hardware and language (throws
///         cfg::ConfigError for a hardware/language combination with no emitter).
std::unique_ptr<TwoCenterEmitter>
make_two_center_emitter(const cfg::RunConfiguration& run_config);

#endif /* two_center_emitters_hpp */
