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

#include "two_center_generators.hpp"

#include <iostream>

#include <string>

#include "config.hpp"
#include "operator.hpp"
#include "tensor.hpp"

#include "v2i_ovl_driver.hpp"
#include "v2i_kin_driver.hpp"
#include "v2i_eri_driver.hpp"

namespace {  // diagnostic helpers

/// Formats a two-center integral as bra[operator]ket, e.g. "S[T]P", so that
/// integrals sharing angular momenta but differing in their integrand (overlap
/// vs kinetic energy) print as distinct terms.
/// @param integral The two-center integral.
/// @return The formatted label.
std::string
format_integral(const I2CIntegral& integral)
{
    return Tensor(integral[0]).label() + "[" + integral.integrand().name() + "]" +
           Tensor(integral[1]).label();
}

}  // namespace

I2CIntegral
TwoCenterGenerator::_get_integral(const cfg::RunConfiguration& run_config,
                                  const std::array<int, 2>&    ang_moms) const
{
    // bra (A) and ket (B) expansion centers

    const auto bra = I1CPair("GA", ang_moms[0]);

    const auto ket = I1CPair("GB", ang_moms[1]);

    // select the integrand operator from the configured operator type

    switch (run_config.operator_type)
    {
        case cfg::OperatorType::overlap:
            return I2CIntegral(bra, ket, Operator("1"), 0, {});

        case cfg::OperatorType::kinetic_energy:
            return I2CIntegral(bra, ket, Operator("T"), 0, {});

        case cfg::OperatorType::electron_repulsion:
            return I2CIntegral(bra, ket, Operator("1/|r-r'|"), 0, {});

        // remaining operators are not handled by this routine yet
        case cfg::OperatorType::nuclear_potential:
        case cfg::OperatorType::dipole_momentum:
        case cfg::OperatorType::linear_momentum:
        case cfg::OperatorType::three_center_overlap:
        case cfg::OperatorType::three_center_r2:
        case cfg::OperatorType::three_center_r_dot_r2:
        case cfg::OperatorType::local_ecp:
        case cfg::OperatorType::projected_ecp:
            break;
    }

    throw cfg::ConfigError("two-center generator: operator '" +
                           cfg::to_string(run_config.operator_type) +
                           "' is not supported for two-center integrals");
}

SI2CIntegrals
TwoCenterGenerator::_generate_vrr_base_integral_group(const cfg::RunConfiguration& run_config,
                                                      const I2CIntegral&           integral) const
{
    const auto a = integral[0];  // bra (A) angular momentum

    const auto b = integral[1];  // ket (B) angular momentum

    SI2CIntegrals tints;

    if (a <= b)
    {
        // HRR builds the bra side, so the seeds keep the bra at 0 and ladder the
        // ket from b to a + b: (0|o|b), ..., (0|o|a+b).

        for (int k = b; k <= a + b; k++)
        {
            tints.insert(_get_integral(run_config, {0, k}));
        }
    }
    else
    {
        // HRR builds the ket side, so the seeds keep the ket at 0 and ladder the
        // bra from a to a + b: (a|o|0), ..., (a+b|o|0).

        for (int k = a; k <= a + b; k++)
        {
            tints.insert(_get_integral(run_config, {k, 0}));
        }
    }

    return tints;
}

SI2CIntegrals
TwoCenterGenerator::_generate_vrr_integral_group(const cfg::RunConfiguration& run_config,
                                                 const SI2CIntegrals&         base) const
{
    // recurse the VRR base integrals down to their primitives. The driver chain
    // mirrors T2CCPUGenerator::_generate_integral_group: kinetic and the like
    // close on an overlap recursion.

    switch (run_config.operator_type)
    {
        case cfg::OperatorType::overlap:
        {
            return V2IOverlapDriver().create_recursion(base);
        }

        case cfg::OperatorType::kinetic_energy:
        {
            const auto tints = V2IKineticEnergyDriver().create_recursion(base);

            return V2IOverlapDriver().create_recursion(tints);
        }

        case cfg::OperatorType::electron_repulsion:
        {
            return V2IElectronRepulsionDriver().create_recursion(base);
        }

        // remaining operators are not handled by this routine yet
        case cfg::OperatorType::nuclear_potential:
        case cfg::OperatorType::dipole_momentum:
        case cfg::OperatorType::linear_momentum:
        case cfg::OperatorType::three_center_overlap:
        case cfg::OperatorType::three_center_r2:
        case cfg::OperatorType::three_center_r_dot_r2:
        case cfg::OperatorType::local_ecp:
        case cfg::OperatorType::projected_ecp:
            break;
    }

    throw cfg::ConfigError("two-center generator: operator '" +
                           cfg::to_string(run_config.operator_type) +
                           "' is not supported for two-center integrals");
}

SI2CIntegrals
TwoCenterGenerator::_generate_hrr_integral_group(const cfg::RunConfiguration& run_config,
                                                 const I2CIntegral&           integral) const
{
    const auto a = integral[0];  // bra (A) angular momentum

    const auto b = integral[1];  // ket (B) angular momentum

    SI2CIntegrals tints;

    if (a <= b)
    {
        // grow the bra from 1 to a; at bra order i the ket spans b to a + b - i,
        // closing on the target (a|o|b). The i = 0 row is the VRR seed group.

        for (int i = 1; i <= a; i++)
        {
            for (int j = b; j <= a + b - i; j++)
            {
                tints.insert(_get_integral(run_config, {i, j}));
            }
        }
    }
    else
    {
        // grow the ket from 1 to b; at ket order j the bra spans a to a + b - j,
        // closing on the target (a|o|b). The j = 0 row is the VRR seed group.

        for (int j = 1; j <= b; j++)
        {
            for (int i = a; i <= a + b - j; i++)
            {
                tints.insert(_get_integral(run_config, {i, j}));
            }
        }
    }

    return tints;
}

void
TwoCenterGenerator::generate(const cfg::RunConfiguration& run_config) const
{
    // loop over the angular-momentum range on the bra (A) and ket (B) sides; for
    // each target integral split the work into three groups: the HRR transfer
    // integrals, the VRR base integrals consumed by HRR, and the remaining VRR
    // integrals produced to evaluate that base. The C++ emission is wired in later.

    for (int i = run_config.min_ang_mom; i <= run_config.max_ang_mom; i++)
    {
        for (int j = run_config.min_ang_mom; j <= run_config.max_ang_mom; j++)
        {
            const auto integral = _get_integral(run_config, {i, j});

            const auto hrr_ints = _generate_hrr_integral_group(run_config, integral);

            const auto vrr_base_ints = _generate_vrr_base_integral_group(run_config, integral);

            const auto vrr_ints = _generate_vrr_integral_group(run_config, vrr_base_ints);

            // the base group is, by construction, a subset of the full VRR group;
            // the remainder is what VRR generates to evaluate the base.

            SI2CIntegrals vrr_rest_ints;

            for (const auto& tint : vrr_ints)
            {
                if (vrr_base_ints.count(tint) == 0) vrr_rest_ints.insert(tint);
            }

            std::cout << "(" << format_integral(integral) << ") : HRR {";

            for (const auto& tint : hrr_ints) std::cout << " " << format_integral(tint);

            std::cout << " } VRR-by-HRR {";

            for (const auto& tint : vrr_base_ints) std::cout << " " << format_integral(tint);

            std::cout << " } VRR-rest {";

            for (const auto& tint : vrr_rest_ints) std::cout << " " << format_integral(tint);

            std::cout << " }" << std::endl;
        }
    }
}
