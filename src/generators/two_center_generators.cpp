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

#include "config.hpp"
#include "operator.hpp"

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

void
TwoCenterGenerator::generate(const cfg::RunConfiguration& run_config) const
{
    // loop over the angular-momentum range on the bra (A) and ket (B) sides and
    // build the base integral for each combination. The recursion generation and
    // C++ emission are wired in later steps.

    for (int i = run_config.min_ang_mom; i <= run_config.max_ang_mom; i++)
    {
        for (int j = run_config.min_ang_mom; j <= run_config.max_ang_mom; j++)
        {
            const auto integral = _get_integral(run_config, {i, j});

            std::cout << "two-center integral: " << integral.label() << std::endl;
        }
    }
}
