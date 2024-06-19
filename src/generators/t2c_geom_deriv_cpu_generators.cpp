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

#include "t2c_geom_deriv_cpu_generators.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "file_stream.hpp"

#include "t2c_utils.hpp"

void
T2CGeomDerivCPUGenerator::generate(const int                 max_ang_mom,
                                   const std::array<int, 3>& geom_drvs) const
{
    if (geom_drvs[2] == 0)
    {
        for (int i = 0; i <= max_ang_mom; i++)
        {
            const auto integral = _get_integral({i, 0}, geom_drvs);
            
            const auto geom_integrals = t2c::get_geom_integrals(integral);
            
            std::cout << " *** REFERENCE: " << integral.prefix_label() << " | " << integral.label() << " : " << geom_integrals.size() << std::endl;
           
            for (const auto& tint : geom_integrals)
            {
                std::cout << " <>" << tint.prefix_label() << " | " << tint.label() << std::endl;
            }
        }
    }
    else
    {
        for (int i = 0; i <= max_ang_mom; i++)
        {
            for (int j = 0; j <= max_ang_mom; j++)
            {
                const auto integral = _get_integral({i, j}, geom_drvs);
                
                const auto geom_integrals = t2c::get_geom_integrals(integral);
                
                std::cout << " *** REFERENCE: " << integral.prefix_label() << " | " << integral.label() << " : " << geom_integrals.size() << std::endl;
               
                for (const auto& tint : geom_integrals)
                {
                    std::cout << " <>" << tint.prefix_label() << " | " << tint.label() << std::endl;
                }
            }
        }
    }
}

I2CIntegral
T2CGeomDerivCPUGenerator::_get_integral(const std::array<int, 2>& ang_moms,
                                        const std::array<int, 3>& geom_drvs) const
{
    // bra and ket sides
    
    const auto bra = I1CPair("GA", ang_moms[0]);
    
    const auto ket = I1CPair("GB", ang_moms[1]);

    VOperators prefixes;

    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[0])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[2])));

    return I2CIntegral(bra, ket, Operator("1"), 0, prefixes);
}
