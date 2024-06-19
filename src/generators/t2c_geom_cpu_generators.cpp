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

#include "t2c_geom_cpu_generators.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "file_stream.hpp"

void
T2CGeomCPUGenerator::generate(const std::string&        label,
                              const int                 max_ang_mom,
                              const std::array<int, 3>& geom_drvs) const
{
    if (_is_available(label))
    {
        for (int i = 0; i <= max_ang_mom; i++)
        {
            for (int j = i; j <= max_ang_mom; j++)
            {
                for (int k = 0; k <= max_ang_mom; k++)
                {
                    for (int l = k; l <= max_ang_mom; l++)
                    {
//                        const auto integral = _get_integral(label, {i, j, k, l}, geom_drvs);
//                        
//                        //const auto geom_integrals = _generate_geom_integral_group(integral);
//                        
//                        const auto geom_integrals = t4c::get_geom_integrals(integral);
//                        
//                        const auto vrr_integrals = _generate_vrr_integral_group(integral, geom_integrals);
//                        
//                        _write_cpp_header(geom_integrals, vrr_integrals, integral);
//                        
//                        std::cout << " *** REFERENCE: " << integral.prefix_label() << " | " << integral.label() << std::endl;
//                        
//                        for (const auto& tint : geom_integrals)
//                        {
//                            std::cout << " <>" << tint.prefix_label() << " | " << tint.label() << std::endl;
//                        }
//                        
//                        std::cout << " --- VRR --- " << std::endl;
//                        
//                        for (const auto& tint : vrr_integrals)
//                        {
//                            std::cout << " <>" << tint.prefix_label() << " | " << tint.label() << "_"  << tint.order() << std::endl;
//                        }
                    }
                }
            }
        }
    }
    else
    {
        std::cerr << "*** ERROR *** Unsupported type of four-center integral: ";
        
        std::cerr << label << " !!!" << std::endl;
        
        std::exit(EXIT_FAILURE);
    }
}

bool
T2CGeomCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "overlap") return true;
    
    if (fstr::lowercase(label) == "kinetic energy") return true;
    
    return false;
}
