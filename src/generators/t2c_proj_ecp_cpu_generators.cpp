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

#include "t2c_proj_ecp_cpu_generators.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "file_stream.hpp"
#include "t2c_utils.hpp"
#include "v2i_proj_ecp_driver.hpp"

void
T2CProjECPCPUGenerator::generate(const std::string& label,
                                 const int          max_ang_mom,
                                 const int          proj_ang_mom) const
{
    if (_is_available(label))
    {
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                for (int l = 0; l <= proj_ang_mom; l++)
                {
                    for (int i = 0; i <= max_ang_mom; i++)
                    {
                        for (int j = 0; j <= max_ang_mom; j++)
                        {
                            #pragma omp task firstprivate(i,j)
                            {
                                const auto integral = _get_integral(label, {i, j}, l);
                                
                                const auto integrals = _generate_integral_group(integral);
                                
                                std::cout << " *** " << integral.second.label() << "_" << integral.second.order() << " *** " << std::endl;
                                
                                for (const auto& [order, tint] : integrals)
                                {
                                    std::cout << "> " << tint.label() << "_" << tint.order() << " : " ;
                                    
                                    std::cout << "(" << order[0] << ",";
                                    
                                    std::cout << order[1] << ","  << order[2] << ")" << std::endl;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    else
    {
        std::cerr << "*** ERROR *** Unsupported type of two-center ECP integral: ";
        
        std::cerr << label << " !!!" << std::endl;
        
        std::exit(EXIT_FAILURE);
    }
}

bool
T2CProjECPCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "projected") return true;
    
    return false;
}

M2Integral
T2CProjECPCPUGenerator::_get_integral(const std::string&        label,
                                      const std::array<int, 2>& ang_moms,
                                      const int                 proj_ang_mom) const
{
    // bra and ket sides
    
    const auto bra = I1CPair("GA", ang_moms[0]);
    
    const auto ket = I1CPair("GB", ang_moms[1]);
    
    // projected core potential
    
    if (fstr::lowercase(label) == "projected")
    {
        M2Integral tint = {{0, 0, 0}, I2CIntegral(bra, ket, Operator("U_l"), proj_ang_mom, {})};
        
        return tint;
    }
    
    M2Integral tint = {{0, 0, 0}, I2CIntegral()};
    
    return tint;
}

SM2Integrals
T2CProjECPCPUGenerator::_generate_integral_group(const M2Integral& integral) const
{
    SM2Integrals tints;

    // Projected potential integrals
    
    if (integral.second.integrand() == Operator("U_l"))
    {
        V2IProjectedECPDriver ecp_drv;
        
        if (integral.second.is_simple())
        {
            tints = ecp_drv.create_recursion({integral,});
        }
        else
        {
            tints = ecp_drv.create_recursion(tints);
        }
    }
    
    return tints;
}
