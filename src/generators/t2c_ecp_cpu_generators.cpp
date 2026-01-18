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

#include "t2c_ecp_cpu_generators.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "file_stream.hpp"

void
T2CECPCPUGenerator::generate(const std::string& label,
                             const int          max_ang_mom) const
{
    if (_is_available(label))
    {
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                for (int i = 0; i <= max_ang_mom; i++)
                {
                    for (int j = 0; j <= max_ang_mom; j++)
                    {
                        #pragma omp task firstprivate(i,j)
                        {
                            const auto integral = _get_integral(label, {i, j});

                            const auto integrals = _generate_integral_group(integral);
                            
                            const auto hrr_integrals = _filter_hrr_integrals(integrals, integral);
                            
                            const auto vrr_integrals = _filter_vrr_integrals(integrals, integral);
                            
                            std::cout << "HRR Integrals : " << integral.label() << " : " << hrr_integrals.size() << std::endl;
                            
                            std::cout << "VRR Integrals : " << integral.label() << " : " << vrr_integrals.size() << std::endl;
                            
//
//                            _write_cpp_header(integrals, integral, rec_form, use_rs);
//                            
//                            if (((i + j) >= 0) && (!use_rs))
//                            {
//                                _write_prim_cpp_header(integral, rec_form);
//                                    
//                                _write_prim_cpp_file(integral);
//                            }
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
T2CECPCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "local") return true;

    if (fstr::lowercase(label) == "projected") return true;
    
    return false;
}

I2CIntegral
T2CECPCPUGenerator::_get_integral(const std::string&        label,
                                  const std::array<int, 2>& ang_moms) const
{
    // bra and ket sides
    
    const auto bra = I1CPair("GA", ang_moms[0]);
    
    const auto ket = I1CPair("GB", ang_moms[1]);
    
    // local core potential
    
    if (fstr::lowercase(label) == "local")
    {
        return I2CIntegral(bra, ket, Operator("U_L"), 0, {});
    }
    
    // projected core potential
    
    if (fstr::lowercase(label) == "projected")
    {
        return I2CIntegral(bra, ket, Operator("U_l"), 0, {});
    }
    
    return I2CIntegral();
}

SI2CIntegrals
T2CECPCPUGenerator::_generate_integral_group(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;

    // Local core potential integrals
    
    if (integral.integrand() == Operator("U_L"))
    {
        if ((integral[0] > 0) && (integral[1] == 0))
        {
            // generate all VRR terms for bra side
            
            for (int i = 0; i <= integral[0]; i++)
            {
                tints.insert(_get_integral("local", {i, 0}));
            }
        }
        
        if ((integral[0] == 0) && (integral[1] > 0))
        {
            // generate all VRR terms for ket side
            
            for (int i = 0; i <= integral[1]; i++)
            {
                tints.insert(_get_integral("local", {0, i}));
            }
        }
        
        if ((integral[0] > 0) && (integral[1] > 0))
        {
            // generate all VRR terms
            
            if (integral[0] > integral[1])
            {
                for (int i = 0; i <= integral[0] + integral[1]; i++)
                {
                    tints.insert(_get_integral("local", {i, 0}));
                }
            } else {
                for (int i = 0; i <= integral[0] + integral[1]; i++)
                {
                    tints.insert(_get_integral("local", {0, i}));
                }
            }
            
            // generate all HRR terms
            
            if (integral[0] > integral[1])
            {
                for (int i = 0; i <= integral[1]; i++)
                {
                    for (int j = integral[0]; j < integral[0] + integral[1] - i;  j++)
                    {
                        tints.insert(_get_integral("local", {j, i}));
                    }
                }
            } else {
                for (int i = 0; i <= integral[0]; i++)
                {
                    for (int j = integral[1]; j < integral[0] + integral[1] - i;  j++)
                    {
                        tints.insert(_get_integral("local", {i, j}));
                    }
                }
            }
            
        }
    }
    
    // Projected core potential 

    return tints;
}

SI2CIntegrals
T2CECPCPUGenerator::_filter_hrr_integrals(const SI2CIntegrals& integrals,
                                          const I2CIntegral&   integral) const
{
    SI2CIntegrals tints;
    
    if (integral.integrand().name() == "U_L")
    {
        if ((integral[0] > 0) && (integral[1] > 0))
        {
            for (const auto& tint : integrals)
            {
                if ((tint[0] > 0) && (tint[1] > 0))
                {
                    tints.insert(tint);
                }
            }
            
            if (integral[0] > integral[1])
            {
                for (const auto& tint : integrals)
                {
                    if ((tint[0] > 0) && (tint[1] > 0))
                    {
                        tints.insert(tint);
                    }
                }
            }
            else
            {
                
            }
        }
    }
    
    return tints;
}

SI2CIntegrals
T2CECPCPUGenerator::_filter_vrr_integrals(const SI2CIntegrals& integrals,
                                          const I2CIntegral&   integral) const
{
    SI2CIntegrals tints;
    
    if (integral.integrand().name() == "U_L")
    {
        for (const auto& tint : integrals)
        {
            if ((tint[0] == 0) || (tint[1] == 0))
            {
                tints.insert(tint);
            }
        }
    }
    
    return tints;
}
