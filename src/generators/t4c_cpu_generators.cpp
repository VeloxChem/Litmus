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

#include "t4c_cpu_generators.hpp"

#include <iostream>

#include "string_formater.hpp"

#include "v4i_eri_driver.hpp"

void
T4CCPUGenerator::generate(const std::string&        label,
                          const int                 max_ang_mom,
                          const std::array<int, 5>& geom_drvs) const
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
                        const auto integral = _get_integral(label, {i, j, k, l}, geom_drvs);
                        
                        const auto bra_integrals = _generate_bra_hrr_integral_group(integral);
                        
                        const auto ket_integrals = _generate_ket_hrr_integral_group(integral, bra_integrals);
                        
                        auto hrr_integrals = bra_integrals;
                        
                        hrr_integrals.insert(ket_integrals.begin(), ket_integrals.end());
                        
                        const auto vrr_integrals = _generate_vrr_integral_group(integral, hrr_integrals);
                        
                        std::cout << "***  Integral *** " << integral.label() << " : " << integral.integrand().name() << std::endl;
                        
                        std::cout << "-> bra hrr outcome:" << std::endl;
                        
                        for (const auto& tint : bra_integrals)
                        {
                            std::cout << tint.label() << " : " << tint.integrand().name() << std::endl;
                        }
                        
                        std::cout << "-> ket hrr outcome:" << std::endl;
                        
                        for (const auto& tint : ket_integrals)
                        {
                            std::cout << tint.label() << " : " << tint.integrand().name() << std::endl;
                        }
                        
                        std::cout << "-> vrr outcome:" << std::endl;
                        
                        for (const auto& tint : vrr_integrals)
                        {
                            std::cout << tint.label() << " : "  << tint.order() << " : " << tint.integrand().name() << std::endl;
                        }
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
T4CCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "electron repulsion") return true;
    
    return false;
}

I4CIntegral
T4CCPUGenerator::_get_integral(const std::string&        label,
                               const std::array<int, 4>& ang_moms,
                               const std::array<int, 5>& geom_drvs) const
{
    // TODO: Add prefixes support.
    
    // bra and ket sides
    
    const auto bpair = I2CPair("GA", ang_moms[0], "GB", ang_moms[1]);
    
    const auto kpair = I2CPair("GC", ang_moms[2], "GD", ang_moms[3]);
    
    // electron repulsion integrals
    
    if (fstr::lowercase(label) == "electron repulsion")
    {
        return I4CIntegral(bpair, kpair, Operator("1/|r-r'|"));
    }
    
    return I4CIntegral();
}

SI4CIntegrals
T4CCPUGenerator::_generate_bra_hrr_integral_group(const I4CIntegral& integral) const
{
    SI4CIntegrals tints;
    
    // Electron repulsion integrals
    
    if (integral.integrand() == Operator("1/|r-r'|"))
    {
        V4IElectronRepulsionDriver eri_drv;
        
        if (integral.is_simple())
        {
            tints = eri_drv.create_bra_hrr_recursion({integral,});
        }
        else
        {
            /// TODO: ...
        }
    }
    
    return tints;
}

SI4CIntegrals
T4CCPUGenerator::_generate_ket_hrr_integral_group(const I4CIntegral&   integral,
                                                  const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    // Electron repulsion integrals
    
    if (integral.integrand() == Operator("1/|r-r'|"))
    {
        V4IElectronRepulsionDriver eri_drv;
        
        for (const auto& tint : integrals)
        {
            if ((tint[0] == 0) && (tint[2] > 0))
            {
                const auto ctints = eri_drv.create_ket_hrr_recursion({tint, });
                
                tints.insert(ctints.cbegin(), ctints.cend()); 
            }
        }
    }
    
    return tints;
}

SI4CIntegrals
T4CCPUGenerator::_generate_vrr_integral_group(const I4CIntegral&   integral,
                                              const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    // Electron repulsion integrals
    
    if (integral.integrand() == Operator("1/|r-r'|"))
    {
        V4IElectronRepulsionDriver eri_drv;
        
        for (const auto& tint : integrals)
        {
            if ((tint[0] == 0) && (tint[2] == 0))
            {
                const auto ctints = eri_drv.create_vrr_recursion({tint, });
                
                tints.insert(ctints.cbegin(), ctints.cend());
            }
        }
    }
    
    return tints;
}
