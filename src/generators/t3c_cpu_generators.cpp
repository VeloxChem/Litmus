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

#include "t3c_cpu_generators.hpp"

#include "string_formater.hpp"
#include "file_stream.hpp"

void
T3CCPUGenerator::generate(const std::string& label,
                          const int          max_ang_mom,
                          const int          max_aux_ang_mom) const
{
    if (_is_available(label))
    {
        for (int i = 0; i <= max_aux_ang_mom; i++)
        {
            for (int j = 0; j <= max_ang_mom; j++)
            {
                for (int k = j; k <= max_ang_mom; k++)
                {
                    const auto integral = _get_integral(label, {i, j, k});
                    
                    const auto ket_integrals = _generate_ket_hrr_integral_group(integral);
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
T3CCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "electron repulsion") return true;
    
    return false;
}

I3CIntegral
T3CCPUGenerator::_get_integral(const std::string&        label,
                               const std::array<int, 3>& ang_moms) const
{
    // bra and ket sides
    
    const auto bpair = I1CPair("GA", ang_moms[0]);
    
    const auto kpair = I2CPair("GC", ang_moms[1], "GD", ang_moms[2]);
    
    // electron repulsion integrals
    
    if (fstr::lowercase(label) == "electron repulsion")
    {
        return I3CIntegral(bpair, kpair, Operator("1/|r-r'|"));
    }
    
    return I3CIntegral();
}

SI3CIntegrals
T3CCPUGenerator::_generate_ket_hrr_integral_group(const I3CIntegral& integral) const
{
    SI3CIntegrals tints;
    
    // Electron repulsion integrals
    
    if (integral.integrand() == Operator("1/|r-r'|"))
    {
//        V3IElectronRepulsionDriver eri_drv;
//        
//        if (integral.is_simple())
//        {
//            tints = eri_drv.create_ket_hrr_recursion({integral,});
//        }
//        else
//        {
//            /// TODO: ...
//        }
    }
    
    return tints;
}
