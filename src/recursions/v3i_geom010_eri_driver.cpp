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

#include "v3i_geom010_eri_driver.hpp"

bool
V3IGeom010ElectronRepulsionDriver::is_electron_repulsion(const I3CIntegral& integral) const
{
    if (integral.prefixes_order() != std::vector<int>({0, 1, 0}))
    {
        return false;
    }
    
    if (integral.integrand() != Operator("1/|r-r'|"))
    {
        return false;
    }
    else
    {
        return true;
    }
}

SI3CIntegrals
V3IGeom010ElectronRepulsionDriver::ket_hrr(const I3CIntegral& integral) const
{
    SI3CIntegrals tints;
    
    if (!is_electron_repulsion(integral)) return tints;
    
    if (const auto tval = integral.shift(-1, 1))
    {
        // first recursion term

        tints.insert(*tval);
        
        // second recursion term

        if (const auto r1val = tval->shift_prefix(-1, 1, false))
        {
            if (r1val->prefixes_order() == std::vector<int>({0, 0, 0}))
            {
                tints.insert(r1val->base());
            }
            else
            {
                tints.insert(*r1val);
            }
        }
        
        // third recursion term
        
        if (const auto r2val = tval->shift(1, 2))
        {
            tints.insert(*r2val);
        }
    }
    
    return tints;
}

SI3CIntegrals
V3IGeom010ElectronRepulsionDriver::ket_aux_hrr(const I3CIntegral& integral) const
{
    SI3CIntegrals tints;
    
    if (!is_electron_repulsion(integral)) return tints;
    
    if (integral[1] > 0) return tints;
    
    if (integral.prefixes_order() == std::vector<int>({0, 1, 0}))
    {
        // first recursion term

        tints.insert(integral.base());
        
        // second recursion term
        
        tints.insert(integral.shift(1, 2)->base());
    }
    
    return tints;
}

SI3CIntegrals
V3IGeom010ElectronRepulsionDriver::apply_ket_hrr_recursion(const I3CIntegral& integral) const
{
    SI3CIntegrals tints;
    
    if (integral[1] > 0)
    {
        SI3CIntegrals rtints({integral, });
                
        while (!rtints.empty())
        {
            SI3CIntegrals new_rtints;
                
            for (const auto& rtint : rtints)
            {
                if (rtint[1] != 0)
                {
                   const auto ctints = ket_hrr(rtint);
                    
                   for (const auto& ctint : ctints)
                   {
                       tints.insert(ctint);
                       
                       if ((ctint[1] != 0) && (!ctint.prefixes().empty()))
                       {
                           new_rtints.insert(ctint);
                       }
                   }
                }
                else
                {
                    tints.insert(rtint);
                }
            }
            
            rtints = new_rtints;
        }
    }
    
    return tints;
}



