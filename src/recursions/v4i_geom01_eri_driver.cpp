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

#include "v4i_geom01_eri_driver.hpp"

bool
V4IGeom01ElectronRepulsionDriver::is_electron_repulsion(const I4CIntegral& integral) const
{
    if (integral.prefixes_order() != std::vector<int>({0, 1, 0, 0}))
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

SI4CIntegrals
V4IGeom01ElectronRepulsionDriver::bra_hrr(const I4CIntegral& integral) const
{
    SI4CIntegrals tints;
    
    if (!is_electron_repulsion(integral)) return tints;
    
    if (const auto tval = integral.shift(-1, 0))
    {
        // first recursion term

        tints.insert(*tval);
        
        // second recursion term

        tints.insert(tval->base());
        
        // third recursion term
        
        if (const auto r2val = tval->shift(1, 1))
        {
            tints.insert(*r2val);
        }
    }
        
    return tints;
}

SI4CIntegrals
V4IGeom01ElectronRepulsionDriver::apply_bra_hrr_recursion(const I4CIntegral& integral) const
{
    SI4CIntegrals tints;
    
    if (integral[0] > 0)
    {
        SI4CIntegrals rtints({integral, });
                
        while (!rtints.empty())
        {
            SI4CIntegrals new_rtints;
                
            for (const auto& rtint : rtints)
            {
                if (rtint[0] != 0)
                {
                   const auto ctints = bra_hrr(rtint);
                    
                   for (const auto& ctint : ctints)
                   {
                       tints.insert(ctint);
                       
                       if ((ctint[0] != 0) && (!ctint.prefixes().empty()))
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
