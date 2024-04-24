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

#include "v4i_eri_driver.hpp"

bool
V4IElectronRepulsionDriver::is_electron_repulsion(const I4CIntegral& integral) const
{
    if (!(integral.prefixes()).empty())
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
V4IElectronRepulsionDriver::bra_hrr(const I4CIntegral& integral) const
{
    SI4CIntegrals tints;
    
    if (!is_electron_repulsion(integral)) return tints;
    
    if (const auto tval = integral.shift(-1, 0))
    {
        // first recursion term

        tints.insert(*tval);
        
        // second recursion term
        
        if (const auto r2val = tval->shift(1, 1))
        {
            tints.insert(*r2val);
        }
    }
        
    return tints;
}

SI4CIntegrals
V4IElectronRepulsionDriver::apply_bra_hrr_recursion(const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    for (const auto& integral : integrals)
    {
        tints.insert(integral);
        
        if (integral[0] > 0)
        {
            
        }
        else
        {
            tints.insert(integral);
        }
    }
    
    return tints;
}

SI4CIntegrals
V4IElectronRepulsionDriver::create_bra_hrr_recursion(const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    for (const auto& integral : integrals)
    {
        if (is_electron_repulsion(integral))
        {
            const auto ctints = apply_bra_hrr_recursion({integral, });
            
            tints.insert(ctints.cbegin(), ctints.cend());
        }
        else
        {
            tints.insert(integral);
        }
    }
    
    return tints;
}
