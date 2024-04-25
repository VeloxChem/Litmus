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

#include <iostream>

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
V4IElectronRepulsionDriver::ket_hrr(const I4CIntegral& integral) const
{
    SI4CIntegrals tints;
    
    if (!is_electron_repulsion(integral)) return tints;
    
    if (const auto tval = integral.shift(-1, 2))
    {
        // first recursion term

        tints.insert(*tval);
        
        // second recursion term
        
        if (const auto r2val = tval->shift(1, 3))
        {
            tints.insert(*r2val);
        }
    }
        
    return tints;
}

SI4CIntegrals
V4IElectronRepulsionDriver::bra_vrr(const I4CIntegral& integral) const
{
    SI4CIntegrals tints;
    
    if (!is_electron_repulsion(integral)) return tints;
    
    if (const auto tval = integral.shift(-1, 1))
    {
        // first recursion term
        
        tints.insert(*tval);
        
        // second recursion term
        
        if (const auto r2val = tval->shift_order(1))
        {
            tints.insert(*r2val);
        }
        
        // third and fourth recursion terms
        
        if (const auto r3val = tval->shift(-1, 1))
        {
            tints.insert(*r3val);
            
            if (const auto r4val = r3val->shift_order(1))
            {
                tints.insert(*r4val);
            }
        }
        
        // fifth recursion term
        
        if (const auto xval = tval->shift(-1, 3))
        {
            if (const auto r5val = xval->shift_order(1))
            {
                tints.insert(*r5val);
            }
        }
    }
        
    return tints;
}

SI4CIntegrals
V4IElectronRepulsionDriver::ket_vrr(const I4CIntegral& integral) const
{
    SI4CIntegrals tints;
    
    if (!is_electron_repulsion(integral)) return tints;
    
    if (const auto tval = integral.shift(-1, 3))
    {
        // first recursion term
        
        tints.insert(*tval);
        
        // second recursion term
        
        if (const auto r2val = tval->shift_order(1))
        {
            tints.insert(*r2val);
        }
        
        // third and fourth recursion terms
        
        if (const auto r3val = tval->shift(-1, 3))
        {
            tints.insert(*r3val);
            
            if (const auto r4val = r3val->shift_order(1))
            {
                tints.insert(*r4val);
            }
        }
    }
        
    return tints;
}

SI4CIntegrals
V4IElectronRepulsionDriver::apply_bra_hrr_recursion(const I4CIntegral& integral) const
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
                       
                       if (ctint[0] != 0)
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

SI4CIntegrals
V4IElectronRepulsionDriver::apply_ket_hrr_recursion(const I4CIntegral& integral) const
{
    SI4CIntegrals tints;
    
    if (integral[2] > 0)
    {
        SI4CIntegrals rtints({integral, });
                
        while (!rtints.empty())
        {
            SI4CIntegrals new_rtints;
                
            for (const auto& rtint : rtints)
            {
                if (rtint[2] != 0)
                {
                   const auto ctints = ket_hrr(rtint);
                    
                   for (const auto& ctint : ctints)
                   {
                       tints.insert(ctint);
                       
                       if (ctint[2] != 0)
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

SI4CIntegrals
V4IElectronRepulsionDriver::apply_bra_vrr_recursion(const I4CIntegral& integral) const
{
    SI4CIntegrals tints;
    
    if (integral[1] > 0)
    {
        SI4CIntegrals rtints({integral, });
                
        while (!rtints.empty())
        {
            SI4CIntegrals new_rtints;
                
            for (const auto& rtint : rtints)
            {
                if (rtint[1] != 0)
                {
                   const auto ctints = bra_vrr(rtint);
                    
                   for (const auto& ctint : ctints)
                   {
                       tints.insert(ctint);
                       
                       if (ctint[1] != 0)
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

SI4CIntegrals
V4IElectronRepulsionDriver::apply_ket_vrr_recursion(const I4CIntegral& integral) const
{
    SI4CIntegrals tints;
    
    if (integral[3] > 0)
    {
        SI4CIntegrals rtints({integral, });
                
        while (!rtints.empty())
        {
            SI4CIntegrals new_rtints;
                
            for (const auto& rtint : rtints)
            {
                if (rtint[3] != 0)
                {
                   const auto ctints = ket_vrr(rtint);
                    
                   for (const auto& ctint : ctints)
                   {
                       tints.insert(ctint);
                       
                       if (ctint[3] != 0)
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


SI4CIntegrals
V4IElectronRepulsionDriver::create_bra_hrr_recursion(const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    for (const auto& integral : integrals)
    {
        tints.insert(integral);
        
        if (is_electron_repulsion(integral))
        {
            const auto ctints = apply_bra_hrr_recursion(integral);
            
            tints.insert(ctints.cbegin(), ctints.cend());
        }
    }
    
    return tints;
}

SI4CIntegrals
V4IElectronRepulsionDriver::create_ket_hrr_recursion(const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    for (const auto& integral : integrals)
    {
        tints.insert(integral);
        
        if (is_electron_repulsion(integral))
        {
            const auto ctints = apply_ket_hrr_recursion(integral);
            
            tints.insert(ctints.cbegin(), ctints.cend());
        }
    }
    
    return tints;
}

SI4CIntegrals
V4IElectronRepulsionDriver::create_vrr_recursion(const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    for (const auto& integral : integrals)
    {
        tints.insert(integral);
        
        if (is_electron_repulsion(integral))
        {
            if (integral[1] > 0)
            {
                for (const auto& bintegral : apply_bra_vrr_recursion(integral))
                {
                    tints.insert(bintegral);
                    
                    if (bintegral[1] == 0)
                    {
                        const auto ctints = apply_ket_vrr_recursion(bintegral);
                        
                        tints.insert(ctints.cbegin(), ctints.cend());
                    }
                }
            }
            else
            {
                const auto ctints = apply_ket_vrr_recursion(integral);
                
                tints.insert(ctints.cbegin(), ctints.cend());
            }
        }
    }
    
    return tints;
}
