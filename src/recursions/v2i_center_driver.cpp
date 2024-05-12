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

#include "v2i_center_driver.hpp"

bool
V2ICenterDriver::is_auxiliary(const I2CIntegral& integral, const int index) const
{
    if (integral.prefixes()[index].shape() == Tensor(0))
    {
        return true;
    }
    else
    {
        return false;
    }
}

SI2CIntegrals
V2ICenterDriver::bra_ket_vrr(const I2CIntegral& integral, const int index) const
{
    SI2CIntegrals tints;

    if (is_auxiliary(integral, index)) return tints;

    if (auto tval = integral.shift_prefix(-1, index, false))
    {
        tval->reduce_prefixes(); 
        
        if (const auto r1val = tval->shift(1, index))
        {
            tints.insert(*r1val);
        }

        if (const auto r2val = tval->shift(-1, index))
        {
            tints.insert(*r2val);
        }
    }
    
    return tints;
}

SI2CIntegrals
V2ICenterDriver::apply_bra_ket_vrr(const I2CIntegral& integral,
                                   const int          index) const
{
    SI2CIntegrals tints;
    
    if (integral.is_simple()) return tints;
    
    if (!is_auxiliary(integral, index))
    {
        SI2CIntegrals rtints({integral, });
                
        while (!rtints.empty())
        {
            SI2CIntegrals new_rtints;
                
            for (const auto& rtint : rtints)
            {
                if (!rtint.is_simple())
                {
                    if (!is_auxiliary(rtint, index))
                    {
                        const auto ctints = bra_ket_vrr(rtint, index);
                        
                        for (const auto& ctint : ctints)
                        {
                            tints.insert(ctint);
                            
                            if (!ctint.is_simple())
                            {
                                if (!is_auxiliary(ctint, index))
                                {
                                    new_rtints.insert(ctint);
                                }
                            }
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
   
    tints.insert(integral);
    
    return tints;
}

SI2CIntegrals
V2ICenterDriver::apply_recursion(const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints;


    
    for (const auto& integral : integrals)
    {
        if (!is_auxiliary(integral, 1))
        {
            for (const auto& bintegral : apply_bra_ket_vrr(integral, 1))
            {
                if (bintegral.prefixes().size() == 2)
                {
                    if (!is_auxiliary(bintegral, 0))
                    {
                        const auto ctints = apply_bra_ket_vrr(bintegral, 0);
                        
                        tints.insert(ctints.cbegin(), ctints.cend());
                    }
                    
                    tints.insert(bintegral);
                }
            }
        }
        
        if ((!is_auxiliary(integral, 0)) && (is_auxiliary(integral, 1)))
        {
            for (const auto& bintegral : apply_bra_ket_vrr(integral, 0))
            {
                tints.insert(bintegral);
            }
        }
    }
        
    return tints;
}
