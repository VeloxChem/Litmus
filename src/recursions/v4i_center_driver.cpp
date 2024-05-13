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

#include "v4i_center_driver.hpp"

bool
V4ICenterDriver::is_auxiliary(const I4CIntegral& integral, const int index) const
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

SI4CIntegrals
V4ICenterDriver::bra_ket_vrr(const I4CIntegral& integral, const int index) const
{
    SI4CIntegrals tints;

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

/// Applies vertical recursion to bra side of overlap integral.
/// @param integral The  overlap integral.
/// @return The recursion expansion of integral.
SI4CIntegrals
V4ICenterDriver::apply_bra_ket_vrr(const I4CIntegral& integral) const
{
    SI4CIntegrals tints;
    
    tints.insert(integral);
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        // apply recursion on A center
        
        auto rtints = tints;
        
        SI4CIntegrals a_tints;
        
        for (int i = 0; i < prefixes[0].shape().order(); i++)
        {
            SI4CIntegrals new_rtints;
            
            for (const auto& rtint : rtints)
            {
                const auto ctints = bra_ket_vrr(rtint, 0);
                
                new_rtints.insert(ctints.cbegin(), ctints.cend());
            }
            
            rtints = new_rtints;
            
            a_tints.insert(new_rtints.begin(), new_rtints.end());
        }
        
        tints.insert(a_tints.begin(), a_tints.end());
       
        if (a_tints.empty()) a_tints = tints;
        
        // apply recursion on B center
        
        rtints = a_tints;
        
        SI4CIntegrals b_tints;
        
        for (int i = 0; i < prefixes[1].shape().order(); i++)
        {
            SI4CIntegrals new_rtints;
            
            for (const auto& rtint : rtints)
            {
                const auto ctints = bra_ket_vrr(rtint, 1);
                
                new_rtints.insert(ctints.cbegin(), ctints.cend());
            }
            
            rtints = new_rtints;
            
            b_tints.insert(new_rtints.begin(), new_rtints.end());
        }
        
        tints.insert(b_tints.begin(), b_tints.end());
        
        if (b_tints.empty()) b_tints = a_tints;
        
        // apply recursion on C center
        
        rtints = b_tints;
        
        SI4CIntegrals c_tints;
        
        for (int i = 0; i < prefixes[2].shape().order(); i++)
        {
            SI4CIntegrals new_rtints;
            
            for (const auto& rtint : rtints)
            {
                const auto ctints = bra_ket_vrr(rtint, 2);
                
                new_rtints.insert(ctints.cbegin(), ctints.cend());
            }
            
            rtints = new_rtints;
            
            c_tints.insert(new_rtints.begin(), new_rtints.end());
        }
        
        tints.insert(c_tints.begin(), c_tints.end());
        
        if (c_tints.empty()) c_tints = b_tints;
        
        // apply recursion on D center
        
        rtints = c_tints;
        
        SI4CIntegrals d_tints;
        
        for (int i = 0; i < prefixes[3].shape().order(); i++)
        {
            SI4CIntegrals new_rtints;
            
            for (const auto& rtint : rtints)
            {
                const auto ctints = bra_ket_vrr(rtint, 3);
                
                new_rtints.insert(ctints.cbegin(), ctints.cend());
            }
            
            rtints = new_rtints;
            
            d_tints.insert(new_rtints.begin(), new_rtints.end());
        }
        
        tints.insert(d_tints.begin(), d_tints.end());
    }
    
    return tints;
}
