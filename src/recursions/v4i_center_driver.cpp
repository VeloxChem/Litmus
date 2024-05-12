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
    
    auto rtints = tints;
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        for (int i = 0; i < static_cast<int>(prefixes.size()); i++)
        {
            for (int j = 0; j < prefixes[i].shape().order(); j++)
            {
                SI4CIntegrals new_rtints;
                
                for (const auto& rtint : rtints)
                {
                    const auto ctints = bra_ket_vrr(rtint, i);
                    
                    new_rtints.insert(ctints.cbegin(), ctints.cend());
                }
                
                rtints = new_rtints;
                
                tints.insert(new_rtints.begin(), new_rtints.end());
            }
        }
    }
    
    return tints;
}
