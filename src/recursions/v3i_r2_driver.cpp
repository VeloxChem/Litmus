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

#include "v3i_r2_driver.hpp"

bool
V3IR2Driver::is_r2(const I2CIntegral& integral) const
{
    if (!(integral.prefixes()).empty())
    {
        return false;
    }
    
    if (integral.integrand() != Operator("GR2(r)"))
    {
        return false;
    }
    else
    {
        return true;
    }
}

SI2CIntegrals
V3IR2Driver::aux_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;
    
    if (is_r2(integral))
    {
        auto rint = integral.replace(Operator("G(r)"));
        
        tints.insert(rint);
        
        if (const auto tval = rint.shift(-1, 0))
        {
            tints.insert(*tval);
        }
        
        if (const auto tval = rint.shift(-1, 1))
        {
            tints.insert(*tval);
            
            if (const auto bkval = tval->shift(-1, 0))
            {
                tints.insert(*bkval);
            }
        }
        
        if (const auto tval = rint.shift(-2, 0))
        {
            tints.insert(*tval);
        }
        
        if (const auto tval = rint.shift(-2, 1))
        {
            tints.insert(*tval);
        }
    }
    
    return tints;
}
