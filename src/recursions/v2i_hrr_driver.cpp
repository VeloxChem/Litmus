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

#include "v2i_hrr_driver.hpp"

SI2CIntegrals
V2IHRRDriver::bra_hrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;
     
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

SI2CIntegrals
V2IHRRDriver::ket_hrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;
     
    if (const auto tval = integral.shift(-1, 1))
    {
        // first recursion term

        tints.insert(*tval);
        
        // second recursion term
        
        if (const auto r2val = tval->shift(1, 0))
        {
            tints.insert(*r2val);
        }
    }
        
    return tints;
}

