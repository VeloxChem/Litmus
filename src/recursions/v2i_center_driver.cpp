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
#include <iostream>

bool
V2ICenterDriver::is_auxiliary(const I2CIntegral& integral, const int index) const
{
    if (const auto nprefixes = integral.prefixes().size(); index >= nprefixes)
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

    if (!is_auxiliary(integral, index)) return tints;

    if (const auto tval = integral.shift_prefix(-1, index, true))
    {

        if (const auto r1val = tval->shift(1, index))
        {
            auto x1val = *r1val;

            tints.insert(x1val);
        }

        if (const auto r2val = tval->shift(axis, -1, index))
        {
            auto x2val = *r2val;

            tints.insert(x2val);
        }

    }
           return tints;
}

