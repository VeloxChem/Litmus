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

#include "v2i_linmom_driver.hpp"
#include <iostream>

bool
V2ILinearMomentumDriver::is_linmom(const I2CIntegral& integral) const
{
    if (!(integral.prefixes()).empty())
    {
        return false;
    }

    if (integral.integrand() != Operator("p", Tensor(1)))
    {
        return false;
    }
    return true;

}

SI2CIntegrals
V2ILinearMomentumDriver::op_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;

    if (!is_linmom(integral)) return tints;

    integral.replace(Operator("1"));

    // This is a template recursion to generate the overall recursion steps (non-axis-specific recursion steps)
    if (const auto tval = integral.shift(1, 1))
    {
        // first recursion term
        auto x1val = *tval;

        // Stores a representation of an unique integral (insert only actually registers a new entry if one like it wast not there from before)
        tints.insert(x1val);

        // second recursion term
        if (const auto r2val = integral.shift(-1, 1))
        {
            tints.insert(*r2val);
        }

    }

    return tints;
}


SI2CIntegrals
V2ILinearMomentumDriver::apply_op_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;

    if (integral[0] > 0)
    {
        SI2CIntegrals rtints({integral, });

        while (!rtints.empty())
        {
            SI2CIntegrals new_rtints;

            for (const auto& rtint : rtints)
            {
                if ((rtint[0] != 0) && is_linmom(rtint))
                {
                   const auto ctints = op_vrr(rtint);

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

    tints.insert(integral);

    return tints;
}


SI2CIntegrals
V2ILinearMomentumDriver::apply_recursion(const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints;

    for (const auto& integral : integrals)
    {
        tints.insert(integral);

        for (const auto& bintegral : apply_op_vrr(integral))
        {
            tints.insert(bintegral);
        }
    }

    return tints;
}

SI2CIntegrals
V2ILinearMomentumDriver::create_recursion(const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints;

    for (const auto& integral : integrals)
    {
        if (is_linmom(integral))
        {
            const auto ctints = apply_recursion({integral, });

            tints.insert(ctints.cbegin(), ctints.cend());
        }
        else
        {
            tints.insert(integral);
        }
    }

    return tints;
}
