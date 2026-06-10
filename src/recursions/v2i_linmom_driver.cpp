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

    // The linear-momentum operator p lowers to the overlap operator "1" on the
    // ket-shifted integrals. Integral::replace is const and returns a new value,
    // so the replacement must be applied to each shifted term.
    if (const auto tval = integral.shift(1, 1))
    {
        // first recursion term
        tints.insert(tval->replace(Operator("1")));

        // second recursion term
        if (const auto r2val = integral.shift(-1, 1))
        {
            tints.insert(r2val->replace(Operator("1")));
        }
    }

    return tints;
}


SI2CIntegrals
V2ILinearMomentumDriver::apply_op_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;

    if (is_linmom(integral))
    {
        SI2CIntegrals rtints({integral, });

        while (!rtints.empty())
        {
            SI2CIntegrals new_rtints;

            for (const auto& rtint : rtints)
            {
                if (is_linmom(rtint))
                {
                   const auto ctints = op_vrr(rtint);

                   for (const auto& ctint : ctints)
                   {
                       tints.insert(ctint);

                       if (is_linmom(ctint))
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
