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

#include "v2i_el_field_driver.hpp"

#include <iostream>

bool
V2IElectricFieldDriver::is_electric_field(const I2CIntegral& integral) const
{
    if (!(integral.prefixes()).empty())
    {
        return false;
    }

    if (integral.integrand() != Operator("A1"))
    {
        return false;
    }
    else
    {
        return true;
    }
}

SI2CIntegrals
V2IElectricFieldDriver::op_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;

    if (!is_electric_field(integral)) return tints;

    if (const auto tval = integral.shift(-1, 0))
    {
        // first recursion term

        tints.insert(*tval);

        // second recursion term

        if (const auto r2val = tval->shift_order(1))
        {
            tints.insert(*r2val);
        }

        // NB: Seventh recursion term
        // MR: Need something else here? This doesn't work here
        // if (const auto r7val = tval->shift_operator(axis, -1))
        // {
        //    tints.insert(*r7val);
        // }

        // third and fourth recursion terms

        if (const auto r3val = tval->shift(-1, 0))
        {
            tints.insert(*r3val);

            if (const auto r4val = r3val->shift_order(1))
            {
                tints.insert(*r4val);
            }
        }

        // fifth and sixth recursion terms

        if (const auto r5val = tval->shift(-1, 1))
        {
            tints.insert(*r5val);

            if (const auto r6val = r5val->shift_order(1))
            {
                tints.insert(*r6val);
            }
        }
    }

    return tints;
}

// MR: Figure this out
I2CIntegral
V2IElectricFieldDriver::aux_vrr(const I2CIntegral& integral) const
{
    if ((integral[0] + integral[1]) == 0)
    {
        auto xint = integral.replace(Operator("1"));

        xint.set_order(0);

        return xint;
    }
    else
    {
        return integral;
    }
}

SI2CIntegrals
V2IElectricFieldDriver::apply_op_vrr(const I2CIntegral& integral) const
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
                if (rtint[0] != 0)
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
V2IElectricFieldDriver::apply_recursion(const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints;

    for (const auto& integral : integrals)
    {
        tints.insert(integral);

        for (const auto& bintegral : apply_op_vrr({integral, }))
        {
            tints.insert(bintegral);
        }
    }

    return tints;
}

SI2CIntegrals
V2IElectricFieldDriver::create_recursion(const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints;

    for (const auto& integral : integrals)
    {
        if (is_electric_field(integral))
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