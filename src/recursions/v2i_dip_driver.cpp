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

#include "v2i_dip_driver.hpp"
#include <iostream>

bool
V2IDipoleDriver::is_dipole(const I2CIntegral& integral) const
{
    if (!(integral.prefixes()).empty())
    {
        return false;
    }

    if (integral.integrand() != Operator("r", Tensor(1)))
    {
        return false;
    }
    return true;

}

SI2CIntegrals
V2IDipoleDriver::bra_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;

    if (!is_dipole(integral)) return tints;

    // This is a template recursions to generate the overall recursion steps (non-axis-specific recursion steps)
    if (const auto tval = integral.shift(-1, 0))
    {
        // first recursion term
        // a-1,
        auto x1val = *tval;

        // Stores a represetation of an unique integral (insert only actually registers a new entry if one like it wast not there from before)
        tints.insert(x1val);

        // second recursion term
// a-2,
        if (const auto r2val = tval->shift(-1, 0))
        {
            tints.insert(*r2val);
        }

        // third recursion term
// a-1, b -1
        if (const auto r3val = tval->shift(-1, 1))
        {
            tints.insert(*r3val);
        }
        
        // fourth recursion term
        
        tints.insert(integral.replace(Operator("1")));
    }

    return tints;
}

// Sim to bra but here special case for all bra ang moms at zero
SI2CIntegrals
V2IDipoleDriver::ket_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;

    if (!is_dipole(integral)) return tints;

    if (const auto tval = integral.shift(-1, 1))
    {
        // first recursion term

        auto x1val = *tval;

        tints.insert(x1val);

        // second recursion term

        if (const auto r2val = tval->shift(-1, 1))
        {
            tints.insert(*r2val);
        }

        // third recursion term

        tints.insert(integral.replace(Operator("1")));

    }

    return tints;
}

SI2CIntegrals
V2IDipoleDriver::apply_bra_vrr(const I2CIntegral& integral) const
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
                if ((rtint[0] != 0) && is_dipole(rtint))
                {
                   const auto ctints = bra_vrr(rtint);

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
V2IDipoleDriver::apply_ket_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;

    if (integral[1] > 0)
    {
        SI2CIntegrals rtints({integral, });

        while (!rtints.empty())
        {
            SI2CIntegrals new_rtints;

            for (const auto& rtint : rtints)
            {
                if ((rtint[1] != 0) && is_dipole(rtint))
                {
                   const auto ctints = ket_vrr(rtint);

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

    tints.insert(integral);

    return tints;
}

SI2CIntegrals
V2IDipoleDriver::apply_recursion(const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints;

    for (const auto& integral : integrals)
    {
        tints.insert(integral);

        for (const auto& bintegral : apply_bra_vrr(integral))
        {
            if (bintegral[0] == 0)
            {
                if (bintegral[1] != 0)
                {
                    const auto ctints = apply_ket_vrr(bintegral);

                    tints.insert(ctints.cbegin(), ctints.cend());
                }
                else
                {
                    tints.insert(bintegral);

                    if (is_dipole(bintegral))
                    {
                        tints.insert(bintegral.replace(Operator("1")));
                    }
                }
            }
            else
            {
                tints.insert(bintegral);
            }
        }
    }

    return tints;
}

SI2CIntegrals
V2IDipoleDriver::create_recursion(const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints;

    for (const auto& integral : integrals)
    {
        if (is_dipole(integral))
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
