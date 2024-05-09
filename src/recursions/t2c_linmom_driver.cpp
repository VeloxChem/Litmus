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

#include "t2c_linmom_driver.hpp"

#include "axes.hpp"
#include "t2c_ovl_driver.hpp"

bool
T2CLinearMomentumDriver::is_linear_momentum(const R2CTerm& rterm) const
{
    /// Is this not geometrically differentiated?
    if (!(rterm.prefixes()).empty())
    {
        return false;
    }

    /// Is this not a multipole integral? (of which dipole is one)
    if (const auto integrand = rterm.integrand(); integrand.name() != "p")
    {
        return false;
    }
    else
    {
        // Is this not a scalar?
        if (integrand.shape() != TensorComponent(0, 0, 0))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

// This is the general case for non-zero (or arbitrary) ang moms for bra and ket sides
// Set up a recursion (same for dipole, quadrupole, ...)
std::optional<R2CDist>
T2CLinearMomentumDriver::op_vrr(const R2CTerm& rterm) const
{
    // Are we requesting a recursion that corresponds to what this routine can do?
    if (!is_linear_momentum(rterm)) return std::nullopt;

    const auto axis = rterm.integrand().shape().primary();

    // This is the recursion
    // You need to do this on paper first and then translate it into code as the one seen here

    R2CDist t2crt(rterm);

    // Recursion term is some integral K * (a | M(m) | b)
    // This if line makes that into K * (a - 1 [from "-1"]_(axis) | M(m) | b)
    // Here: Step down in ang mom until reaching an invalid integral (i.e. attempting to generate ang mom < 0)
    // tval is the coefficient resulting from the ang mom step
    if (const auto tval = rterm.shift(axis, 1, 1)) // 0: bra; 1: ket (operator treated separately) (for 2-el: 0:a, 1:b, NB: 2:c (no skip for operator), 3:d)
    {

        // first recursion term
        auto x1val = *tval;

        x1val.add(Factor("eta", "fz"), Fraction(2));
        x1val = x1val.replace(OperatorComponent("1"));

        t2crt.add(x1val);

    }

    if (const auto tval = rterm.shift(axis, -1, 1))
    {

        // second recursion term
        auto x2val = *tval;

        const auto nb = x2val[1][axis];

        x2val.add(Factor("1", "1"), Fraction(-nb));
        x2val = x2val.replace(OperatorComponent("1"));

        t2crt.add(x2val);

    }

    return t2crt;

}

R2CDist
T2CLinearMomentumDriver::apply_op_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;

    if (const auto trec = op_vrr(rterm))
    {
        // If this is the current leader have it be the current candidate
        t2crt = *trec;

    }

    return t2crt;
}