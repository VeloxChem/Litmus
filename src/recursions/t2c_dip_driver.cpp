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

#include "t2c_dip_driver.hpp"

#include "axes.hpp"
#include "t2c_ovl_driver.hpp"

/// Needs to be in every file like this for any integral
T2CMultipoleDriver::T2CMultipoleDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T2CMultipoleDriver::is_multipole(const R2CTerm& rterm) const
{
    /// Is this not geometrically differentiatied?
    if (!(rterm.prefixes()).empty())
    {
        return false;
    }

    /// Is this not a multipole integral? (of which dipole is one)
    if (const auto integrand = rterm.integrand(); integrand.name() != "r")
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
T2CMultipoleDriver::bra_vrr(const R2CTerm& rterm,
                            const char     axis) const
{
    // Are we requesting a recursion that corresponds to what this routine can do?
    if (!is_multipole(rterm)) return std::nullopt;

    // This is the recursion
    // You need to do this on paper first and then translate it into code as the one seen here

    // Recursion term is some integral K * (a | M(m) | b)
    // This if line makes that into K * (a - 1 [from "-1"]_(axis) | M(m) | b)
    // Here: Step down in ang mom until reaching an invalid integral (i.e. attempting to generate ang mom < 0)
    // tval is the coefficient resulting from the ang mom step
    if (const auto tval = rterm.shift(axis, -1, 0)) // 0: bra; 1: ket (operator treated separately) (for 2-el: 0:a, 1:b, NB: 2:c (no skip for operator), 3:d)
    {
        R2CDist t2crt(rterm);

        // first recursion term

        auto x1val = *tval;

        const auto coord = _rxyz[axes::to_index(axis)];

        // The answer is generated here for the first recursion term
        // Factor arguments: Run-time here label of factor, name as used in autogen code, rank of data (scalar, tensor etc.) (here: rank 1)
        x1val.add(Factor("PA", "pa", coord), Fraction(1));
// Example: (a | O | b) = C1(a -1, O | b ) [first recursion term] + C2 (a | O | b-1) [second recursion term]
        t2crt.add(x1val);

        // second recursion term

        // This means that this specific recursion involves an "a-2" term
        if (const auto r2val = tval->shift(axis, -1, 0))
        {
            auto x2val = *r2val;

            // Here: Value of the angular momentum along a specific axis
            // here: 0 for bra side, 1 for ket side
            const auto na = x1val[0][axis];

            // Right here this is standard Obara--Saika notation but we choose the names
            x2val.add(Factor("1/eta", "fe"), Fraction(na));

            t2crt.add(x2val);
        }

        // third recursion term

        // a-1, b-1
        if (const auto r3val = tval->shift(axis, -1, 1))
        {
            auto x3val = *r3val;

            const auto nb = x1val[1][axis];

            x3val.add(Factor("1/eta", "fe"), Fraction(nb));

            t2crt.add(x3val);
        }

        // fourth recursion term
        // a-1, operator - 1
        // Here: Special for this operator (possibly other operators too): Step downwards in orders of multipole expansion (e.g. quadrupole to dipole, dipole to identity, identity to fail/stop)
        if (const auto r4val = rterm.shift_operator(axis, -1))
        {
            auto x4val = *r4val;

            // if the operator downwards step reduced it to "scalar", then make the operator be the identity (and then next downwards steo will be the end condition)
            if (x4val.integrand().shape() == TensorComponent(0, 0, 0))
            {
                x4val = x4val.replace(OperatorComponent("1"));
            }

            // This returns the order of the operator (in multipole expansion) along a specific axis
            const auto nc = tval->integrand()[axis];

            x4val.add(Factor("1/eta", "fe"), Fraction(nc));

            t2crt.add(x4val);
        }

        // What is returned is (indirectly through class specification) a tuple of entries (integral, coefficient for recursion associated with that term)

        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}

// THis is a limiting case for when you have stepped down to zero ang mom on the bra side
std::optional<R2CDist>
T2CMultipoleDriver::ket_vrr(const R2CTerm& rterm,
                            const char     axis) const
{
    if (!is_multipole(rterm)) return std::nullopt;

    // This does not change rterm but gives me a tval holder under which to accumulate "terms" (products of factors) (or it gives nothing if the step is invalid)
    if (const auto tval = rterm.shift(axis, -1, 1))
    {
        R2CDist t2crt(rterm);

        // first recursion term

        auto x1val = *tval;

        const auto coord = _rxyz[axes::to_index(axis)];

        // Different center for step, therefore different value/name for coefficient than for bra_vrr
        x1val.add(Factor("PB", "pb", coord), Fraction(1));

        t2crt.add(x1val);

        // second recursion term

        if (const auto r2val = tval->shift(axis, -1, 1))
        {
            auto x2val = *r2val;

            const auto nb = x1val[1][axis];

            x2val.add(Factor("1/eta", "fe"), Fraction(nb));

            t2crt.add(x2val);
        }

        // third recursion term

        if (const auto r3val = rterm.shift_operator(axis, -1))
        {
            auto x3val = *r3val;

            if (x3val.integrand().shape() == TensorComponent(0, 0, 0))
            {
                x3val = x3val.replace(OperatorComponent("1"));
            }

            const auto nc = tval->integrand()[axis];

            x3val.add(Factor("1/eta", "fe"), Fraction(nc));

            t2crt.add(x3val);
        }

        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}


R2CDist
T2CMultipoleDriver::apply_bra_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;

    // Gauge this from pen and paper, can be larger but never smaller
    // At least one more than number of integrals to be returned for a single recursion step (a single use of (here) bra_vrr)
    size_t nints = 5;

    // Depending on which config of ang mom on bra, ket you got in, decide which axis of recursion is most economical and choose to use that one
    for (const auto axis : "xyz")
    {
        if (const auto trec = bra_vrr(rterm, axis))
        {
            if (const auto nterms = trec->terms(); nterms < nints)
            {
                // If this is the current leader have it be the current candidate
                t2crt = *trec;

                // this now becomes the score to beat
                nints = nterms;
            }
        }
    }

    return t2crt;
}

// Same functioning as apply_bra_vrr but this time for ket
R2CDist
T2CMultipoleDriver::apply_ket_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;

    size_t nints = 6;

    for (const auto axis : "xyz")
    {
        if (const auto trec = ket_vrr(rterm, axis))
        {
            if (const auto nterms = trec->terms(); nterms < nints)
            {
                t2crt = *trec;

                nints = nterms;
            }
        }
    }

    return t2crt;
}

// Stops here for May 8 new code but this is still in use on main branch of that date
void
T2CMultipoleDriver::apply_recursion(R2CDist& rdist) const
{
    // vertical recursions on bra side

    apply_bra_vrr(rdist);

    // vertical recursions on ket side

    apply_ket_vrr(rdist);
}

void
T2CMultipoleDriver::apply_bra_vrr(R2CDist& rdist) const
{
    if (!rdist.auxilary(0))
    {
        R2CDist new_dist(rdist.root());

        V2CTerms rec_terms;

        // set up initial terms for recursion expansion

        if (const auto nterms = rdist.terms(); nterms > 0)
        {
            for (size_t i = 0; i < nterms; i++)
            {
                // This is pretty much the only place where a routine like this would change over the different integral cases
                if (const auto rterm = rdist[i]; is_multipole(rterm))
                {
                    if (rterm.auxilary(0))
                    {
                        new_dist.add(rterm);
                    }
                    else
                    {
                        rec_terms.push_back(rterm);
                    }
                }
                else
                {
                    new_dist.add(rterm);
                }
            }
        }
        else
        {
            if (const auto rterm = rdist.root(); is_multipole(rterm))
            {
                rec_terms.push_back(rterm);
            }
        }

        // apply recursion until only elementary integrals are left

        while (!rec_terms.empty())
        {
            V2CTerms new_terms;

            for (size_t i = 0; i < rec_terms.size(); i++)
            {
                const auto cdist = apply_bra_vrr(rec_terms[i]);

                if (const auto nterms = cdist.terms(); nterms > 0)
                {
                    for (size_t j = 0; j < nterms; j++)
                    {
                        if (const auto rterm = cdist[j]; (rterm.auxilary(0) || (!is_multipole(rterm))))
                        {
                            new_dist.add(rterm);
                        }
                        else
                        {
                            new_terms.push_back(rterm);
                        }
                    }
                }
            }

            rec_terms = new_terms;
        }

        // update recursion distribution

        rdist = new_dist;
    }
}

void
T2CMultipoleDriver::apply_ket_vrr(R2CDist& rdist) const
{
    if (!rdist.auxilary(1))
    {
        R2CDist new_dist(rdist.root());

        V2CTerms rec_terms;

        // set up initial terms for recursion expansion

        if (const auto nterms = rdist.terms(); nterms > 0)
        {
            for (size_t i = 0; i < nterms; i++)
            {
                if (const auto rterm = rdist[i]; is_multipole(rterm))
                {
                    if (rterm.auxilary(1))
                    {
                        new_dist.add(rterm);
                    }
                    else
                    {
                        rec_terms.push_back(rterm);
                    }
                }
                else
                {
                    new_dist.add(rterm);
                }
            }
        }
        else
        {
            if (const auto rterm = rdist.root(); is_multipole(rterm))
            {
                rec_terms.push_back(rterm);
            }
        }

        // apply recursion until only

        while (!rec_terms.empty())
        {
            V2CTerms new_terms;

            for (size_t i = 0; i < rec_terms.size(); i++)
            {
                const auto cdist = apply_ket_vrr(rec_terms[i]);

                if (const auto nterms = cdist.terms(); nterms > 0)
                {
                    for (size_t j = 0; j < nterms; j++)
                    {
                        if (const auto rterm = cdist[j]; (rterm.auxilary(1) || (!is_multipole(rterm))))
                        {
                            new_dist.add(rterm);
                        }
                        else
                        {
                            new_terms.push_back(rterm);
                        }
                    }
                }
            }

            rec_terms = new_terms;
        }

        // update recursion distribution

        rdist = new_dist;
    }
}

// Same structure over all integral cases
R2Group
T2CMultipoleDriver::create_recursion(const VT2CIntegrals& vints) const
{
    // create nuclear potential integrals driver

    T2COverlapDriver ovl_drv;

    R2Group r2group;

    for (const auto& tcomp : vints)
    {
        auto rdist = R2CDist(R2CTerm(tcomp));

        apply_recursion(rdist);

        // apply overlap recursion

        ovl_drv.apply_recursion(rdist);

        r2group.add(rdist);
    }

    r2group.simplify();

    return r2group;
}

