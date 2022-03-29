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

#include "eri_driver.hpp"

#include <array>

#include "axes.hpp"
#include "integral.hpp"
#include "two_center_pair.hpp"

using I2CPair = TwoCenterPair;

using I4CIntegral = Integral<I2CPair, I2CPair>;

EriDriver::EriDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

std::optional<R4CDist>
EriDriver::bra_hrr(const R4CTerm& rterm,
                   const char     axis) const
{
    if (const auto tval = rterm.shift(axis, -1, 0))
    {
        R4CDist t4crt(rterm);
        
        // first recursion term
        
        if (const auto r1val = tval->shift(axis, 1, 1))
        {
            t4crt.add(*r1val);
        }
        
        // second recursion term
        
        auto r2val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
    
        r2val.add(Factor("AB", "rab", coord), Fraction(-1));
        
        t4crt.add(r2val);
        
        // done with recursion terms
        
        return t4crt;
    }
    else
    {
        return std::nullopt;
    }
}

std::optional<R4CDist>
EriDriver::ket_hrr(const R4CTerm& rterm,
                   const char     axis) const
{
    if (const auto tval = rterm.shift(axis, -1, 2))
    {
        R4CDist t4crt(rterm);
        
        // first recursion term
        
        if (const auto r1val = tval->shift(axis, 1, 3))
        {
            t4crt.add(*r1val);
        }
        
        // second recursion term
        
        auto r2val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
    
        r2val.add(Factor("CD", "rcd", coord), Fraction(-1));
        
        t4crt.add(r2val);
        
        // done with recursion terms
        
        return t4crt;
    }
    else
    {
        return std::nullopt;
    }
}

std::optional<R4CDist>
EriDriver::bra_vrr(const R4CTerm& rterm,
                   const char     axis) const
{
    if (const auto tval = rterm.shift(axis, -1, 1))
    {
        R4CDist t4crt(rterm);
        
        // first recursion term
        
        auto r1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        r1val.add(Factor("PB", "rpb", coord), Fraction(1));
        
        t4crt.add(r1val);
        
        // second recursion term
        
        if (const auto r2val = tval->shift_order(1))
        {
            auto x2val = *r2val;
            
            x2val.add(Factor("WP", "rwp", coord), Fraction(1));
            
            t4crt.add(x2val);
        }
        
        // third and fourth recursion terms
        
        if (const auto r3val = tval->shift(axis, -1, 1))
        {
            auto x3val = *r3val;
            
            const auto nb = r1val[1][axis];
            
            x3val.add(Factor("1/zeta", "fz"), Fraction(nb, 2));
            
            t4crt.add(x3val);
            
            if (const auto r4val = r3val->shift_order(1))
            {
                auto x4val = *r4val;
                
                x4val.add(Factor("rho/zeta^2", "frz2"), Fraction(-nb, 2));
                
                t4crt.add(x4val);
            }
        }
        
        // fifth recursion term
        
        if (const auto xval = tval->shift(axis, -1, 3))
        {
            if (const auto r5val = xval->shift_order(1))
            {
                auto x5val = *r5val;
               
                const auto nd = r1val[3][axis];
               
                x5val.add(Factor("1/(zeta+eta)", "fze"), Fraction(nd, 2));
               
                t4crt.add(x5val);
            }
        }
        
        return t4crt;
    }
    else
    {
        return std::nullopt;
    }
}

std::optional<R4CDist>
EriDriver::ket_vrr(const R4CTerm& rterm,
                   const char     axis) const
{
    if (const auto tval = rterm.shift(axis, -1, 3))
    {
        R4CDist t4crt(rterm);
        
        // first recursion term
        
        auto r1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        r1val.add(Factor("QD", "rqd", coord), Fraction(1));
        
        t4crt.add(r1val);
        
        // second recursion term
        
        if (const auto r2val = tval->shift_order(1))
        {
            auto x2val = *r2val;
            
            x2val.add(Factor("WQ", "rwq", coord), Fraction(1));
            
            t4crt.add(x2val);
        }
        
        // third and fourth recursion terms
        
        if (const auto r3val = tval->shift(axis, -1, 3))
        {
            auto x3val = *r3val;
            
            const auto nd = r1val[3][axis];
            
            x3val.add(Factor("1/eta", "fe"), Fraction(nd, 2));
            
            t4crt.add(x3val);
            
            if (const auto r4val = r3val->shift_order(1))
            {
                auto x4val = *r4val;
                
                x4val.add(Factor("rho/eta^2", "fre2"), Fraction(-nd, 2));
                
                t4crt.add(x4val);
            }
        }
        
        return t4crt;
    }
    else
    {
        return std::nullopt;
    }
}

R4CDist
EriDriver::apply_bra_hrr(const R4CTerm&       rterm,
                               ST4CIntegrals& sints) const
{
    R4CDist t4crt;
    
    int nints = 3;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = bra_hrr(rterm, axis))
        {
            const auto nterms = trec->count_new_integrals(sints);
            
            if (nterms < nints)
            {
                t4crt = *trec;
                
                nints = nterms;
            }
        }
    }
    
    const auto vints = t4crt.unique_integrals();
    
    sints.insert(vints.cbegin(), vints.cend());

    return t4crt;
}

R4CDist
EriDriver::apply_ket_hrr(const R4CTerm&       rterm,
                               ST4CIntegrals& sints) const
{
    R4CDist t4crt;
    
    int nints = 3;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = ket_hrr(rterm, axis))
        {
            const auto nterms = trec->count_new_integrals(sints);
            
            if (nterms < nints)
            {
                t4crt = *trec;
                
                nints = nterms;
            }
        }
    }
    
    const auto vints = t4crt.unique_integrals();
    
    sints.insert(vints.cbegin(), vints.cend());

    return t4crt;
}

R4CDist
EriDriver::apply_bra_vrr(const R4CTerm&       rterm,
                                ST4CIntegrals& sints) const
{
    R4CDist t4crt;
    
    int nints = 6;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = bra_vrr(rterm, axis))
        {
            const auto nterms = trec->count_new_integrals(sints);
            
            if (nterms < nints)
            {
                t4crt = *trec;
                
                nints = nterms;
            }
        }
    }
    
    const auto vints = t4crt.unique_integrals();
    
    sints.insert(vints.cbegin(), vints.cend());

    return t4crt;
}

R4CDist
EriDriver::apply_ket_vrr(const R4CTerm&       rterm,
                                ST4CIntegrals& sints) const
{
    R4CDist t4crt;
    
    int nints = 6;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = ket_vrr(rterm, axis))
        {
            const auto nterms = trec->count_new_integrals(sints);
            
            if (nterms < nints)
            {
                t4crt = *trec;
                
                nints = nterms;
            }
        }
    }
    
    const auto vints = t4crt.unique_integrals();
    
    sints.insert(vints.cbegin(), vints.cend());

    return t4crt;
}

R4Group
EriDriver::apply_bra_hrr(const V4CTerms&      rterms,
                               ST4CIntegrals& sints) const
{
    R4Group rgroup;
    
    for (const auto& rterm : rterms)
    {
        if (const auto tval = apply_bra_hrr(rterm, sints); tval.terms() > 0)
        {
            rgroup.add(tval);
        }
    }
    
    return rgroup;
}

R4Group
EriDriver::apply_ket_hrr(const V4CTerms&      rterms,
                               ST4CIntegrals& sints) const
{
    R4Group rgroup;
    
    for (const auto& rterm : rterms)
    {
        if (const auto tval = apply_ket_hrr(rterm, sints); tval.terms() > 0)
        {
            rgroup.add(tval);
        }
    }
    
    return rgroup;
}

R4Group
EriDriver::apply_bra_vrr(const V4CTerms&      rterms,
                               ST4CIntegrals& sints) const
{
    R4Group rgroup;
    
    for (const auto& rterm : rterms)
    {
        if (const auto tval = apply_bra_vrr(rterm, sints); tval.terms() > 0)
        {
            rgroup.add(tval);
        }
    }
    
    return rgroup;
}

R4Group
EriDriver::apply_ket_vrr(const V4CTerms&      rterms,
                               ST4CIntegrals& sints) const
{
    R4Group rgroup;
    
    for (const auto& rterm : rterms)
    {
        if (const auto tval = apply_ket_vrr(rterm, sints); tval.terms() > 0)
        {
            rgroup.add(tval);
        }
    }
    
    return rgroup;
}

void
EriDriver::apply_bra_hrr(R4Graph&       rgraph,
                         ST4CIntegrals& sints) const
{
    // special cases: single vertices without expansion terms
    
    for (const auto i : rgraph.orphans())
    {
        if (rgraph[i].empty() && (!rgraph[i].auxilary(0)))
        {
            const auto rgroup = apply_bra_hrr(rgraph[i].roots(), sints);

            rgraph.replace(rgroup, i);
        }
    }
    
    // loop over orphaned vertices
    
    bool use_hrr = true;
    
    while (use_hrr)
    {
        int cnt = 0;
        
        const auto vidx = rgraph.orphans();
        
        for (const auto i : vidx)
        {
            if (!rgraph[i].auxilary(0))
            {
                for (const auto& vterms : rgraph[i].split_terms<I4CIntegral>())
                {
                    auto rgroup = apply_bra_hrr(vterms, sints);
                    
                    if (rgroup.expansions() == 0)
                    {
                        for (const auto& tval : vterms)
                        {
                            rgroup.add(R4CDist(tval));
                        }
                    }
                    
                    rgraph.add(rgroup, i);
                    
                    cnt++;
                }
            }
        }
        
        if (cnt == 0) use_hrr = false;
    }
    
    rgraph.reduce();
    
    rgraph.sort<I4CIntegral>(true);
}

void
EriDriver::apply_ket_hrr(R4Graph&       rgraph,
                         ST4CIntegrals& sints) const
{
    // special cases: single vertices without expansion terms
    
    for (const auto i : rgraph.orphans())
    {
        if (rgraph[i].empty() && (!rgraph[i].auxilary(2)))
        {
            const auto rgroup = apply_ket_hrr(rgraph[i].roots(), sints);

            rgraph.replace(rgroup, i);
        }
    }
    
    // loop over orphaned vertices
    
    bool use_hrr = true;
    
    while (use_hrr)
    {
        int cnt = 0;
        
        const auto vidx = rgraph.orphans();
        
        for (const auto i : vidx)
        {
            if (!rgraph[i].auxilary(2))
            {
                for (const auto& vterms : rgraph[i].split_terms<I4CIntegral>())
                {
                    auto rgroup = apply_ket_hrr(vterms, sints);
                    
                    if (rgroup.expansions() == 0)
                    {
                        for (const auto& tval : vterms)
                        {
                            rgroup.add(R4CDist(tval));
                        }
                    }
                    
                    rgraph.add(rgroup, i);
                    
                    cnt++;
                }
            }
        }
        
        if (cnt == 0) use_hrr = false;
    }
    
    rgraph.reduce();
    
    rgraph.sort<I4CIntegral>(true); 
}

void
EriDriver::apply_bra_vrr(R4Graph&       rgraph,
                         ST4CIntegrals& sints) const
{
    // special cases: single vertices without expansion terms
    
    for (const auto i : rgraph.orphans())
    {
        if (rgraph[i].empty() && (!rgraph[i].auxilary(1)))
        {
            const auto rgroup = apply_bra_vrr(rgraph[i].roots(), sints);

            rgraph.replace(rgroup, i);
        }
    }
    
    // loop over orphaned vertices
    
    bool use_hrr = true;
    
    while (use_hrr)
    {
        int cnt = 0;
        
        const auto vidx = rgraph.orphans();
        
        for (const auto i : vidx)
        {
            if (!rgraph[i].auxilary(1))
            {
                for (const auto& vterms : rgraph[i].split_terms<I4CIntegral>())
                {
                    auto rgroup = apply_bra_vrr(vterms, sints);
                    
                    if (rgroup.expansions() == 0)
                    {
                        for (const auto& tval : vterms)
                        {
                            rgroup.add(R4CDist(tval));
                        }
                    }
                    
                    rgraph.add(rgroup, i);
                    
                    cnt++;
                }
            }
        }
        
        if (cnt == 0) use_hrr = false;
    }
    
    rgraph.reduce();
    
    rgraph.sort<I4CIntegral>(true);
}

void
EriDriver::apply_ket_vrr(R4Graph&       rgraph,
                         ST4CIntegrals& sints) const
{
    // special cases: single vertices without expansion terms
    
    for (const auto i : rgraph.orphans())
    {
        if (rgraph[i].empty() && (!rgraph[i].auxilary(3)))
        {
            const auto rgroup = apply_ket_vrr(rgraph[i].roots(), sints);

            rgraph.replace(rgroup, i);
        }
    }
    
    // loop over orphaned vertices
    
    bool use_hrr = true;
    
    while (use_hrr)
    {
        int cnt = 0;
        
        const auto vidx = rgraph.orphans();
        
        for (const auto i : vidx)
        {
            if (!rgraph[i].auxilary(3))
            {
                for (const auto& vterms : rgraph[i].split_terms<I4CIntegral>())
                {
                    auto rgroup = apply_ket_vrr(vterms, sints);
                    
                    if (rgroup.expansions() == 0)
                    {
                        for (const auto& tval : vterms)
                        {
                            rgroup.add(R4CDist(tval));
                        }
                    }
                    
                    rgraph.add(rgroup, i);
                    
                    cnt++;
                }
            }
        }
        
        if (cnt == 0) use_hrr = false;
    }
    
    rgraph.reduce();
    
    rgraph.sort<I4CIntegral>(true);
}