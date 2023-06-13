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

#include "t2c_ovl_driver.hpp"

#include <iostream>

#include "axes.hpp"

T2COverlapDriver::T2COverlapDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T2COverlapDriver::is_overlap(const R2CTerm& rterm) const
{
    if (!(rterm.prefixes()).empty())
    {
        return false;
    }
    
    if (rterm.integrand() != OperatorComponent("1"))
    {
        return false;
    }
    else
    {
        return true;
    }
}

std::optional<R2CDist>
T2COverlapDriver::bra_vrr(const R2CTerm& rterm,
                          const char     axis) const
{
    if (!is_overlap(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 0))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        auto r1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        r1val.add(Factor("PA", "rpa", coord), Fraction(1));
        
        t2crt.add(r1val);
        
        // second recursion term
        
        if (const auto r2val = tval->shift(axis, -1, 0))
        {
            auto x2val = *r2val;
            
            const auto na = r1val[0][axis];
            
            x2val.add(Factor("1/eta", "fe"), Fraction(na, 2));
            
            t2crt.add(x2val);
        }
        
        // third recursion term
        
        if (const auto r3val = tval->shift(axis, -1, 1))
        {
            auto x3val = *r3val;
            
            const auto nb = r1val[1][axis];
            
            x3val.add(Factor("1/eta", "fe"), Fraction(nb, 2));
            
            t2crt.add(x3val);
        }
        
        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}

std::optional<R2CDist>
T2COverlapDriver::ket_vrr(const R2CTerm& rterm,
                          const char     axis) const
{
    if (!is_overlap(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 1))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        auto r1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        r1val.add(Factor("PB", "rpb", coord), Fraction(1));
        
        t2crt.add(r1val);
        
        // second recursion term
        
        if (const auto r2val = tval->shift(axis, -1, 1))
        {
            auto x2val = *r2val;
            
            const auto nb = r1val[1][axis];
            
            x2val.add(Factor("1/eta", "fe"), Fraction(nb, 2));
            
            t2crt.add(x2val);
        }
        
        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}

R2CDist
T2COverlapDriver::apply_bra_vrr(const R2CTerm& rterm,
                                      R2CMap&  sints) const
{
    R2CDist t2crt;
    
    size_t nints = 3;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = bra_vrr(rterm, axis))
        {
            const auto nterms = trec->count_new_integrals(sints);
            
            if (nterms < nints)
            {
                t2crt = *trec;
                
                nints = nterms;
            }
        }
    }
        
    sints.add(t2crt.unique_integrals());

    return t2crt;
}

R2CDist
T2COverlapDriver::apply_bra_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    size_t nints = 4;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = bra_vrr(rterm, axis))
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

R2CDist
T2COverlapDriver::apply_ket_vrr(const R2CTerm& rterm,
                                      R2CMap&  sints) const
{
    R2CDist t2crt;
    
    size_t nints = 2;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = ket_vrr(rterm, axis))
        {
            const auto nterms = trec->count_new_integrals(sints);
            
            if (nterms < nints)
            {
                t2crt = *trec;
                
                nints = nterms;
            }
        }
    }
        
    sints.add(t2crt.unique_integrals());
    
    return t2crt;
}

R2CDist
T2COverlapDriver::apply_ket_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    size_t nints = 3;
    
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

R2Group
T2COverlapDriver::apply_bra_vrr(const V2CTerms& rterms,
                                      R2CMap&   sints) const
{
    R2Group rgroup;
    
    for (const auto& rterm : rterms)
    {
        if (const auto tval = apply_bra_vrr(rterm, sints); tval.terms() > 0)
        {
            rgroup.add(tval);
        }
    }
    
    return rgroup;
}

R2Group
T2COverlapDriver::apply_ket_vrr(const V2CTerms& rterms,
                                      R2CMap&   sints) const
{
    R2Group rgroup;
    
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
T2COverlapDriver::apply_bra_vrr(R2GroupContainer& rgroups,
                                R2CMap&           sints) const
{
    // special cases: recursion groups without expansion terms
    
    auto ngroups = rgroups.recursion_groups();
    
    for (size_t i = 0; i < ngroups; i++)
    {
        if (rgroups[i].empty() && (!rgroups[i].auxilary(0)))
        {
            const auto rgroup = apply_bra_vrr(rgroups[i].roots(), sints);

            rgroups.replace(rgroup, i);
        }
    }
    
    // loop over unprocessed recursion groups
    
    size_t gstart = 0; size_t gend = ngroups;
    
    while (gend != gstart)
    {
        for (size_t i = gstart; i < gend; i++)
        {
            if (!rgroups[i].auxilary(0))
            {
                for (const auto& vterms : rgroups[i].split_terms<I2CIntegral>())
                {
                    auto rgroup = apply_bra_vrr(vterms, sints);
                    
                    if (rgroup.expansions() == 0)
                    {
                        for (const auto& tval : vterms)
                        {
                            rgroup.add(R2CDist(tval));
                        }
                    }
                    
                    rgroups.add(rgroup);
                }
            }
        }
        
        gstart = gend;
        
        gend = rgroups.recursion_groups();
    }
}

void
T2COverlapDriver::apply_ket_vrr(R2GroupContainer& rgroups,
                                R2CMap&           sints) const
{
    // special cases: recursion groups without expansion terms
    
    auto ngroups = rgroups.recursion_groups();
    
    for (size_t i = 0; i < ngroups; i++)
    {
        if (rgroups[i].empty() && (!rgroups[i].auxilary(1)))
        {
            const auto rgroup = apply_ket_vrr(rgroups[i].roots(), sints);

            rgroups.replace(rgroup, i);
        }
    }
    
    // loop over unprocessed recursion groups
    
    size_t gstart = 0; size_t gend = ngroups;
    
    while (gend != gstart)
    {
        for (size_t i = gstart; i < gend; i++)
        {
            if (!rgroups[i].auxilary(1))
            {
                for (const auto& vterms : rgroups[i].split_terms<I2CIntegral>())
                {
                    auto rgroup = apply_ket_vrr(vterms, sints);
                    
                    if (rgroup.expansions() == 0)
                    {
                        for (const auto& tval : vterms)
                        {
                            rgroup.add(R2CDist(tval));
                        }
                    }
                    
                    rgroups.add(rgroup);
                }
            }
        }
        
        gstart = gend;
        
        gend = rgroups.recursion_groups();
    }
}

void
T2COverlapDriver::apply_recursion(R2GroupContainer& rgroups,
                                  R2CMap&           sints) const
{
    // vertical recursions
    
    apply_bra_vrr(rgroups, sints);
    
    apply_ket_vrr(rgroups, sints);
}

void
T2COverlapDriver::apply_recursion(R2CDist& rdist) const
{
    // vertical recursions on bra side
    
    apply_bra_vrr(rdist);
    
    // vertical recursions on ket side
    
    apply_ket_vrr(rdist);
}

void
T2COverlapDriver::apply_bra_vrr(R2CDist& rdist) const
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
                if (const auto rterm = rdist[i]; is_overlap(rterm))
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
            if (const auto rterm = rdist.root(); is_overlap(rterm))
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
                const auto cdist = apply_bra_vrr(rec_terms[i]);
                
                if (const auto nterms = cdist.terms(); nterms > 0)
                {
                    for (size_t j = 0; j < nterms; j++)
                    {
                        if (const auto rterm = cdist[j]; rterm.auxilary(0))
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
T2COverlapDriver::apply_ket_vrr(R2CDist& rdist) const
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
                if (const auto rterm = rdist[i]; is_overlap(rterm))
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
            if (const auto rterm = rdist.root(); is_overlap(rterm))
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
                        if (const auto rterm = cdist[j]; rterm.auxilary(1))
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

R2GroupContainer
T2COverlapDriver::create_container(const int anga,
                                   const int angb) const
{
    // reference integral
    
    const auto operi = Operator("1");
    
    const auto bra = I1CPair("GA", anga);
    
    const auto ket = I1CPair("GB", angb);
    
    const auto refint = I2CIntegral(bra, ket, operi);
    
    // create refrence integral components
    
    auto refcomps = refint.components<T1CPair, T1CPair>();
    
    // create reference group
    
    R2Group r2group;
    
    for (const auto& tcomp : refcomps)
    {
        r2group.add(R2CDist(R2CTerm(tcomp)));
    }
    
    // apply Obara-Saika recursion
    
    R2CMap sints;
    
    R2GroupContainer rcont({r2group,});
    
    apply_recursion(rcont, sints);
    
    std::cout << "Create recursion groups container for " << refint.label() << " : " << sints.size() << std::endl;
    
    return rcont;
}

V2GroupContainers
T2COverlapDriver::create_containers(const int mang) const
{
    V2GroupContainers vconts;
    
    for (int i = 0; i <= mang; i++)
    {
        for (int j = 0; j <= mang; j++)
        {
            vconts.push_back(create_container(i, j));
        }
   }
    
    return vconts;
}

R2Group
T2COverlapDriver::create_recursion(const VT2CIntegrals& vints) const
{
    // create reference group
    
    R2Group r2group;
    
    for (const auto& tcomp : vints)
    {
        auto rdist = R2CDist(R2CTerm(tcomp));
        
        apply_recursion(rdist);
        
//        std::cout << "*** RECURSION FOR INTEGRAL COMPONENT: " << rdist.root().label() << std::endl;
//
//        std::cout << " NUMBER OF TERMS:" << rdist.terms() << std::endl;
//
//        for (size_t i = 0; i < rdist.terms(); i++)
//        {
//            std::cout << " RECURSION TERM (" << i << "): " << rdist[i].label() << std::endl;
//        }
//
//        std::cout << std::endl << std::endl;
        
        r2group.add(rdist);
    }
    
    r2group.simplify(); 
    
    return r2group;
}
