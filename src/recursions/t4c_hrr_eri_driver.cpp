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

#include "t4c_hrr_eri_driver.hpp"

#include "axes.hpp"
#include "t4c_utils.hpp"

T4CHrrElectronRepulsionDriver::T4CHrrElectronRepulsionDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T4CHrrElectronRepulsionDriver::is_electron_repulsion(const R4CTerm& rterm) const
{
    if (!(rterm.prefixes()).empty())
    {
        return false;
    }
    
    if (rterm.integrand() != OperatorComponent("1/|r-r'|"))
    {
        return false;
    }
    else
    {
        return true;
    }
}

std::optional<R4CDist>
T4CHrrElectronRepulsionDriver::bra_hrr_a(const R4CTerm& rterm,
                                         const char     axis) const
{
    if (!is_electron_repulsion(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 0))
    {
        R4CDist t4crt(rterm);
        
        // first recursion term
        
        auto x1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("PA", "rpa", coord), Fraction(1));
        
        t4crt.add(x1val);
        
        // second recursion term
        
        if (const auto r2val = tval->shift_order(1))
        {
            auto x2val = *r2val;
            
            x2val.add(Factor("WP", "rwp", coord), Fraction(1));
            
            t4crt.add(x2val);
        }
        
        // third and fourth recursion terms
        
        if (const auto r3val = tval->shift(axis, -1, 0))
        {
            auto x3val = *r3val;
            
            const auto na = x1val[0][axis];
            
            x3val.add(Factor("1/eta", "fi_ab"), Fraction(na, 2));
            
            t4crt.add(x3val);
            
            if (const auto r4val = r3val->shift_order(1))
            {
                auto x4val = *r4val;
                
                x4val.add(Factor("1/eta^2", "fti_ab"), Fraction(-na, 2));
                
                t4crt.add(x4val);
            }
        }
        
        // fifth and sixth recursion terms
        
        if (const auto r5val = tval->shift(axis, -1, 1))
        {
            auto x5val = *r5val;
            
            const auto nb = x1val[1][axis];
            
            x5val.add(Factor("1/eta", "fi_ab"), Fraction(nb, 2));
            
            t4crt.add(x5val);
            
            if (const auto r6val = r5val->shift_order(1))
            {
                auto x6val = *r6val;
                
                x6val.add(Factor("1/eta^2", "fti_ab"), Fraction(-nb, 2));
                
                t4crt.add(x6val);
            }
        }
        
        // seventh recursion term
        
        if (const auto xval = tval->shift(axis, -1, 2))
        {
            if (const auto r7val = xval->shift_order(1))
            {
                auto x7val = *r7val;
               
                const auto nc = x1val[2][axis];
               
                x7val.add(Factor("1/(eta+nu)", "fi_abcd"), Fraction(nc, 2));
               
                t4crt.add(x7val);
            }
        }
        
        // eigth recursion term
        
        if (const auto xval = tval->shift(axis, -1, 3))
        {
            if (const auto r8val = xval->shift_order(1))
            {
                auto x8val = *r8val;
               
                const auto nd = x1val[3][axis];
               
                x8val.add(Factor("1/(eta+nu) ", "fi_abcd"), Fraction(nd, 2));
               
                t4crt.add(x8val);
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
T4CHrrElectronRepulsionDriver::ket_hrr_c(const R4CTerm& rterm,
                                         const char     axis) const
{
    if (!is_electron_repulsion(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 2))
    {
        R4CDist t4crt(rterm);
        
        // first recursion term
        
        auto x1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("QC", "rqc", coord), Fraction(1));
        
        t4crt.add(x1val);
        
        // second recursion term
        
        if (const auto r2val = tval->shift_order(1))
        {
            auto x2val = *r2val;
            
            x2val.add(Factor("WQ", "rwq", coord), Fraction(1));
            
            t4crt.add(x2val);
        }
        
        // third and fourth recursion terms
        
        if (const auto r3val = tval->shift(axis, -1, 2))
        {
            auto x3val = *r3val;
            
            const auto nc = x1val[2][axis];
            
            x3val.add(Factor("1/nu", "fi_cd"), Fraction(nc, 2));
            
            t4crt.add(x3val);
            
            if (const auto r4val = r3val->shift_order(1))
            {
                auto x4val = *r4val;
                
                x3val.add(Factor("1/nu^2", "fti_cd"), Fraction(-nc, 2));
                
                t4crt.add(x3val);
            }
        }
        
        // fifth and sixth recursion terms
        
        if (const auto r5val = tval->shift(axis, -1, 3))
        {
            auto x5val = *r5val;
            
            const auto nd = x1val[3][axis];
            
            x5val.add(Factor("1/nu", "fi_cd"), Fraction(nd, 2));
            
            t4crt.add(x5val);
            
            if (const auto r6val = r5val->shift_order(1))
            {
                auto x6val = *r6val;
                
                x6val.add(Factor("1/nu^2", "fti_cd"), Fraction(-nd, 2));
                
                t4crt.add(x6val);
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
T4CHrrElectronRepulsionDriver::apply_bra_hrr_a(const R4CTerm& rterm) const
{
    R4CDist t4crt;
    
    size_t nints = 9;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = bra_hrr_a(rterm, axis))
        {
            if (const auto nterms = trec->terms(); nterms < nints)
            {
                t4crt = *trec;
                
                nints = nterms;
            }
        }
    }
    
    return t4crt;
}


R4CDist
T4CHrrElectronRepulsionDriver::apply_ket_hrr_c(const R4CTerm& rterm) const
{
    R4CDist t4crt;
    
    size_t nints = 7;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = ket_hrr_c(rterm, axis))
        {
            if (const auto nterms = trec->terms(); nterms < nints)
            {
                t4crt = *trec;
                
                nints = nterms;
            }
        }
    }
    
    return t4crt;
}

void
T4CHrrElectronRepulsionDriver::apply_recursion(R4CDist& rdist) const
{
    // vertical recursions on bra side center A
    
    apply_bra_hrr_a(rdist);
    
    // vertical recursions on ket side center C
    
    apply_ket_hrr_c(rdist);
}

void
T4CHrrElectronRepulsionDriver::apply_bra_hrr_a(R4CDist& rdist) const
{
    if (!rdist.auxilary(0))
    {
        R4CDist new_dist(rdist.root());
            
        V4CTerms rec_terms;
            
        // set up initial terms for recursion expansion
            
        if (const auto nterms = rdist.terms(); nterms > 0)
        {
            for (size_t i = 0; i < nterms; i++)
            {
                if (const auto rterm = rdist[i]; is_electron_repulsion(rterm))
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
            if (const auto rterm = rdist.root(); is_electron_repulsion(rterm))
            {
                rec_terms.push_back(rterm);
            }
        }
            
        // apply recursion until only
                
        while (!rec_terms.empty())
        {
            V4CTerms new_terms;
                
            for (size_t i = 0; i < rec_terms.size(); i++)
            {
                const auto cdist = apply_bra_hrr_a(rec_terms[i]);
                    
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
T4CHrrElectronRepulsionDriver::apply_ket_hrr_c(R4CDist& rdist) const
{
    if (!rdist.auxilary(2))
    {
        R4CDist new_dist(rdist.root());
            
        V4CTerms rec_terms;
            
        // set up initial terms for recursion expansion
            
        if (const auto nterms = rdist.terms(); nterms > 0)
        {
            for (size_t i = 0; i < nterms; i++)
            {
                if (const auto rterm = rdist[i]; is_electron_repulsion(rterm))
                {
                    if (rterm.auxilary(2))
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
            if (const auto rterm = rdist.root(); is_electron_repulsion(rterm))
            {
                rec_terms.push_back(rterm);
            }
        }
            
        // apply recursion until only
                
        while (!rec_terms.empty())
        {
            V4CTerms new_terms;
                
            for (size_t i = 0; i < rec_terms.size(); i++)
            {
                const auto cdist = apply_ket_hrr_c(rec_terms[i]);
                    
                if (const auto nterms = cdist.terms(); nterms > 0)
                {
                    for (size_t j = 0; j < nterms; j++)
                    {
                        if (const auto rterm = cdist[j]; rterm.auxilary(2))
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

R4Group
T4CHrrElectronRepulsionDriver::create_recursion(const VT4CIntegrals& vints) const
{
    // create reference group
    
    R4Group r4group;
    
    for (const auto& tcomp : vints)
    {
        auto rdist = R4CDist(R4CTerm(tcomp));
        
        apply_recursion(rdist);
                
        r4group.add(rdist);
    }
    
    r4group.simplify();
    
    return r4group;
}

void
T4CHrrElectronRepulsionDriver::apply_recursion(R4Group& rgroup) const
{
    if (const auto nterms = rgroup.expansions(); nterms > 0)
    {
        R4Group mgroup;
        
        for (size_t i = 0; i < nterms; i++)
        {
            auto rdist = rgroup[i];
            
            apply_recursion(rdist);
            
            mgroup.add(rdist);
        }
        
        rgroup = mgroup;
    }
}
