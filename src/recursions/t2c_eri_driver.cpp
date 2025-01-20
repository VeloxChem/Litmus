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

#include "t2c_eri_driver.hpp"

#include "axes.hpp"

T2CElectronRepulsionDriver::T2CElectronRepulsionDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T2CElectronRepulsionDriver::is_electron_repulsion(const R2CTerm& rterm) const
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

std::optional<R2CDist>
T2CElectronRepulsionDriver::bra_vrr(const R2CTerm& rterm,
                                    const char     axis) const
{
    if (!is_electron_repulsion(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 0))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        if (const auto r1val = tval->shift_order(1))
        {
            auto x1val = *r1val;
            
            const auto coord = _rxyz[axes::to_index(axis)];
            
            x1val.add(Factor("PA", "pa", coord), Fraction(1));
            
            t2crt.add(x1val);
            
            // second recursion term
            
            if (const auto r2val = tval->shift(axis, -1, 0))
            {
                auto x2val = *r2val;
                
                const auto na = x1val[0][axis];
                
                x2val.add(Factor("1/b_e", "fbe"), Fraction(na));
                
                t2crt.add(x2val);
                
                // third recursion term
                
                if (const auto r3val = r2val->shift_order(1))
                {
                    auto x3val = *r3val;
                    
                    x3val.add(Factor("zeta/b_e^2", "fz_be"), Fraction(-na));
                    
                    t2crt.add(x3val);
                }
            }
            
            // fourth recursion term
            
            if (const auto r4val = tval->shift(axis, -1, 1))
            {
                if (const auto r4mval = r4val->shift_order(1))
                {
                    auto x4val = *r4mval;
                    
                    const auto nb = x1val[1][axis];
                    
                    x4val.add(Factor("1/eta", "fe"), Fraction(nb));
                    
                    t2crt.add(x4val);
                }
            }
        }
        
        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}

std::optional<R2CDist>
T2CElectronRepulsionDriver::ket_vrr(const R2CTerm& rterm,
                          const char     axis) const
{
    if (!is_electron_repulsion(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 1))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        if (const auto r1val = tval->shift_order(1))
        {
            auto x1val = *r1val;
            
            const auto coord = _rxyz[axes::to_index(axis)];
            
            x1val.add(Factor("PB", "pb", coord), Fraction(1));
            
            t2crt.add(x1val);
            
            // second recursion term
            
            if (const auto r2val = tval->shift(axis, -1, 1))
            {
                auto x2val = *r2val;
                
                const auto nb = x1val[1][axis];
                
                x2val.add(Factor("1/k_e", "fke"), Fraction(nb));
                
                t2crt.add(x2val);
                
                // third recursion term
                
                if (const auto r3val = r2val->shift_order(1))
                {
                    auto x3val = *r3val;
                    
                    x3val.add(Factor("zeta/k_e^2", "fz_ke"), Fraction(-nb));
                    
                    t2crt.add(x3val);
                }
            }
        }
        
        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}

R2CDist
T2CElectronRepulsionDriver::apply_bra_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    size_t nints = 5;
    
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
T2CElectronRepulsionDriver::apply_ket_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    size_t nints = 4;
    
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

void
T2CElectronRepulsionDriver::apply_recursion(R2CDist& rdist) const
{

    // vertical recursions on bra side
    
    apply_bra_vrr(rdist);
    
    // vertical recursions on ket side
    
    apply_ket_vrr(rdist);
}

void
T2CElectronRepulsionDriver::apply_bra_vrr(R2CDist& rdist) const
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
T2CElectronRepulsionDriver::apply_ket_vrr(R2CDist& rdist) const
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
                if (const auto rterm = rdist[i]; is_electron_repulsion(rterm))
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
            if (const auto rterm = rdist.root(); is_electron_repulsion(rterm))
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

R2Group
T2CElectronRepulsionDriver::create_recursion(const VT2CIntegrals& vints) const
{
    // create reference group

    R2Group r2group;
    
    for (const auto& tcomp : vints)
    {
        auto rdist = R2CDist(R2CTerm(tcomp));
        
        apply_recursion(rdist);
                
        r2group.add(rdist);
    }
    
    r2group.simplify();
    
    return r2group;
}

void
T2CElectronRepulsionDriver::apply_recursion(R2Group& rgroup) const
{
    if (const auto nterms = rgroup.expansions(); nterms > 0)
    {
        R2Group mgroup;
        
        for (size_t i = 0; i < nterms; i++)
        {
            auto rdist = rgroup[i];
            
            apply_recursion(rdist);
            
            mgroup.add(rdist);
        }
        
        rgroup = mgroup;
    }
}
