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


#include "cold_kin_driver.hpp"

#include "axes.hpp"
#include "cold_ovl_driver.hpp"

ColdKineticEnergyDriver::ColdKineticEnergyDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
ColdKineticEnergyDriver::is_kinetic_energy(const R2CTerm& rterm) const
{
    if (!(rterm.prefixes()).empty())
    {
        return false;
    }
    
    if (rterm.integrand() != OperatorComponent("T"))
    {
        return false;
    }
    else
    {
        return true;
    }
}

std::optional<R2CDist>
ColdKineticEnergyDriver::bra_vrr(const R2CTerm& rterm,
                                 const char     axis) const
{
    if (!is_kinetic_energy(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 0))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        auto x1val = rterm.replace(OperatorComponent("1"));
            
        x1val.add(Factor("N", "n"), Fraction(1));
        
        x1val.add(Factor("M", "m"), Fraction(1));
        
        x1val.add(Factor("T", "t"), Fraction(2));
            
        t2crt.add(x1val);
        
        // second recursion term
        
        auto x2val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x2val.add(Factor("AB", "rab", coord), Fraction(-1));
        
        x2val.add(Factor("M", "m"), Fraction(1));
        
        x2val.add(Factor("T", "t"), Fraction(1));
        
        t2crt.add(x2val);
        
        // third recursion term
        
        if (const auto r3val = tval->shift(axis, -1, 0))
        {
            auto x3val = *r3val;
            
            const auto na = x2val[0][axis];
            
            x3val.add(Factor("T", "t"), Fraction(na, 2));
            
            t2crt.add(x3val);
        }
        
        // fourth recursion term
        
        if (const auto r4val = tval->shift(axis, -1, 1))
        {
            auto x4val = *r4val;
            
            const auto nb = x2val[1][axis];
            
            x4val.add(Factor("T", "t"), Fraction(nb, 2));
            
            t2crt.add(x4val);
        }
        
        // fifth recursion term
        
        if (const auto r5val = tval->shift(axis, -1, 0))
        {
            auto x5val = *r5val;
            
            const auto na = x2val[0][axis];
            
            x5val = x5val.replace(OperatorComponent("1"));
            
            x5val.add(Factor("M", "m"), Fraction(1));
            
            x5val.add(Factor("T", "t"), Fraction(-na));
            
            t2crt.add(x5val);
        }
        
        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}

std::optional<R2CDist>
ColdKineticEnergyDriver::ket_vrr(const R2CTerm& rterm,
                                 const char     axis) const
{
    if (!is_kinetic_energy(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 1))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        auto x1val = rterm.replace(OperatorComponent("1"));
            
        x1val.add(Factor("N", "n"), Fraction(1));
        
        x1val.add(Factor("M", "m"), Fraction(1));
        
        x1val.add(Factor("T", "t"), Fraction(2));
        
        t2crt.add(x1val);
        
        // second recursion term
        
        auto x2val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x2val.add(Factor("AB", "rab", coord), Fraction(1));
        
        x2val.add(Factor("N", "n"), Fraction(1));
        
        x2val.add(Factor("T", "t"), Fraction(1));
        
        t2crt.add(x2val);
        
        // third recursion term
        
        if (const auto r3val = tval->shift(axis, -1, 1))
        {
            auto x3val = *r3val;
            
            const auto nb = x2val[1][axis];
            
            x3val.add(Factor("T", "t"), Fraction(nb, 2));
            
            t2crt.add(x3val);
        }
        
        // fourth recursion term
        
        if (const auto r4val = tval->shift(axis, -1, 1))
        {
            auto x4val = *r4val;
            
            const auto nb = x2val[1][axis];
            
            x4val = x4val.replace(OperatorComponent("1"));
            
            x4val.add(Factor("N", "n"), Fraction(1));
            
            x4val.add(Factor("T", "t"), Fraction(-nb));
            
            t2crt.add(x4val);
        }
        
        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}

R2CDist
ColdKineticEnergyDriver::aux_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    if (!is_kinetic_energy(rterm)) return t2crt;
    
    if (rterm.auxilary(0) && rterm.auxilary(1))
    {
        // first recursion term
        
        auto x1val = rterm.replace(OperatorComponent("1"));
        
        x1val.add(Factor("N", "n"), Fraction(1));
        
        x1val.add(Factor("M", "m"), Fraction(1));
        
        x1val.add(Factor("T", "t"), Fraction(3));
            
        t2crt.add(x1val);
        
        // second recursion term
        
        x1val.add(Factor("N", "n"), Fraction(1));
        
        x1val.add(Factor("M", "m"), Fraction(1));
        
        x1val.add(Factor("T", "t"), Fraction(2, 3));
        
        x1val.add(Factor("AB2", "r2ab"), Fraction(-1));
        
        t2crt.add(x1val);
    }
    
    return t2crt;
}

R2CDist
ColdKineticEnergyDriver::apply_bra_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    size_t nints = 6;
    
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
ColdKineticEnergyDriver::apply_ket_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    size_t nints = 5;
    
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
ColdKineticEnergyDriver::apply_recursion(R2CDist& rdist) const
{
    // vertical recursions on bra side
    
    apply_bra_vrr(rdist);
    
    // vertical recursions on ket side
    
    apply_ket_vrr(rdist);
    
    // auxilary recursion
    
    apply_aux_vrr(rdist);
}

void
ColdKineticEnergyDriver::apply_bra_vrr(R2CDist& rdist) const
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
                if (const auto rterm = rdist[i]; is_kinetic_energy(rterm))
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
            if (const auto rterm = rdist.root(); is_kinetic_energy(rterm))
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
                        if (const auto rterm = cdist[j]; (rterm.auxilary(0) || (!is_kinetic_energy(rterm))))
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
ColdKineticEnergyDriver::apply_ket_vrr(R2CDist& rdist) const
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
                if (const auto rterm = rdist[i]; is_kinetic_energy(rterm))
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
            if (const auto rterm = rdist.root(); is_kinetic_energy(rterm))
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
                        if (const auto rterm = cdist[j]; (rterm.auxilary(1) || (!is_kinetic_energy(rterm))))
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
ColdKineticEnergyDriver::apply_aux_vrr(R2CDist& rdist) const
{
    R2CDist new_dist(rdist.root());
    
    // update recursion distribution
    
    if (const auto nterms = rdist.terms(); nterms > 0)
    {
        for (size_t i = 0; i < nterms; i++)
        {
            if (const auto rterm = rdist[i]; is_kinetic_energy(rterm))
            {
                const auto cdist = aux_vrr(rterm);
                
                if (const auto nterms = cdist.terms(); nterms > 0)
                {
                    for (size_t i = 0; i < nterms; i++ )
                    {
                        new_dist.add(cdist[i]);
                    }
                }
                else
                {
                    new_dist.add(rterm);
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
        if (const auto rterm = rdist.root(); is_kinetic_energy(rterm))
        {
            const auto cdist = aux_vrr(rterm);
            
            if (const auto nterms = cdist.terms(); nterms > 0)
            {
                for (size_t i = 0; i < nterms; i++ )
                {
                    new_dist.add(cdist[i]);
                }
            }
        }
    }
   
    rdist = new_dist;
}

R2Group
ColdKineticEnergyDriver::create_recursion(const VT2CIntegrals& vints) const
{
    // create overlap integrals driver
    
    ColdOverlapDriver ovl_drv;
    
    // create reference group
    
    R2Group r2group;
    
    for (const auto& tcomp : vints)
    {
        auto rdist = R2CDist(R2CTerm(tcomp));
        
        // apply kinetic energy recursion
        
        apply_recursion(rdist);
        
        // apply overlap recursion
        
        ovl_drv.apply_recursion(rdist);
        
        r2group.add(rdist);
    }
    
    r2group.simplify();
    
    return r2group;
}

void
ColdKineticEnergyDriver::apply_recursion(R2Group& rgroup) const
{
    if (const auto nterms = rgroup.expansions(); nterms > 0)
    {
        // create overlap integrals driver
        
        ColdOverlapDriver ovl_drv;
        
        R2Group mgroup;
        
        for (size_t i = 0; i < nterms; i++)
        {
            auto rdist = rgroup[i];
            
            apply_recursion(rdist);
            
            // apply overlap recursion
            
            ovl_drv.apply_recursion(rdist);
            
            mgroup.add(rdist);
        }
        
        rgroup = mgroup;
    }
}
