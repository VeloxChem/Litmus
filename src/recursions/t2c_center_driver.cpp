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

#include "t2c_center_driver.hpp"

#include "axes.hpp"

T2CCenterDriver::T2CCenterDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T2CCenterDriver::is_auxilary(const R2CTerm& rterm,
                             const int      index) const
{
    if (const auto nprefixes = rterm.prefixes().size(); index > nprefixes)
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::optional<R2CDist>
T2CCenterDriver::bra_ket_vrr(const R2CTerm& rterm,
                             const char     axis,
                             const int      index) const
{
    if (is_auxilary(rterm, index)) return std::nullopt;
    
    if (const auto tval = rterm.shift_prefix(axis, -1, index, true))
    {
        R2CDist t2crt(rterm);
        
        if (const auto r1val = tval->shift(axis, 1, index))
        {
            auto x1val = *r1val;
            
            if (index == 0)
            {
                x1val.add(Factor("b_e", "tbe"), Fraction(2));
            }
            else
            {
                x1val.add(Factor("k_e", "tke"), Fraction(2));
            }
            
            t2crt.add(x1val); 
        }
        
        if (const auto r2val = tval->shift(axis, -1, index))
        {
            auto x2val = *r2val;
            
            x2val.scale(Fraction(-(*tval)[index][axis]));
            
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
T2CCenterDriver::apply_bra_ket_vrr(const R2CTerm& rterm,
                                   const int      index) const
{
    R2CDist t2crt;
    
    size_t nints = 3;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = bra_ket_vrr(rterm, axis, index))
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
T2CCenterDriver::apply_recursion(R2CDist& rdist) const
{
    // vertical recursions on ket side
    
    apply_bra_ket_vrr(rdist, 1);
    
    // vertical recursions on bra side
    
    apply_bra_ket_vrr(rdist, 0);
}

void
T2CCenterDriver::apply_bra_ket_vrr(      R2CDist& rdist,
                                   const int      index) const
{
    if (!is_auxilary(rdist.root(), index))
    {
        R2CDist new_dist(rdist.root());
        
        V2CTerms rec_terms;
        
        // set up initial terms for recursion expansion
        
        if (const auto nterms = rdist.terms(); nterms > 0)
        {
            for (size_t i = 0; i < nterms; i++)
            {
                if (const auto rterm = rdist[i]; is_auxilary(rterm, index))
                {
                    new_dist.add(rterm);
                }
                else
                {
                    rec_terms.push_back(rterm);
                }
            }
        }
        else
        {
            if (const auto rterm = rdist.root(); !is_auxilary(rterm, index))
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
                const auto cdist = apply_bra_ket_vrr(rec_terms[i], index);
                
                if (const auto nterms = cdist.terms(); nterms > 0)
                {
                    for (size_t j = 0; j < nterms; j++)
                    {
                        if (const auto rterm = cdist[j]; is_auxilary(rterm, index))
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
T2CCenterDriver::create_recursion(const VT2CIntegrals& vints) const
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
