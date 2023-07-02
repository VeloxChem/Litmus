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

#include "t2c_npot_geom_driver.hpp"

#include <iostream>

#include "axes.hpp"
#include "t2c_npot_driver.hpp"

T2CNuclearPotentialGeometryDriver::T2CNuclearPotentialGeometryDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T2CNuclearPotentialGeometryDriver::is_nuclear_potential_geom(const R2CTerm& rterm) const
{
    if (!(rterm.prefixes()).empty())
    {
        return false;
    }
    
    if (const auto integrand = rterm.integrand(); integrand.name() != "AG")
    {
        return false;
    }
    else
    {
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

std::optional<R2CDist>
T2CNuclearPotentialGeometryDriver::bra_vrr(const R2CTerm& rterm,
                                           const char     axis) const
{
    if (!is_nuclear_potential_geom(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 0))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        auto x1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("PA", "rpa", coord), Fraction(1));
        
        t2crt.add(x1val);
        
        // second recursion term
        
        if (const auto r2val = tval->shift_order(1))
        {
            auto x2val = *r2val;
            
            x2val.add(Factor("PC", "rpc", coord), Fraction(-1));
            
            t2crt.add(x2val);
        }
        
        // third and fourth recursion terms
        
        if (const auto r3val = tval->shift(axis, -1, 0))
        {
            auto x3val = *r3val;
            
            const auto na = x1val[0][axis];
            
            x3val.add(Factor("1/eta", "fe"), Fraction(na, 2));
            
            t2crt.add(x3val);
            
            if (const auto r4val = r3val->shift_order(1))
            {
                auto x4val = *r4val;
                
                x4val.add(Factor("1/eta", "fe"), Fraction(-na, 2));
                
                t2crt.add(x4val);
            }
        }
        
        // fifth and sixth recursion terms
        
        if (const auto r5val = tval->shift(axis, -1, 1))
        {
            auto x5val = *r5val;
            
            const auto nb = x1val[1][axis];
            
            x5val.add(Factor("1/eta", "fe"), Fraction(nb, 2));
            
            t2crt.add(x5val);
            
            if (const auto r6val = r5val->shift_order(1))
            {
                auto x6val = *r6val;
                
                x6val.add(Factor("1/eta", "fe"), Fraction(-nb, 2));
                
                t2crt.add(x6val);
            }
        }
        
        // seventh recursion term
        
        if (const auto r7val = tval->shift_operator(axis, -1))
        {
            auto x7val = *(r7val->shift_order(1));
            
            const auto nc = tval->integrand()[axis];
            
            if (x7val.integrand().shape() == TensorComponent(0, 0, 0))
            {
                x7val = x7val.replace(OperatorComponent("A"));
            }
            
            x7val.scale(Fraction(nc));
            
            t2crt.add(x7val);
        }
        
        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}

std::optional<R2CDist>
T2CNuclearPotentialGeometryDriver::ket_vrr(const R2CTerm& rterm,
                                           const char     axis) const
{
    if (!is_nuclear_potential_geom(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 1))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        auto x1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("PB", "rpb", coord), Fraction(1));
        
        t2crt.add(x1val);
        
        // second recursion term
        
        if (const auto r2val = tval->shift_order(1))
        {
            auto x2val = *r2val;
            
            x2val.add(Factor("PC", "rpc", coord), Fraction(-1));
            
            t2crt.add(x2val);
        }
        
        // third and fourth recursion terms
        
        if (const auto r3val = tval->shift(axis, -1, 1))
        {
            auto x3val = *r3val;
            
            const auto nb = x1val[1][axis];
            
            x3val.add(Factor("1/eta", "fe"), Fraction(nb, 2));
            
            t2crt.add(x3val);
            
            if (const auto r4val = r3val->shift_order(1))
            {
                auto x4val = *r4val;
                
                x4val.add(Factor("1/eta", "fe"), Fraction(-nb, 2));
                
                t2crt.add(x4val);
            }
        }
        
        // fifth recursion term
        
        if (const auto r5val = tval->shift_operator(axis, -1))
        {
            auto x5val = *(r5val->shift_order(1));
            
            const auto nc = tval->integrand()[axis];
            
            if (x5val.integrand().shape() == TensorComponent(0, 0, 0))
            {
                x5val = x5val.replace(OperatorComponent("A"));
            }
            
            x5val.scale(Fraction(nc));
            
            t2crt.add(x5val);
        }
        
        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}

R2CDist
T2CNuclearPotentialGeometryDriver::apply_bra_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    size_t nints = 8;
    
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
T2CNuclearPotentialGeometryDriver::apply_ket_vrr(const R2CTerm& rterm) const
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

void
T2CNuclearPotentialGeometryDriver::apply_recursion(R2CDist& rdist) const
{
    // vertical recursions on bra side
    
    apply_bra_vrr(rdist);
    
    // vertical recursions on ket side
    
    apply_ket_vrr(rdist);
}

void
T2CNuclearPotentialGeometryDriver::apply_bra_vrr(R2CDist& rdist) const
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
                if (const auto rterm = rdist[i]; is_nuclear_potential_geom(rterm))
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
            if (const auto rterm = rdist.root(); is_nuclear_potential_geom(rterm))
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
                        if (const auto rterm = cdist[j]; (rterm.auxilary(0) || (!is_nuclear_potential_geom(rterm))))
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
T2CNuclearPotentialGeometryDriver::apply_ket_vrr(R2CDist& rdist) const
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
                if (const auto rterm = rdist[i]; is_nuclear_potential_geom(rterm))
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
            if (const auto rterm = rdist.root(); is_nuclear_potential_geom(rterm))
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
                        if (const auto rterm = cdist[j]; (rterm.auxilary(1) || (!is_nuclear_potential_geom(rterm))))
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
T2CNuclearPotentialGeometryDriver::create_recursion(const VT2CIntegrals& vints) const
{
    // create nuclear potential integrals driver
    
    T2CNuclearPotentialDriver npot_drv;
    
    R2Group r2group;
    
    for (const auto& tcomp : vints)
    {
        auto rdist = R2CDist(R2CTerm(tcomp));
        
        apply_recursion(rdist);
        
        // apply nuclear potential recursion
        
        npot_drv.apply_recursion(rdist);
        
//                std::cout << "*** RECURSION FOR INTEGRAL COMPONENT: " << rdist.root().label() << std::endl;
//        
//                std::cout << " NUMBER OF TERMS:" << rdist.terms() << std::endl;
//                
//                for (size_t i = 0; i < rdist.terms(); i++)
//               {
//                    std::cout << " RECURSION TERM (" << i << "): " << rdist[i].integrand().name() << " Factors: " ;
//        
//                  for (const auto& fact : rdist[i].factors())
//                    {
//                        std::cout << fact.label() << " , ";
//                    }
//        
//                    std::cout << std::endl;
//                }
//    
//               std::cout << std::endl << std::endl;
        
        r2group.add(rdist);
    }
    
    r2group.simplify();
    
    return r2group;
}
