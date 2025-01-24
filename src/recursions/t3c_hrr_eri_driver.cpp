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

#include "t3c_hrr_eri_driver.hpp"

#include "axes.hpp"
#include "t4c_utils.hpp"

T3CHrrElectronRepulsionDriver::T3CHrrElectronRepulsionDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T3CHrrElectronRepulsionDriver::is_electron_repulsion(const R3CTerm& rterm) const
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

std::optional<R3CDist>
T3CHrrElectronRepulsionDriver::ket_hrr(const R3CTerm& rterm,
                                       const char     axis) const
{
    if (!is_electron_repulsion(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 1))
    {
        R3CDist t3crt(rterm);
        
        // first recursion term
        
        auto x1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("DC", "cd", coord), Fraction(-1));
        
        t3crt.add(x1val);
        
        // second recursion term
        
        if (const auto r2val = tval->shift(axis, 1, 2))
        {
            t3crt.add(*r2val);
        }
        
        return t3crt;
    }
    else
    {
        return std::nullopt;
    }
}

R3CDist
T3CHrrElectronRepulsionDriver::apply_ket_hrr(const R3CTerm& rterm) const
{
    R3CDist t3crt;
    
    size_t nints = 3;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = ket_hrr(rterm, axis))
        {
            if (const auto nterms = trec->terms(); nterms < nints)
            {
                t3crt = *trec;
                
                nints = nterms;
            }
        }
    }
    
    return t3crt;
}

void
T3CHrrElectronRepulsionDriver::apply_recursion(R3CDist& rdist) const
{
    // vertical recursions on ket side center C
    
    apply_ket_hrr(rdist);
}

void
T3CHrrElectronRepulsionDriver::apply_ket_hrr(R3CDist& rdist) const
{
    if (!rdist.auxilary(1))
    {
        R3CDist new_dist(rdist.root());
            
        V3CTerms rec_terms;
            
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
            V3CTerms new_terms;
                
            for (size_t i = 0; i < rec_terms.size(); i++)
            {
                const auto cdist = apply_ket_hrr(rec_terms[i]);
                    
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

R3Group
T3CHrrElectronRepulsionDriver::create_recursion(const VT3CIntegrals& vints) const
{
    // create reference group
    
    R3Group r4group;
    
    for (const auto& tcomp : vints)
    {
        auto rdist = R3CDist(R3CTerm(tcomp));
        
        apply_recursion(rdist);
                
        r4group.add(rdist);
    }
    
    r4group.simplify();
    
    return r4group;
}

void
T3CHrrElectronRepulsionDriver::apply_recursion(R3Group& rgroup) const
{
    if (const auto nterms = rgroup.expansions(); nterms > 0)
    {
        R3Group mgroup;
        
        for (size_t i = 0; i < nterms; i++)
        {
            auto rdist = rgroup[i];
            
            apply_recursion(rdist);
            
            mgroup.add(rdist);
        }
        
        rgroup = mgroup;
    }
}
