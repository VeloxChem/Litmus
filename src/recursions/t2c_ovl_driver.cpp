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

#include "axes.hpp"

T2COverlapDriver::T2COverlapDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

std::optional<R2CDist>
T2COverlapDriver::bra_vrr(const R2CTerm& rterm,
                          const char     axis) const
{
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
T2COverlapDriver::apply_bra_vrr(const R2CTerm&       rterm,
                                      ST2CIntegrals& sints) const
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
    
    const auto vints = t2crt.unique_integrals();
        
    sints.insert(vints.cbegin(), vints.cend());

    return t2crt;
}

R2CDist
T2COverlapDriver::apply_ket_vrr(const R2CTerm&       rterm,
                                      ST2CIntegrals& sints) const
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
    
    const auto vints = t2crt.unique_integrals();
        
    sints.insert(vints.cbegin(), vints.cend());
    
    return t2crt;
}

R2Group
T2COverlapDriver::apply_bra_vrr(const V2CTerms&      rterms,
                                      ST2CIntegrals& sints) const
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
T2COverlapDriver::apply_ket_vrr(const V2CTerms&      rterms,
                                      ST2CIntegrals& sints) const
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
