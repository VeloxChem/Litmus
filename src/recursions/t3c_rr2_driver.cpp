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

#include "t3c_rr2_driver.hpp"

#include "axes.hpp"
#include "t2c_utils.hpp"
#include "t3c_r2_driver.hpp"

T3CRR2Driver::T3CRR2Driver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T3CRR2Driver::is_rr2(const R2CTerm& rterm) const
{
    if (!(rterm.prefixes()).empty())
    {
        return false;
    }
    
    if (rterm.integrand().name() != "GR.R2(r)")
    {
        return false;
    }
    else
    {
        return true;
    }
}

std::optional<R2CDist>
T3CRR2Driver::aux_vrr(const R2CTerm& rterm,
                      const char     axis) const
{
    if (!is_rr2(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift_operator(axis, -1))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        auto r1val = tval->replace(OperatorComponent("GR2(r)"));
        
        auto x1val = r1val;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("GC", "gc", coord), Fraction(1));
                
        t2crt.add(x1val);
        
        // second recursion term
        
        if (const auto r2val = r1val.shift(axis, -1, 0))
        {
            auto x2val = *r2val;
            
            const auto na = x1val[0][axis];
            
            x2val.add(Factor("1/geta", "gfe"), Fraction(na));
            
            t2crt.add(x2val);
        }
        
        // third recursion term
        
        if (const auto r3val = r1val.shift(axis, -1, 1))
        {
            auto x3val = *r3val;
            
            const auto nb = x1val[1][axis];
            
            x3val.add(Factor("1/geta", "gfe"), Fraction(nb));
            
            t2crt.add(x3val);
        }
        
        // fourth recursion term
        
        auto r4val = tval->replace(OperatorComponent("G(r)"));
        
        auto x4val = r4val;
    
        x4val.add(Factor("GC", "gc", coord), Fraction(1));
        
        x4val.add(Factor("1/geta", "gfe"), Fraction(1));
                
        t2crt.add(x4val);
        
        // fifth recursion term
        
        if (const auto r5val = r4val.shift(axis, -1, 0))
        {
            auto x5val = *r5val;
            
            const auto na = x4val[0][axis];
            
            x5val.add(Factor("1/geta2", "gfe2"), Fraction(na));
            
            t2crt.add(x5val);
        }
        
        // sixth recursion term
        
        if (const auto r6val = r4val.shift(axis, -1, 1))
        {
            auto x6val = *r6val;
            
            const auto nb = x4val[1][axis];
            
            x6val.add(Factor("1/geta2", "gfe2"), Fraction(nb));
            
            t2crt.add(x6val);
        }
        
        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}

R2CDist
T3CRR2Driver::apply_aux_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    size_t nints = 9;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = aux_vrr(rterm, axis))
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
