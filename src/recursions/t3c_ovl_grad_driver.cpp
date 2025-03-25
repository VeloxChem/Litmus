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

#include "t3c_ovl_grad_driver.hpp"

#include "axes.hpp"
#include "t2c_utils.hpp"

T3COverlapGradientDriver::T3COverlapGradientDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T3COverlapGradientDriver::is_overlap_gradient(const R2CTerm& rterm) const
{
    if (!(rterm.prefixes()).empty())
    {
        return false;
    }
    
    if (rterm.integrand().name() != "GX(r)")
    {
        return false;
    }
    else
    {
        return true;
    }
}

std::optional<R2CDist>
T3COverlapGradientDriver::aux_vrr(const R2CTerm& rterm,
                                  const char     axis) const
{
    if (!is_overlap_gradient(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift_operator(axis, -1))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        auto r1val = tval->replace(OperatorComponent("G(r)"));
        
        auto x1val = r1val;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("GC", "gc", coord), Fraction(1));
        
        x1val.add(Factor("c_e", "tce"), Fraction(2));
        
        t2crt.add(x1val);
        
        // second recursion term
        
        if (const auto r2val = r1val.shift(axis, -1, 0))
        {
            auto x2val = *r2val;
            
            const auto na = x1val[0][axis];
            
            x2val.add(Factor("1/geta", "gfe"), Fraction(na));
            
            x2val.add(Factor("c_e", "tce"), Fraction(2));
            
            t2crt.add(x2val);
        }
        
        // third recursion term
        
        if (const auto r3val = r1val.shift(axis, -1, 1))
        {
            auto x3val = *r3val;
            
            const auto nb = x1val[1][axis];
            
            x3val.add(Factor("1/geta", "gfe"), Fraction(nb));
            
            x3val.add(Factor("c_e", "tce"), Fraction(2));
            
            t2crt.add(x3val);
        }
        
        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}

R2CDist
T3COverlapGradientDriver::apply_aux_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    size_t nints = 4;
    
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
