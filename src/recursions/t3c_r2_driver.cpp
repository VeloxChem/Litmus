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

#include "t3c_r2_driver.hpp"

#include "axes.hpp"

T3CR2Driver::T3CR2Driver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T3CR2Driver::is_r2(const R2CTerm& rterm) const
{
    if (!(rterm.prefixes()).empty())
    {
        return false;
    }
    
    if (rterm.integrand().name() != "GR2(r)")
    {
        return false;
    }
    else
    {
        return true;
    }
}

R2CDist
T3CR2Driver::aux_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt(rterm);
    
    if (is_r2(rterm))
    {
        auto tval = rterm.replace(OperatorComponent("G(r)"));
        
        // R2(GC) contribution
        
        auto r1term = tval;
        
        r1term.add(Factor("r2gc", "rgc2"), Fraction(1));
        
        t2crt.add(r1term);
        
        // bra - 1 contributions
        
        for (const auto axis : "xyz")
        {
            if (const auto r2val = tval.shift(axis, -1, 0))
            {
                const auto coord = _rxyz[axes::to_index(axis)];
                
                const auto na = tval[0][axis];
                
                auto x2val = *r2val;
                
                x2val.add(Factor("GC", "gc", coord), Fraction(2));
                
                x2val.add(Factor("1/geta", "gfe"), Fraction(na));
                
                t2crt.add(x2val);
            }
        }
        
        // ket - 1 contributions
        
        for (const auto axis : "xyz")
        {
            if (const auto r3val = tval.shift(axis, -1, 1))
            {
                const auto coord = _rxyz[axes::to_index(axis)];
                
                const auto nb = tval[1][axis];
                
                auto x3val = *r3val;
                
                x3val.add(Factor("GC", "gc", coord), Fraction(2));
                
                x3val.add(Factor("1/geta", "gfe"), Fraction(nb));
                
                t2crt.add(x3val);
            }
        }
        
        // bra - 1, ket - 1 contributions
        
        for (const auto axis : "xyz")
        {
            if (const auto r4val = tval.shift(axis, -1, 0))
            {
                if (const auto r5val = r4val->shift(axis, -1, 1))
                {
                    const auto na = tval[0][axis];
                    
                    const auto nb = tval[1][axis];
                    
                    auto x5val = *r5val;
                    
                    x5val.add(Factor("1/geta2", "gfe2"), Fraction(2 * na * nb));
                    
                    t2crt.add(x5val);
                }
            }
        }
        
        // bra - 2 contributions
        
        for (const auto axis : "xyz")
        {
            if (const auto r6val = tval.shift(axis, -2, 0))
            {
                const auto na = tval[0][axis];
                
                auto x6val = *r6val;
                
                x6val.add(Factor("1/geta2", "gfe2"), Fraction(na * (na - 1)));
                
                t2crt.add(x6val);
            }
        }
        
        // ket - 2 contributions
        
        for (const auto axis : "xyz")
        {
            if (const auto r7val = tval.shift(axis, -2, 1))
            {
                const auto nb = tval[1][axis];
                
                auto x7val = *r7val;
                
                x7val.add(Factor("1/geta2", "gfe2"), Fraction(nb * (nb - 1)));
                
                t2crt.add(x7val);
            }
        }
    
        // operator delta contribution
        
        auto r6term = tval;
        
        r6term.add(Factor("1/geta", "gfe"), Fraction(3));
        
        t2crt.add(r6term);
    }
    
    return t2crt;
}
