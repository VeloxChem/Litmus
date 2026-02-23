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

#include "t2c_proj_ecp_driver.hpp"

#include "axes.hpp"
#include "t2c_utils.hpp"

T2CProjectedECPDriver::T2CProjectedECPDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T2CProjectedECPDriver::is_projected_ecp(const R2CTerm& rterm) const
{
    if (!(rterm.prefixes()).empty())
    {
        return false;
    }
    
    if (rterm.integrand() != OperatorComponent("U_l"))
    {
        return false;
    }
    else
    {
        return true;
    }
}

std::optional<R2CDist>
T2CProjectedECPDriver::bra_vrr(const R2CTerm& rterm,
                               const char     axis) const
{
    if (!is_projected_ecp(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 0))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        auto x1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("RA", "ra", coord), Fraction(1));
        
        x1val.add(Factor("a-z/z", "faz"), Fraction(1));
        
        t2crt.add(x1val);
        
        // second recursion term
        
        auto x2val = *tval;
        
        x2val.add(Factor("RA", "ra", coord), Fraction(1));
        
        x2val.add(Factor("a", "fa"), Fraction(2));
        
        x2val.add(Factor("b/z", "fbzi"), Fraction(1));
        
        x2val.add(Factor("b/z", "fbzi"), Fraction(1));
        
        x2val.add(Factor("m", "m"), Fraction(1));
        
        t2crt.add(x2val);
        
        // third and fourth recursion term
        
        if (const auto rval = tval->shift(axis, -1, 0))
        {
            auto x3val = *rval;
            
            const auto na = x1val[0][axis];
            
            x3val.add(Factor("1/2z", "fzi"), Fraction(na));
            
            t2crt.add(x3val);
            
            auto x4val = *rval;
            
            x4val.add(Factor("b/z", "fbzi"), Fraction(na));
            
            x4val.add(Factor("b/z", "fbzi"), Fraction(1));
            
            x4val.add(Factor("m", "m"), Fraction(1));
            
            t2crt.add(x4val);
        }
        
        // set up angular momentum values
        
        const auto l = rterm.order();
        
        const int l1p = (int)(std::floor(0.5 * (l - 1)));
        
        const int l2p = (int)(std::floor(0.5 * (l - 2)));
        
        // (l - 1) / 2 terms
        
        for (int k = 0; k <= l1p; k++)
        {
            if (const auto rkval = tval->shift_order(- 2 * k - 1))
            {
                auto x5val = *rkval;
                
                x5val.add(Factor("RB", "rb", coord), Fraction(1));
                
                x5val.add(Factor("b/z", "fbzi"), Fraction(2 * l + 1));
                
                for (int i = 0; i < 2 * k; i++)
                {
                    x5val.add(Factor("2ab/z", "f2abz"), Fraction(1));
                }
                
                for (int m = 1; m <= k; m++)
                {
                    x5val.add(Factor("m", "m"), Fraction(1));
                }
                
                for (int p = 1; p <= k; p++)
                {
                    x5val.add(Factor("p", "p"), Fraction(1));
                }
                
                x5val.add(Factor("q", "q"), Fraction(1));
                
                t2crt.add(x5val);
                
                if (const auto rmval = rkval->shift(axis, -1, 1))
                {
                    auto x6val = *rkval;
                    
                    const auto nb = x1val[1][axis];
                    
                    x6val.add(Factor("1/b", "fbi"), Fraction(nb, 2));
                    
                    x6val.add(Factor("b/z", "fbzi"), Fraction(2 * l + 1));
                    
                    for (int i = 0; i < 2 * k; i++)
                    {
                        x6val.add(Factor("2ab/z", "f2abz"), Fraction(1));
                    }
                    
                    for (int m = 1; m <= k; m++)
                    {
                        x6val.add(Factor("m", "m"), Fraction(1));
                    }
                    
                    for (int p = 1; p <= k; p++)
                    {
                        x6val.add(Factor("p", "p"), Fraction(1));
                    }
                    
                    x6val.add(Factor("q", "q"), Fraction(1));
                    
                    t2crt.add(x6val);
                }
            }
        }
        
        // (l - 2) / 2 terms
        
        for (int k = 0; k <= l2p; k++)
        {
            if (const auto rkval = tval->shift_order(- 2 * k - 2))
            {
                auto x7val = *rkval;
                
                x7val.add(Factor("RA", "ra", coord), Fraction(-1));
                
                x7val.add(Factor("b/z", "fbzi"), Fraction(2 * l + 1));
                
                for (int i = 0; i < 2 * k + 1; i++)
                {
                    x7val.add(Factor("2ab/z", "f2abz"), Fraction(1));
                }
                
                for (int m = 1; m <= (k + 1); m++)
                {
                    x7val.add(Factor("m", "m"), Fraction(1));
                }
                
                for (int p = 1; p <= k; p++)
                {
                    x7val.add(Factor("p", "p"), Fraction(1));
                }
                
                x7val.add(Factor("q", "q"), Fraction(1));
                
                t2crt.add(x7val);
                
                if (const auto rmval = rkval->shift(axis, -1, 0))
                {
                    auto x8val = *rkval;
                    
                    const auto na = x1val[0][axis];
                    
                    x8val.add(Factor("1/a", "fai"), Fraction(-na, 2));
                    
                    x8val.add(Factor("b/z", "fbzi"), Fraction(2 * l + 1));
                    
                    for (int i = 0; i < 2 * k + 1; i++)
                    {
                        x8val.add(Factor("2ab/z", "f2abz"), Fraction(1));
                    }
                    
                    for (int m = 1; m <= (k + 1); m++)
                    {
                        x8val.add(Factor("m", "m"), Fraction(1));
                    }
                    
                    for (int p = 1; p <= k; p++)
                    {
                        x8val.add(Factor("p", "p"), Fraction(1));
                    }
                    
                    x8val.add(Factor("q", "q"), Fraction(1));
                    
                    t2crt.add(x8val);
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
T2CProjectedECPDriver::ket_vrr(const R2CTerm& rterm,
                           const char     axis) const
{
    if (!is_projected_ecp(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 1))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        auto x1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("RB", "rb", coord), Fraction(1));
        
        x1val.add(Factor("b-z/z", "fbz"), Fraction(1));
        
        t2crt.add(x1val);
        
        // second recursion term
        
        auto x2val = *tval;
        
        x2val.add(Factor("RB", "rb", coord), Fraction(1));
        
        x2val.add(Factor("b", "fb"), Fraction(2));
        
        x2val.add(Factor("a/z", "fazi"), Fraction(1));
        
        x2val.add(Factor("a/z", "fazi"), Fraction(1));
        
        x2val.add(Factor("p", "p"), Fraction(1));
        
        t2crt.add(x2val);
        
        // third and fourth recursion term
        
        if (const auto rval = tval->shift(axis, -1, 1))
        {
            auto x3val = *rval;
            
            const auto nb = x1val[1][axis];
            
            x3val.add(Factor("1/2z", "fzi"), Fraction(nb));
            
            t2crt.add(x3val);
            
            auto x4val = *rval;
            
            x4val.add(Factor("a/z", "fazi"), Fraction(nb));
            
            x4val.add(Factor("a/z", "fazi"), Fraction(1));
            
            x4val.add(Factor("p", "p"), Fraction(1));
            
            t2crt.add(x4val);
        }
        
        // set up angular momentum values
        
        const auto l = rterm.order();
        
        const int l1p = (int)(std::floor(0.5 * (l - 1)));
        
        const int l2p = (int)(std::floor(0.5 * (l - 2)));
        
        // (l - 1) / 2 terms
        
        for (int k = 0; k <= l1p; k++)
        {
            if (const auto rkval = tval->shift_order(- 2 * k - 1))
            {
                auto x5val = *rkval;
                
                x5val.add(Factor("RA", "ra", coord), Fraction(1));
                
                x5val.add(Factor("a/z", "fazi"), Fraction(2 * l + 1));
                
                for (int i = 0; i < 2 * k; i++)
                {
                    x5val.add(Factor("2ab/z", "f2abz"), Fraction(1));
                }
                
                for (int m = 1; m <= k; m++)
                {
                    x5val.add(Factor("m", "m"), Fraction(1));
                }
                
                for (int p = 1; p <= k; p++)
                {
                    x5val.add(Factor("p", "p"), Fraction(1));
                }
                
                x5val.add(Factor("q", "q"), Fraction(1));
                
                t2crt.add(x5val);
            }
        }
        
        // (l - 2) / 2 terms
        
        for (int k = 0; k <= l2p; k++)
        {
            if (const auto rkval = tval->shift_order(- 2 * k - 2))
            {
                auto x6val = *rkval;
                
                x6val.add(Factor("RB", "rb", coord), Fraction(-1));
                
                x6val.add(Factor("a/z", "fazi"), Fraction(2 * l + 1));
                
                for (int i = 0; i < 2 * k + 1; i++)
                {
                    x6val.add(Factor("2ab/z", "f2abz"), Fraction(1));
                }
                
                for (int m = 1; m <= k; m++)
                {
                    x6val.add(Factor("m", "m"), Fraction(1));
                }
                
                for (int p = 1; p <= (k + 1); p++)
                {
                    x6val.add(Factor("p", "p"), Fraction(1));
                }
                
                x6val.add(Factor("q", "q"), Fraction(1));
                
                t2crt.add(x6val);
                
                if (const auto rmval = rkval->shift(axis, -1, 0))
                {
                    auto x8val = *rkval;
                    
                    const auto nb = x1val[1][axis];
                    
                    x8val.add(Factor("1/b", "fbi"), Fraction(-nb, 2));
                    
                    x8val.add(Factor("a/z", "fazi"), Fraction(2 * l + 1));
                    
                    for (int i = 0; i < 2 * k + 1; i++)
                    {
                        x8val.add(Factor("2ab/z", "f2abz"), Fraction(1));
                    }
                    
                    for (int m = 1; m <= k; m++)
                    {
                        x8val.add(Factor("m", "m"), Fraction(1));
                    }
                    
                    for (int p = 1; p <= (k + 1); p++)
                    {
                        x8val.add(Factor("p", "p"), Fraction(1));
                    }
                    
                    x8val.add(Factor("q", "q"), Fraction(1));
                    
                    t2crt.add(x8val);
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
T2CProjectedECPDriver::red_bra_vrr(const R2CTerm& rterm,
                                   const char     axis) const
{
    if (!is_projected_ecp(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 0))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        auto x1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("RA", "ra", coord), Fraction(1));
        
        x1val.add(Factor("1/a", "fai"), Fraction(1));
        
        x1val.add(Factor("fp", "fp"), Fraction(1));
        
        x1val.add(Factor("q", "q"), Fraction(1));
        
        t2crt.add(x1val);
        
        // second term
        
        if (const auto rval = tval->shift(axis, -1, 0))
        {
            auto x2val = *rval;
            
            const auto na = x1val[0][axis];
            
            x2val.add(Factor("1/a", "fai"), Fraction(na, 2));
            
            x2val.add(Factor("1/a", "fai"), Fraction(1));
            
            x2val.add(Factor("fp", "fp"), Fraction(1));
            
            x2val.add(Factor("q", "q"), Fraction(1));
            
            t2crt.add(x2val);
        }
        
        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}

std::optional<R2CDist>
T2CProjectedECPDriver::red_ket_vrr(const R2CTerm& rterm,
                                   const char     axis) const
{
    if (!is_projected_ecp(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 1))
    {
        R2CDist t2crt(rterm);
        
        // first recursion term
        
        auto x1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("RB", "rb", coord), Fraction(1));
        
        x1val.add(Factor("1/b", "fbi"), Fraction(1));
        
        x1val.add(Factor("fm", "fm"), Fraction(1));
        
        x1val.add(Factor("q", "q"), Fraction(1));
        
        t2crt.add(x1val);
        
        // second term
        
        if (const auto rval = tval->shift(axis, -1, 1))
        {
            auto x2val = *rval;
            
            const auto na = x1val[1][axis];
            
            x2val.add(Factor("1/b", "fbi"), Fraction(na, 2));
            
            x2val.add(Factor("1/b", "fbi"), Fraction(1));
            
            x2val.add(Factor("fm", "fm"), Fraction(1));
            
            x2val.add(Factor("q", "q"), Fraction(1));
            
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
T2CProjectedECPDriver::apply_bra_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    const auto l = rterm.order();
    
    size_t nints = 5;
    
    if (l > 0) nints += 4 * l - 2;
    
    if (l > 1) nints += 4 * l - 4;
    
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
T2CProjectedECPDriver::apply_ket_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    const auto l = rterm.order();
    
    size_t nints = 5;
    
    if (l > 0) nints += 2 * l - 1;
    
    if (l > 1) nints += 4 * l - 4;
    
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

R2CDist
T2CProjectedECPDriver::apply_red_bra_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    size_t nints = 3;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = red_bra_vrr(rterm, axis))
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
T2CProjectedECPDriver::apply_red_ket_vrr(const R2CTerm& rterm) const
{
    R2CDist t2crt;
    
    size_t nints = 3;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = red_ket_vrr(rterm, axis))
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
