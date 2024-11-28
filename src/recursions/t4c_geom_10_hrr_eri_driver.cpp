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

#include "t4c_geom_10_hrr_eri_driver.hpp"

#include "axes.hpp"
#include "t4c_utils.hpp"

T4CGeom10HrrElectronRepulsionDriver::T4CGeom10HrrElectronRepulsionDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T4CGeom10HrrElectronRepulsionDriver::is_electron_repulsion(const R4CTerm& rterm) const
{
    if ((rterm.prefixes_order() != std::vector<int>({1, 0, 0, 0})) &&
        (rterm.prefixes_order() != std::vector<int>({0, 0, 1, 0})) &&
        (rterm.prefixes_order() != std::vector<int>({1, 0, 1, 0})))
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

std::optional<R4CDist>
T4CGeom10HrrElectronRepulsionDriver::bra_hrr(const R4CTerm& rterm,
                                             const char     axis) const
{
    if (!is_electron_repulsion(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 0))
    {
        R4CDist t4crt(rterm);
        
        // first recursion term
        
        auto x1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("BA", "ab", coord), Fraction(-1));
        
        t4crt.add(x1val);
        
        // second recursion term
        
        auto x2val = *tval;
        
        if (const auto prefixes = x2val.integral().prefixes(); !prefixes.empty())
        {
            if (prefixes[0].shape().primary() == axis)
            {
                x2val.clear_prefixes();
                
                x2val.scale(Fraction(-1));
                
                t4crt.add(x2val); 
            }
        }
        
        // third recursion term
        
        if (const auto r2val = tval->shift(axis, 1, 1))
        {
            t4crt.add(*r2val);
        }
        
        return t4crt;
    }
    else
    {
        return std::nullopt;
    }
}

std::optional<R4CDist>
T4CGeom10HrrElectronRepulsionDriver::ket_hrr(const R4CTerm& rterm,
                                             const char     axis) const
{
    if (!is_electron_repulsion(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift(axis, -1, 2))
    {
        R4CDist t4crt(rterm);
        
        // first recursion term
        
        auto x1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("DC", "cd", coord), Fraction(-1));
        
        t4crt.add(x1val);
        
        // second recursion term
        
        auto x2val = *tval;
        
        if (const auto prefixes = x2val.integral().prefixes(); !prefixes.empty())
        {
            if (prefixes[2].shape().primary() == axis)
            {
                x2val.clear_prefixes();
                
                x2val.scale(Fraction(-1));
                
                t4crt.add(x2val);
            }
        }
        
        // third recursion term
        
        if (const auto r2val = tval->shift(axis, 1, 3))
        {
            t4crt.add(*r2val);
        }
        
        return t4crt;
    }
    else
    {
        return std::nullopt;
    }
}

R4CDist
T4CGeom10HrrElectronRepulsionDriver::apply_bra_hrr(const R4CTerm& rterm) const
{
    R4CDist t4crt;
    
    size_t nints = 4;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = bra_hrr(rterm, axis))
        {
            if (const auto nterms = trec->terms(); nterms < nints)
            {
                t4crt = *trec;
                
                nints = nterms;
            }
        }
    }
    
    return t4crt;
}

R4CDist
T4CGeom10HrrElectronRepulsionDriver::apply_ket_hrr(const R4CTerm& rterm) const
{
    R4CDist t4crt;
    
    size_t nints = 4;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = ket_hrr(rterm, axis))
        {
            if (const auto nterms = trec->terms(); nterms < nints)
            {
                t4crt = *trec;
                
                nints = nterms;
            }
        }
    }
    
    return t4crt;
}

std::optional<R4CDist>
T4CGeom10HrrElectronRepulsionDriver::bra_aux_hrr(const R4CTerm& rterm,
                                                 const char     axis) const
{
    if (!is_electron_repulsion(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift_prefix(axis, -1, 0))
    {
        R4CDist t4crt(rterm);
        
        // first recursion term
        
        auto x1val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("BA", "ab", coord), Fraction(-1));
        
        t4crt.add(x1val);
        
        // third recursion term
        
        if (const auto r2val = tval->shift(axis, 1, 1))
        {
            auto x2val = *r2val;
            
            //x2val.clear_prefixes();
            
            t4crt.add(x2val);
        }
        
        return t4crt;
    }
    else
    {
        return std::nullopt;
    }
}

std::optional<R4CDist>
T4CGeom10HrrElectronRepulsionDriver::ket_aux_hrr(const R4CTerm& rterm,
                                                 const char     axis) const
{
    if (!is_electron_repulsion(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift_prefix(axis, -1, 2))
    {
        R4CDist t4crt(rterm);
        
        // first recursion term
        
        auto x1val = *tval;
        
        x1val.clear_prefixes();
        
        const auto coord = _rxyz[axes::to_index(axis)];
        
        x1val.add(Factor("DC", "cd", coord), Fraction(-1));
        
        t4crt.add(x1val);
        
        // third recursion term
        
        if (const auto r2val = tval->shift(axis, 1, 3))
        {
            auto x2val = *r2val;
            
            x2val.clear_prefixes();
            
            t4crt.add(x2val);
        }
        
        return t4crt;
    }
    else
    {
        return std::nullopt;
    }
}

R4CDist
T4CGeom10HrrElectronRepulsionDriver::apply_bra_aux_hrr(const R4CTerm& rterm) const
{
    R4CDist t4crt;
    
    if (const auto prefixes = rterm.integral().prefixes(); !prefixes.empty())
    {
        const auto axis = prefixes[0].shape().primary();
        
        if (const auto trec = bra_aux_hrr(rterm, axis)) t4crt = *trec;
    }
   
    return t4crt;
}

R4CDist
T4CGeom10HrrElectronRepulsionDriver::apply_ket_aux_hrr(const R4CTerm& rterm) const
{
    R4CDist t4crt;
    
    if (const auto prefixes = rterm.integral().prefixes(); !prefixes.empty())
    {
        const auto axis = prefixes[2].shape().primary();
        
        if (const auto trec = ket_aux_hrr(rterm, axis)) t4crt = *trec;
    }
   
    return t4crt;
}
