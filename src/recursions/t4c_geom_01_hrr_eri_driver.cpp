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

#include "t4c_geom_01_hrr_eri_driver.hpp"

#include "axes.hpp"
#include "t4c_utils.hpp"

T4CGeom01HrrElectronRepulsionDriver::T4CGeom01HrrElectronRepulsionDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T4CGeom01HrrElectronRepulsionDriver::is_electron_repulsion(const R4CTerm& rterm) const
{
    if (rterm.prefixes_order() != std::vector<int>({0, 1, 0, 0}))
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
T4CGeom01HrrElectronRepulsionDriver::bra_aux_hrr(const R4CTerm& rterm,
                                                 const char     axis) const
{
    if (const auto tval = rterm.shift_prefix(axis, -1, 1, false))
    {
        R4CDist t4crt(rterm);
        
        if (auto r1val = tval->shift(axis, 1, 1))
        {
            auto x1val = *r1val;
            
            x1val.clear_prefixes();
            
            t4crt.add(x1val);
        }
        
        if (auto r2val = tval->shift(axis, -1, 1))
        {
            auto x2val = *r2val;
            
            x2val.clear_prefixes(); 
            
            x2val.scale(Fraction(-(*tval)[1][axis]));
            
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
T4CGeom01HrrElectronRepulsionDriver::apply_bra_aux_hrr(const R4CTerm& rterm) const
{
    R4CDist t4crt;
    
    if (const auto prefixes = rterm.integral().prefixes(); !prefixes.empty())
    {
        const auto axis = prefixes[1].shape().primary();
        
        if (const auto trec = bra_aux_hrr(rterm, axis)) t4crt = *trec;
    }
   
    return t4crt;
}

std::optional<R4CDist>
T4CGeom01HrrElectronRepulsionDriver::bra_hrr(const R4CTerm& rterm,
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
        
        if (auto r2val = tval->shift_prefix(axis, -1, 1, false))
        {
            auto x2val = *r2val;
            
            x2val.clear_prefixes();
            
            t4crt.add(x2val);
        }
        
        // third recursion term
        
        if (const auto r4val = tval->shift(axis, 1, 1))
        {
            t4crt.add(*r4val);
        }
        
        return t4crt;
    }
    else
    {
        return std::nullopt;
    }
}

R4CDist
T4CGeom01HrrElectronRepulsionDriver::apply_bra_hrr(const R4CTerm& rterm) const
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
