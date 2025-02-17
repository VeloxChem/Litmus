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

#include "t3c_geom_100_eri_driver.hpp"

#include "t4c_geom_01_hrr_eri_driver.hpp"

#include "axes.hpp"
#include "t4c_utils.hpp"

T3CGeom100ElectronRepulsionDriver::T3CGeom100ElectronRepulsionDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T3CGeom100ElectronRepulsionDriver::is_electron_repulsion(const R3CTerm& rterm) const
{
    if (rterm.prefixes_order() != std::vector<int>({1, 0, 0}))
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
T3CGeom100ElectronRepulsionDriver::bra_aux_hrr(const R3CTerm& rterm,
                                               const char     axis) const
{
    if (const auto tval = rterm.shift_prefix(axis, -1, 0, false))
    {
        R3CDist t3crt(rterm);
        
        if (auto r1val = tval->shift(axis, 1, 0))
        {
            auto x1val = *r1val;
            
            x1val.clear_prefixes();
            
            t3crt.add(x1val);
        }
        
        if (auto r2val = tval->shift(axis, -1, 0))
        {
            auto x2val = *r2val;
            
            x2val.clear_prefixes();
            
            x2val.scale(Fraction(-(*tval)[0][axis]));
            
            t3crt.add(x2val);
        }
        
        return t3crt;
    }
    else
    {
        return std::nullopt;
    }
}

R3CDist
T3CGeom100ElectronRepulsionDriver::apply_bra_aux_hrr(const R3CTerm& rterm) const
{
    R3CDist t3crt;
    
    if (const auto prefixes = rterm.integral().prefixes(); !prefixes.empty())
    {
        const auto axis = prefixes[1].shape().primary();
        
        if (const auto trec = bra_aux_hrr(rterm, axis)) t3crt = *trec;
    }
   
    return t3crt;
}

std::optional<R3CDist>
T3CGeom100ElectronRepulsionDriver::bra_hrr(const R3CTerm& rterm,
                                           const char     axis) const
{
    if (!is_electron_repulsion(rterm)) return std::nullopt;
    
    if (const auto tval = rterm.shift_prefix(axis, -1, 0, false))
    {
        R3CDist t3crt(rterm);
        
        if (auto r1val = tval->shift(axis, 1, 0))
        {
            auto x1val = *r1val;
            
            x1val.clear_prefixes();
            
            t3crt.add(x1val);
        }
        
        if (auto r2val = tval->shift(axis, -1, 0))
        {
            auto x2val = *r2val;
            
            x2val.clear_prefixes();
            
            x2val.scale(Fraction(-(*tval)[0][axis]));
            
            t3crt.add(x2val);
        }
        
        return t3crt;
    }
    else
    {
        return std::nullopt;
    }
}

R3CDist
T3CGeom100ElectronRepulsionDriver::apply_bra_hrr(const R3CTerm& rterm) const
{
    R3CDist t3crt;
    
    size_t nints = 4;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = bra_hrr(rterm, axis))
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
