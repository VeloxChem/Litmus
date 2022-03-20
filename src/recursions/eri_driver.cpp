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

#include "eri_driver.hpp"

#include <array>

#include "axes.hpp"

EriDriver::EriDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

std::optional<R4CDist>
EriDriver::bra_hrr(const R4CTerm& rterm,
                   const char     axis) const
{
    if (const auto tval = rterm.shift(axis, -1, 0))
    {
        R4CDist t4crt(rterm);
        
        // first recursion term
        
        if (const auto r1val = tval->shift(axis, 1, 1))
        {
            t4crt.add(*r1val);
        }
        
        // second recursion term
        
        auto r2val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
    
        r2val.add(Factor("AB", "rab", coord), Fraction(-1));
        
        t4crt.add(r2val);
        
        // done with recursion terms
        
        return t4crt;
    }
    else
    {
        return std::nullopt;
    }
}

std::optional<R4CDist>
EriDriver::ket_hrr(const R4CTerm& rterm,
                   const char     axis) const
{
    if (const auto tval = rterm.shift(axis, -1, 2))
    {
        R4CDist t4crt(rterm);
        
        // first recursion term
        
        if (const auto r1val = tval->shift(axis, 1, 3))
        {
            t4crt.add(*r1val);
        }
        
        // second recursion term
        
        auto r2val = *tval;
        
        const auto coord = _rxyz[axes::to_index(axis)];
    
        r2val.add(Factor("CD", "rcd", coord), Fraction(-1));
        
        t4crt.add(r2val);
        
        // done with recursion terms
        
        return t4crt;
    }
    else
    {
        return std::nullopt;
    }
}

R4CDist
EriDriver::apply_bra_hrr(const R4CTerm&       rterm,
                               ST4CIntegrals& sints) const
{
    R4CDist t4crt;
    
    int nints = 3;
    
    for (const auto axis : "xyz")
    {
        if (const auto trec = bra_hrr(rterm, axis))
        {
            const auto nterms = trec->count_new_integrals(sints);
            
            if (nterms < nints)
            {
                t4crt = *trec;
                
                nints = nterms;
            }
        }
    }
    
    const auto vints = t4crt.unique_integrals();
    
    sints.insert(vints.cbegin(), vints.cend());

    return t4crt;
}
