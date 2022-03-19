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
    _rxyz = {TensorComponent(1, 0, 0), TensorComponent(0, 1, 0),
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

R4CDist
EriDriver::apply_bra_hrr(const R4CTerm&               rterm,
                               std::set<T4CIntegral>& vints) const
{
    R4CDist t4crt;
    
//    for (const auto axis : "xyz")
//    {
//        if (const auto tval = bra_hrr(rterm, axis))
//        {
//            const auto nterms = tval->unique_terms(vints);
//        }
//    }
//    
//    std::array<R4CDist, 3> vdists({t4crt, t4crt, t4crt});
//
//    for (int i = 0; i < 3; i++)
//    {
//        if (const auto taux = rterm.shift(axes[i], -1, 0))
//        {
//            // first term
//
//            if (const auto r1val = taux->shift(axes[i], 1, 1))
//            {
//                vdists[i].add(*r1val);
//            }
//
//            // second term
//
//            auto r2val = taux.value();
//
//            r2val.add(Factor("AB", "rab", rxyz[i]), Fraction(-1));
//
//            vdists[i].add(r2val);
//        }
//    }
    
    // select recurion variant
    
    //t4crt = vdists[0];
    
//    auto mints = vdists[i].unique_terms(vin)
//
//    for (int i = 1; i < 3; i++)
//    {
//
//    }
    
   return t4crt;
}
