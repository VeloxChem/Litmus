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

#include "geom_center_driver.hpp"

#include "axes.hpp"

T2CCenterDriver::T2CCenterDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T2CCenterDriver::is_auxilary(const R2CTerm& rterm,
                             const int      index) const
{
    if (const auto nprefixes = rterm.prefixes().size(); index > nprefixes)
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::optional<R2CDist>
T2CCenterDriver::bra_ket_vrr(const R2CTerm& rterm,
                             const char     axis,
                             const int      index) const
{
    if (is_auxilary(rterm, index)) return std::nullopt;
    
    if (const auto tval = rterm.shift_prefix(axis, -1, index, true))
    {
        R2CDist t2crt(rterm);
        
        
        
        return t2crt;
    }
    else
    {
        return std::nullopt;
    }
}
