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

#include "t2c_trans_one_driver.hpp"

#include "t2c_center_driver.hpp"
#include "axes.hpp"
#include "t2c_utils.hpp"

T2CTransOneDriver::T2CTransOneDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T2CTransOneDriver::is_auxilary(const R2CTerm& rterm) const
{
    if (rterm.integrand().shape().order() != 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::optional<R2CDist>
T2CTransOneDriver::bra_ket_vrr(const R2CTerm& rterm) const
{
    if (is_auxilary(rterm)) return std::nullopt;
    
    const auto axis = rterm.integrand().shape().primary();
    
    R2CDist t2crt(rterm);

    auto bterm = rterm.replace(OperatorComponent("R"));
    
    bterm.clear_prefixes(); 
        
    if (const auto r1val = bterm.shift(axis, 1, 0))
    {
        auto x1val = *r1val;
        
        x1val.add(Factor("b_e", "tbe"), Fraction(2));
        
        t2crt.add(x1val);
    }
    
    if (const auto r2val = bterm.shift(axis, 1, 1))
    {
        auto x2val = *r2val;
        
        x2val.add(Factor("k_e", "tke"), Fraction(2));
        
        t2crt.add(x2val);
    }
    
    if (const auto r3val = bterm.shift(axis, -1, 0))
    {
        auto x3val = *r3val;
        
        x3val.scale(Fraction(-bterm[0][axis]));

        t2crt.add(x3val);
    }
    
    if (const auto r4val = bterm.shift(axis, -1, 1))
    {
        auto x4val = *r4val;
        
        x4val.scale(Fraction(-bterm[1][axis]));

        t2crt.add(x4val);
    }
    
    return t2crt;
}

R2Group
T2CTransOneDriver::create_recursion(const VT2CIntegrals& vints) const
{
    // create reference group
    
    R2Group r2group;
    
    for (const auto& tcomp : vints)
    {
        auto rdist = bra_ket_vrr(R2CTerm(tcomp));
                
        r2group.add(*rdist);
    }
    
    r2group.simplify();
    
    return r2group;
}
