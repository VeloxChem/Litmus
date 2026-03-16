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

#include "t2c_trans_gen_driver.hpp"

#include "t2c_center_driver.hpp"
#include "axes.hpp"
#include "t2c_utils.hpp"

T2CTransGenDriver::T2CTransGenDriver()
{
    _rxyz = {TensorComponent(1, 0, 0),
             TensorComponent(0, 1, 0),
             TensorComponent(0, 0, 1)};
}

bool
T2CTransGenDriver::is_auxilary(const R2CTerm& rterm) const
{
    if (rterm.integrand().shape().order() > 1)
    {
        return false;
    }
    
    if (rterm.integrand().shape().order() == 1)
    {
        const auto gorders = rterm.prefixes_order();
        
        if (gorders == std::vector<int>({1, 0})) return false;
        
        return true;
    }
    
    return true;
}

std::optional<R2CDist>
T2CTransGenDriver::operator_vrr(const R2CTerm& rterm) const
{
    if (is_auxilary(rterm)) return std::nullopt;
   
    R2CDist t2crt(rterm);
    
    const auto axis = rterm.integrand().shape().primary();

    const auto gorders = rterm.prefixes_order();

    // mixed g110 derivatives
    
    if ((gorders == std::vector<int>({1, 0})) && (rterm.integrand().shape().order() == 1))
    {
        auto bterm = rterm.replace(OperatorComponent("R"));
        
        if (const auto r1val = bterm.shift_prefix(axis, 1, 0))
        {
            t2crt.add(*r1val);
        }
        
        if (const auto r2val = bterm.shift_prefix(axis, 1, 1))
        {
            t2crt.add(*r2val);
        }
    }
    
    // g020 derivatives
    
    if ((gorders == std::vector<int>({0, 0})) && (rterm.integrand().shape().order() == 2))
    {
        if (const auto oterm = rterm.shift_operator(axis, -1))
        {
            const auto taxis = oterm->integrand().shape().primary();
            
            auto bterm = rterm.replace(OperatorComponent("R"));
            
            if (const auto r1val = bterm.shift_prefix(axis, 1, 0))
            {
                if (const auto r1aval = r1val->shift_prefix(taxis, 1, 0))
                {
                    t2crt.add(*r1aval);
                }
                
                if (const auto r2aval = r1val->shift_prefix(taxis, 1, 1))
                {
                    t2crt.add(*r2aval);
                }
            }
            
            if (const auto r1val = bterm.shift_prefix(axis, 1, 1))
            {
                if (const auto r1aval = r1val->shift_prefix(taxis, 1, 0))
                {
                    t2crt.add(*r1aval);
                }
                
                if (const auto r2aval = r1val->shift_prefix(taxis, 1, 1))
                {
                    t2crt.add(*r2aval);
                }
            }
        }
    }
    
    return t2crt;
}
