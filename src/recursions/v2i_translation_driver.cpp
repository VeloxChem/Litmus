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

#include "v2i_translation_driver.hpp"

bool
V2ITranslationDriver::is_auxiliary(const I2CIntegral& integral) const
{
    if (integral.integrand().shape() == Tensor(0))
    {
        return true;
    }
    else
    {
        return false;
    }
}

SI2CIntegrals
V2ITranslationDriver::operator_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;

    if (is_auxiliary(integral)) return tints;
    
    const auto gorder = integral.prefixes_sum_order();
    
    const auto gorders = integral.prefixes_order();
    
    // simple translation : first derivative of operator
    
    if ((gorder == 0) && (integral.integrand().shape() == Tensor(1)))
    {
        if (auto rint = integral.base().shift_operator(-1))
        {
            if (const auto r1val = rint->shift(1, 0))
            {
                tints.insert(*r1val);
            }
            
            if (const auto r2val = rint->shift(-1, 0))
            {
                tints.insert(*r2val);
            }
            
            if (const auto r1val = rint->shift(1, 1))
            {
                tints.insert(*r1val);
            }
            
            if (const auto r2val = rint->shift(-1, 1))
            {
                tints.insert(*r2val);
            }
        }
        
        return tints; 
    }
    
    // simple translation : second derivative of operator
    
    if ((gorder == 0) && (integral.integrand().shape() == Tensor(2)))
    {
        if (auto rint = integral.shift_operator(-2))
        {
            if (const auto r1val = rint->shift_prefix(2, 0, false))
            {
                tints.insert(*r1val);
            }
            
            if (const auto r2val = rint->shift_prefix(1, 1, false))
            {
                if (const auto t2val = r2val->shift_prefix(1, 0, false))
                {
                    tints.insert(*t2val); 
                }
            }
            
            if (const auto r3val = rint->shift_prefix(2, 1, false))
            {
                tints.insert(*r3val);
            }
        }
        
        return tints;
    }
    
    // mixed translation : first derivatives of operator and bra side
    if ((gorders == std::vector<int>({1, 0})) && (integral.integrand().shape() == Tensor(1)))
    {
        if (auto rint = integral.shift_operator(-1))
        {
            if (const auto r1val = rint->shift_prefix(1, 0, false))
            {
                tints.insert(*r1val);
            }
            
            if (const auto r2val = rint->shift_prefix(1, 1, false))
            {
                tints.insert(*r2val);
            }
        }
        
        return tints;
    }
    
    return tints;
}

