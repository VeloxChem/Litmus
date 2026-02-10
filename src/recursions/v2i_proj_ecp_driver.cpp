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

#include "v2i_proj_ecp_driver.hpp"

#include <cmath>

bool
V2IProjectedECPDriver::is_projected_ecp(const M2Integral& integral) const
{
    if (!(integral.second.prefixes()).empty())
    {
        return false;
    }
    
    if (integral.second.integrand() != Operator("U_l"))
    {
        return false;
    }
    else
    {
        return true;
    }
}

SM2Integrals
V2IProjectedECPDriver::bra_vrr(const M2Integral& integral) const
{
    SM2Integrals tints;
    
    if (!is_projected_ecp(integral)) return tints;
    
    auto [order, rint] = integral;
    
    if (const auto tval = rint.shift(-1, 0))
    {
        // set up recursion orders
        
        auto morder = order; morder[0] += 1;
        
        auto pq_order = order; pq_order[1] -= 1; pq_order[2] += 1;
        
        // first recursion term

        tints.insert({order, *tval});
    
        // second recursion term
        
        tints.insert({morder, *tval});
        
        // third recursion term
        
        if (pq_order[1] > 0)
        {
            tints.insert({pq_order, *tval});
        }
        
        // fourth, fifth and sixth recursion terms
        
        if (const auto r2val = tval->shift(-1, 0))
        {
            tints.insert({order, *r2val});
            
            tints.insert({morder, *r2val});
            
            if (pq_order[1] > 0)
            {
                tints.insert({pq_order, *r2val});
            }
        }
        
        // set up lower order projectors
    
        const auto l = rint.order();
        
        const int l1p = (int)(std::floor(0.5 * (l - 1)));
        
        const int l2p = (int)(std::floor(0.5 * (l - 2)));
        
        // (l - 1) / 2 terms
        
        for (int k = 0; k <= l1p; k++)
        {
            // set up adjusted order
            
            auto mpq_order = order;
            
            mpq_order[0] += k;
            
            mpq_order[1] += k;
            
            mpq_order[2] += 1;
            
            // first and second terms
            
            if (const auto r3val = tval->shift_order(-2 * k - 1))
            {
                tints.insert({mpq_order, *r3val});
                
                if (const auto r4val = r3val->shift(-1, 1))
                {
                    tints.insert({mpq_order, *r4val});
                }
            }
        }
        
        // (l - 2) / 2 terms
        
        for (int k = 0; k <= l2p; k++)
        {
            // set up adjusted order
            
            auto mpq_order = order;
            
            mpq_order[0] += k + 1;
            
            mpq_order[1] += k;
            
            mpq_order[2] += 1;
            
            // first and second terms
            
            if (const auto r3val = tval->shift_order(-2 * k - 2))
            {
                tints.insert({mpq_order, *r3val});
                
                if (const auto r4val = r3val->shift(-1, 0))
                {
                    tints.insert({mpq_order, *r4val});
                }
            }
        }
    }
        
    return tints;
}

SM2Integrals
V2IProjectedECPDriver::ket_vrr(const M2Integral& integral) const
{
    SM2Integrals tints;
    
    if (!is_projected_ecp(integral)) return tints;
    
    auto [order, rint] = integral;
    
    if (const auto tval = rint.shift(-1, 1))
    {
        // set up recursion orders
        
        auto morder = order; morder[0] += 1;
        
        auto pq_order = order; pq_order[1] += 1; pq_order[2] -= 1;
        
        // first recursion term

        tints.insert({order, *tval});
    
        // second recursion term
        
        tints.insert({morder, *tval});
        
        // third recursion term
        
        if (pq_order[2] > 0)
        {
            tints.insert({pq_order, *tval});
        }
        
        // fourth, fifth and sixth recursion terms
        
        if (const auto r2val = tval->shift(-1, 1))
        {
            tints.insert({order, *r2val});
            
            tints.insert({morder, *r2val});
            
            if (pq_order[2] > 0)
            {
                tints.insert({pq_order, *r2val});
            }
        }
        
        // set up lower order projectors
    
        const auto l = rint.order();
        
        const int l1p = (int)(std::floor(0.5 * (l - 1)));
        
        const int l2p = (int)(std::floor(0.5 * (l - 2)));
        
        // (l - 1) / 2 terms
        
        for (int k = 0; k <= l1p; k++)
        {
            // set up adjusted order
            
            auto mpq_order = order;
            
            mpq_order[0] += k;
            
            mpq_order[1] += 1;
            
            mpq_order[2] += k;
            
            // first and second terms
            
            if (const auto r3val = tval->shift_order(-2 * k - 1))
            {
                tints.insert({mpq_order, *r3val});
                
                if (const auto r4val = r3val->shift(-1, 0))
                {
                    tints.insert({mpq_order, *r4val});
                }
            }
        }
        
        // (l - 2) / 2 terms
        
        for (int k = 0; k <= l2p; k++)
        {
            // set up adjusted order
            
            auto mpq_order = order;
            
            mpq_order[0] += k + 1;
            
            mpq_order[1] += 1;
            
            mpq_order[2] += k;
            
            // first and second terms
            
            if (const auto r3val = tval->shift_order(-2 * k - 2))
            {
                tints.insert({mpq_order, *r3val});
                
                if (const auto r4val = r3val->shift(-1, 1))
                {
                    tints.insert({mpq_order, *r4val});
                }
            }
        }
    }
        
    return tints;
}

SM2Integrals
V2IProjectedECPDriver::apply_bra_vrr(const M2Integral& integral) const
{
    SM2Integrals tints;
    
    if (integral.second[0] > 0)
    {
        SM2Integrals rtints({integral, });
                
        while (!rtints.empty())
        {
            SM2Integrals new_rtints;
                
            for (const auto& rtint : rtints)
            {
                if (rtint.second[0] != 0)
                {
                   const auto ctints = bra_vrr(rtint);
                    
                   for (const auto& ctint : ctints)
                   {
                       tints.insert(ctint);
                       
                       if (ctint.second[0] != 0)
                       {
                           new_rtints.insert(ctint);
                       }
                   }
                }
                else
                {
                    tints.insert(rtint);
                }
            }
            
            rtints = new_rtints;
        }
    }
   
    tints.insert(integral);
    
    return tints;
}

SM2Integrals
V2IProjectedECPDriver::apply_ket_vrr(const M2Integral& integral) const
{
    SM2Integrals tints;
    
    if (integral.second[1] > 0)
    {
        SM2Integrals rtints({integral, });
                
        while (!rtints.empty())
        {
            SM2Integrals new_rtints;
                
            for (const auto& rtint : rtints)
            {
                if (rtint.second[1] != 0)
                {
                   const auto ctints = ket_vrr(rtint);
                    
                   for (const auto& ctint : ctints)
                   {
                       tints.insert(ctint);
                       
                       if (ctint.second[1] != 0)
                       {
                           new_rtints.insert(ctint);
                       }
                   }
                }
                else
                {
                    tints.insert(rtint);
                }
            }
            
            rtints = new_rtints;
        }
    }
   
    tints.insert(integral);
    
    return tints;
}

SM2Integrals
V2IProjectedECPDriver::apply_recursion(const SM2Integrals& integrals) const
{
    SM2Integrals tints;
    
    for (const auto& integral : integrals)
    {
        tints.insert(integral);
        
        for (const auto& bintegral : apply_bra_vrr(integral))
        {
            if (bintegral.second[0] == 0)
            {
                const auto ctints = apply_ket_vrr(bintegral);

                tints.insert(ctints.cbegin(), ctints.cend());
            }
            else
            {
                tints.insert(bintegral);
            }
        }
    }
    
    return tints;
}

SM2Integrals
V2IProjectedECPDriver::create_recursion(const SM2Integrals& integrals) const
{
    SM2Integrals tints;
    
    for (const auto& integral : integrals)
    {
        if (is_projected_ecp(integral))
        {
            const auto ctints = apply_recursion({integral, });
            
            tints.insert(ctints.cbegin(), ctints.cend());
        }
        else
        {
            tints.insert(integral);
        }
    }
    
    return tints;
}
