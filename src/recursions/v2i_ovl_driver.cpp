#include "v2i_ovl_driver.hpp"

bool
V2IOverlapDriver::is_overlap(const I2CIntegral& integral) const
{
    if (!(integral.prefixes()).empty())
    {
        return false;
    }
    
    if (integral.integrand() != Operator("1"))
    {
        return false;
    }
    else
    {
        return true;
    }
}

SI2CIntegrals
V2IOverlapDriver::apply_bra_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;
    
    if (integral[0] > 0)
    {
        SI2CIntegrals rtints({integral, });
                
        while (!rtints.empty())
        {
            SI2CIntegrals new_rtints;
                
            for (const auto& rtint : rtints)
            {
                
                
//                const auto cdist = _apply_bra_vrr(rec_terms[i]);
//
//                if (const auto nterms = cdist.terms(); nterms > 0)
//                {
//                    for (size_t j = 0; j < nterms; j++)
//                    {
//                        if (const auto rterm = cdist[j]; rterm.auxilary(0))
//                        {
//                            new_dist.add(rterm);
//                        }
//                        else
//                        {
//                            new_terms.push_back(rterm);
//                        }
//                    }
//                }
            }
                
            rtints = new_rtints;
        }
    }
    else
    {
        tints.insert(integral);
    }
    
    return tints;
}

SI2CIntegrals
V2IOverlapDriver::apply_ket_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;
    
    return tints;
}

SI2CIntegrals
V2IOverlapDriver::apply_recursion(const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints;
    
    for (const auto& integral : integrals)
    {
        for (const auto& bintegral : apply_bra_vrr({integral, }))
        {
            const auto ctints = apply_ket_vrr(bintegral);
            
            tints.insert(ctints.cbegin(), ctints.cend());
        }
    }
    
    return tints;
}

SI2CIntegrals
V2IOverlapDriver::create_recursion(const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints;
    
    for (const auto& integral : integrals)
    {
        if (is_overlap(integral))
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
