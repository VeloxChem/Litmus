#include "v2i_ovl_driver.hpp"

#include <iostream>

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
                if (rtint[0] != 0)
                {
                   const auto ctints = _bra_vrr(rtint);
                    
                   for (const auto& ctint : ctints)
                   {
                       tints.insert(ctint);
                       
                       if (ctint[0] != 0)
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

    std::cout << "Bra Expansion: ";
    
    for (const auto& tint : tints) std::cout << tint.label() << " : "; 
    
    std::cout << std::endl;
    
    return tints;
}

SI2CIntegrals
V2IOverlapDriver::apply_ket_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;
    
    if (integral[1] > 0)
    {
        SI2CIntegrals rtints({integral, });
                
        while (!rtints.empty())
        {
            SI2CIntegrals new_rtints;
                
            for (const auto& rtint : rtints)
            {
                const auto ctints = _ket_vrr(rtint);
                
                if (!ctints.empty())
                {
                    for (const auto& ctint : ctints)
                    {
                        if (ctint[1] == 0)
                        {
                            tints.insert(ctint);
                        }
                        else
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
    else
    {
        tints.insert(integral);
    }
    
    return tints;
}

SI2CIntegrals
V2IOverlapDriver::apply_recursion(const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints;
    
    for (const auto& integral : integrals)
    {
        tints.insert(integral);
        
        std::cout << "@integral " << integral.label() << std::endl;
        
        const auto btints = apply_bra_vrr({integral, });
        
//        for (const auto& bintegral : apply_bra_vrr({integral, }))
//        {
//            //std::cout << " bra --> " << bintegral.label() << std::endl;
//
//            const auto ctints = apply_ket_vrr(bintegral);
//
//            tints.insert(ctints.cbegin(), ctints.cend());
//        }
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

SI2CIntegrals
V2IOverlapDriver::_bra_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;
    
    if (!is_overlap(integral)) return tints;
    
    if (const auto tval = integral.shift(-1, 0))
    {
        auto x1val = *tval;

        tints.insert(x1val);
        
        // second recursion term
        
        if (const auto r2val = tval->shift(-1, 0))
        {
            tints.insert(*r2val);
        }
        
        // third recursion term
        
        if (const auto r3val = tval->shift(-1, 1))
        {
            tints.insert(*r3val);
        }
    }
    
    std::cout << "Bra VRR: ";
    
    for (const auto& tint : tints) std::cout << tint.label() << " : ";
    
    std::cout << std::endl;
    
    return tints;
}

SI2CIntegrals
V2IOverlapDriver::_ket_vrr(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;
    
    if (!is_overlap(integral)) return tints;
    
    if (const auto tval = integral.shift(-1, 1))
    {
        auto x1val = *tval;

        tints.insert(x1val);
        
        // second recursion term
        
        if (const auto r2val = tval->shift(-1, 1))
        {
            tints.insert(*r2val);
        }
    }
    
    return tints;
}
