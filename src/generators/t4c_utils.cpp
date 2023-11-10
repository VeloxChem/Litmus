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

#include "t4c_utils.hpp"

#include <iostream>

#include "string_formater.hpp"

namespace t4c { // t2c namespace

std::string
integral_label(const I4CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    const auto prefixes = integral.prefixes();
    
    auto border = std::string("0");
    
    auto korder = std::string("0");
    
    const auto iorder = std::to_string(integrand.shape().order());
    
    if (const auto nterms = prefixes.size(); nterms > 0)
    {
        if (nterms >= 1) border =std::to_string(prefixes[0].shape().order());
        
        if (nterms >= 2) korder = std::to_string(prefixes[1].shape().order());
    }
    
    std::string suffix = "Geom" + border + iorder  + korder;
    
    if (integrand.name() == "1/|r-r'|")
    {
        return (prefixes.empty()) ? "ElectronRepulsion" : "ElectronRepulsion" + suffix;
    }
    
    return std::string();
}

std::pair<size_t, std::string>
diag_compute_func_name(const I4CIntegral& integral)
{
    const auto label = "compDiagonal" + t4c::integral_label(integral) + integral.label();
        
    return {label.size() + 1, label};
}

std::pair<size_t, std::string>
full_compute_func_name(const I4CIntegral& integral)
{
    const auto label = "compFull" + t4c::integral_label(integral) + integral.label();
        
    return {label.size() + 1, label};
}

std::string
integrand_label(const Operator& integrand)
{
    const auto iname = integrand.name();
    
    return iname;
}

std::string
namespace_label(const I4CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "1/|r-r'|")
    {
        return "erirec";
    }
    
    return std::string();
}

std::pair<size_t, std::string>
prim_diag_compute_func_name(const T4CIntegral& component,
                            const I4CIntegral& integral)
{
    auto label = "compPrimitiveDiag" + t4c::integral_label(integral) + integral.label();
    
    label += "_" + fstr::upcase(component.label());
    
    return {label.size() + 1, label};
}

std::pair<size_t, std::string>
prim_full_compute_func_name(const T4CIntegral& component,
                            const I4CIntegral& integral)
{
    auto label = "compPrimitiveFull" + t4c::integral_label(integral) + integral.label();
    
    if ((integral[0] + integral[1] + integral[2] + integral[3]) > 0)
    {
        label += "_" + fstr::upcase(component.label());
    }
    
    return {label.size() + 1, label};
}

std::pair<size_t, std::string>
prim_vrr_compute_func_name(const T4CIntegral& component,
                           const I4CIntegral& integral)
{
    auto label = "compPrimitiveVRR" + t4c::integral_label(integral) + integral.label();
    
    if ((integral[0] + integral[1] + integral[2] + integral[3]) > 0)
    {
        label += "_" + fstr::upcase(component.label());
    }
    
    return {label.size() + 1, label};
}

std::string
diag_prim_file_name(const T4CIntegral& component,
               const I4CIntegral& integral)
{
    auto label = "PrimitiveDiag" + t4c::integral_label(integral) + integral.label();
    
    if ((integral[0] + integral[1] + integral[2] + integral[3]) > 0)
    {
        label += "_" + fstr::upcase(component.label());
    }
    
    return label;
}

std::string
full_prim_file_name(const T4CIntegral& component,
                    const I4CIntegral& integral)
{
    auto label = "PrimitiveFull" + t4c::integral_label(integral) + integral.label();
    
    if ((integral[0] + integral[1] + integral[2] + integral[3]) > 0)
    {
        label += "_" + fstr::upcase(component.label());
    }
    
    return label;
}

std::string
full_vrr_file_name(const T4CIntegral& component,
                   const I4CIntegral& integral)
{
    auto label = "PrimitiveVRR" + t4c::integral_label(integral) + integral.label();
    
    if ((integral[0] + integral[1] + integral[2] + integral[3]) > 0)
    {
        label += "_" + fstr::upcase(component.label());
    }
    
    return label;
}

int
boys_order(const I4CIntegral& integral)
{
    const auto order = integral[0] + integral[1] + integral[2] + integral[3];
        
    return order;
}

std::string
get_factor_label(const R4CTerm&     rterm,
                 const I4CIntegral& integral,
                 const bool         first,
                 const bool         diagonal)
{
    auto mterm = R4CTerm(rterm);
    
    if (diagonal)
    {
        mterm.scale(Fraction(1, 2 * integral.order() + 1));
    }
    
    if (mterm.prefactor() != Fraction(0))
    {
        const auto pre_fact = mterm.prefactor();
        
        auto plabel = pre_fact.label();
        
        if (plabel == "1.0")  plabel = "";
        
        if (plabel == "-1.0") plabel = "-";
        
        if (pre_fact.denominator() != 1)
        {
            if (pre_fact.numerator() < 0) plabel.erase(0, 1);
            
            plabel = "(" + plabel + ")";
            
            if (pre_fact.numerator() < 0) plabel = "-" + plabel;
        }
        
        const auto facts = mterm.factors();
        
        std::string flabel;
        
        for (const auto& fact : facts)
        {
            const auto norder = mterm.factor_order(fact);
            
            for (size_t n = 0; n < norder; n++)
            {
                flabel += " * " + fact.label();
            }
        }
        
        // remove multiplication for special cases
        
        if ((pre_fact == Fraction(1)) || (pre_fact == Fraction(-1)))
        {
            flabel.erase(0, 3);
        }
        
        // merge labels
        
        flabel = plabel + flabel;
        
        if (!first)
        {
            if (flabel[0] == '-')
            {
                flabel.insert(1, " ");
            }
            else
            {
                flabel = "+ " + flabel;
            }
            
            flabel = " " + flabel;
        }
        
        return flabel;
    }
    else
    {
        return std::string();
    }
}

bool
find_factor(const R4CDist&     rdist,
            const std::string& label)
{
    for (const auto& fact : rdist.factors())
    {
        if (fact.label() == label) return true;
    }
    
    return false;
}

void
debug_info(const R4CDist& rdist)
{
    std::cout << "*** RECURSION FOR INTEGRAL COMPONENT: " << rdist.root().label() << std::endl;
        
    std::cout << " NUMBER OF TERMS:" << rdist.terms() << std::endl;
        
    for (size_t i = 0; i < rdist.terms(); i++)
    {
        std::cout << " RECURSION TERM (" << i << "): " << rdist[i].integral().bra().to_string() << " : ";
        
        std::cout << rdist[i].integral().ket().to_string() << " (" << rdist[i].order() << ") -> Factors: ";
        
        for (const auto& fact : rdist[i].factors())
        {
            std::cout << fact.label() << "  ";
        }
        
        std::cout << std::endl;
    }
        
    std::cout << std::endl << std::endl;
}

} // t4c namespace
