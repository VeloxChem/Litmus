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

#include "t2c_utils.hpp"

#include "string_formater.hpp"

namespace t2c { // t2c namespace

std::string
integral_label(const I2CIntegral& integral)
{
    if (integral.is_simple())
    {
        auto labels = TMapOfStrings({ {Operator("1"), "Overlap"},
                                      {Operator("T"), "KineticEnergy"},
        });
        
        return labels[integral.integrand()];
    }
    else
    {
        return std::string();
    }
}

std::string
integrand_label(const Operator& integrand)
{
    auto labels = TMapOfStrings({ {Operator("1"), ""},
                                  {Operator("T"), "T"},
                                 });
        
    return labels[integrand];
}

std::vector<std::string>
integrand_components(const Operator&    integrand,
                     const std::string& label)
{
    if (const auto icomps = integrand.components(); icomps.size() == 1)
    {
        return std::vector<std::string>({label,});
    }
    else
    {
        std::vector<std::string> ilabels;
            
        for (const auto& icomp : icomps)
        {
            ilabels.push_back(label + "_" + icomp.label());
        }
            
        return ilabels;
    }
}

std::vector<std::string>
tensor_components(const Tensor&      tensor,
                  const std::string& label)
{
    if (const auto tcomps = tensor.components(); tcomps.size() == 1)
    {
        return std::vector<std::string>({label,});
    }
    else
    {
        std::vector<std::string> tlabels;
            
        for (const auto& tcomp : tcomps)
        {
            tlabels.push_back(label + "_" + tcomp.label());
        }
            
        return tlabels;
    }
}

std::pair<size_t, std::string>
compute_func_name(const I2CIntegral& integral)
{
    const auto label = "comp" + t2c::integral_label(integral) + integral.label();
        
    return {label.size() + 1, label};
}

std::pair<size_t, std::string>
prim_compute_func_name(const I2CIntegral& integral)
{
    const auto label = "compPrimitive" + t2c::integral_label(integral) + integral.label();
    
    return {label.size() + 1, label};
}

std::pair<size_t, std::string>
prim_compute_func_name(const TensorComponent& component,
                       const I2CIntegral&     integral,
                       const bool             bra_first)
{
    auto label = "compPrimitive" + t2c::integral_label(integral) + integral.label();
    
    if (bra_first)
    {
        label += "_" + fstr::upcase(component.label()) + "_T";
    }
    else
    {
        label += "_T_" + fstr::upcase(component.label());
    }
    
    return {label.size() + 1, label};
}

std::pair<size_t, std::string>
prim_compute_func_name(const TensorComponent& bra_component,
                       const TensorComponent& ket_component,
                       const I2CIntegral&     integral)
{
    auto label = "compPrimitive" + t2c::integral_label(integral) + integral.label();
    
    label += "_" + fstr::upcase(bra_component.label());
   
    label += "_" + fstr::upcase(ket_component.label());
    
    return {label.size() + 1, label};
}

std::string
namespace_label(const I2CIntegral& integral)
{
    auto labels = TMapOfStrings({ {Operator("1"), "ovlrec"},
                                  {Operator("T"), "kinrec"},
                                });
    
    return labels[integral.integrand()];
}

int
tensor_component_index(const TensorComponent& component)
{
    const auto tcomps = Tensor(component).components();
        
    for (int i = 0; i < tcomps.size(); i++)
    {
            if (tcomps[i] == component) return  i;
    }
        
    return -1;
}

std::string
combine_factors(const std::string& bra_factor,
                const std::string& ket_factor)
{
    auto bra_label = bra_factor;
        
    auto ket_label = ket_factor;
        
    // get signs
        
    int bra_sign = (bra_factor[0] == '-') ? -1 : 1;
        
    int ket_sign = (ket_factor[0] == '-') ? -1 : 1;
        
    //  combinne symbolic factor
        
    if (bra_sign < 0) bra_label.erase(0, 1);
        
    if (ket_sign < 0) ket_label.erase(0, 1);
        
    std::string label;
        
    if (bra_label != "1.0") label = bra_label;
        
    if (ket_label != "1.0") label = (label.empty()) ? ket_label : label + " * " + ket_label;
        
    if (label.empty()) label = "1.0";
        
    if (bra_sign * ket_sign < 0) label.insert(0, "-");
        
    return label;
}

} // t2c namespace
