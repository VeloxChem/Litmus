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

#include <iostream>

#include "string_formater.hpp"
#include "v2i_ovl_driver.hpp"

namespace t2c { // t2c namespace

std::string
fraction_label(const Fraction& fraction)
{
    return "fr_" + std::to_string(std::abs(fraction.numerator())) + "_" + std::to_string(fraction.denominator());
}

std::string
integral_label(const I2CIntegral& integral)
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
    
    if (integrand.name() == "AG")
    {
        return "NuclearPotentialGeom" + suffix;
    }
    
    if (integrand.name() == "A")
    {
        return (prefixes.empty()) ? "NuclearPotential" : "NuclearPotential" + suffix;
    }
    
    if (integrand.name() == "T")
    {
        return (prefixes.empty()) ? "KineticEnergy" : "KineticEnergy" + suffix;
    }
    
    if (integrand.name() == "1")
    {
        return (prefixes.empty()) ? "Overlap" : "Overlap" + suffix;
    }
    
    if (integrand.name() == "r")
    {
        const auto iorder = integrand.shape().order();
        
        if (iorder == 1)
        {
            return "Dipole";
        }
        
        if (iorder == 2)
        {
            return "Quadrupole";
        }
        
        if (iorder == 3)
        {
            return "Octupole";
        }
    }
    
    if (integrand.name() == "G(r)")
    {
        const auto iorder = integrand.shape().order();
        
        if (iorder == 0)
        {
            return "ThreeCenterOverlap";
        }
        
        if (iorder == 1)
        {
            return "ThreeCenterOverlapGradient";
        }
    }
    
    return std::string();
}

std::string
integrand_label(const Operator& integrand)
{
    const auto iname = integrand.name();
    
    const auto iorder = std::to_string(integrand.shape().order());
    
    if (iname == "AG")
    {
        return iname + "(" + iorder + ")";
    }
    
    if ((iname == "r") && (iorder != "1"))
    {
        return iname + "^" + iorder;
    }
    
    return iname;
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
integrand_components(const Tensor&      bra_tensor,
                     const Operator&    integrand,
                     const std::string& label)
{
    if (const auto icomps = integrand.components(); icomps.size() == 1)
    {
        return t2c::tensor_components(bra_tensor, label);
    }
    else
    {
        std::vector<std::string> labels;
        
        for (const auto& bcomp : bra_tensor.components())
        {
            for (const auto& icomp : icomps)
            {
                labels.push_back(label + "_" + bcomp.label() + "_" + icomp.label());
            }
        }
            
        return labels;
    }
}

std::vector<std::string>
integrand_components(const Tensor&      bra_tensor,
                     const Tensor&      ket_tensor,
                     const Operator&    integrand,
                     const std::string& label)
{
    if (const auto icomps = integrand.components(); icomps.size() == 1)
    {
        std::vector<std::string> labels;
        
        for (const auto& bcomp : bra_tensor.components())
        {
            for (const auto& kcomp : ket_tensor.components())
            {
                labels.push_back(label + "_" + bcomp.label() + "_" + kcomp.label());
            }
        }
            
        return labels;
    }
    else
    {
        std::vector<std::string> labels;
        
        for (const auto& bcomp : bra_tensor.components())
        {
            for (const auto& kcomp : ket_tensor.components())
            {
                for (const auto& icomp : icomps)
                {
                    labels.push_back(label + "_" + bcomp.label() + "_" + kcomp.label() + "_" + icomp.label());
                }
            }
        }
            
        return labels;
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
compute_func_name(const I2CIntegral& integral,
                  const bool         sum_form)
{
    std::string prefix = (sum_form) ? "comp_sum_" : "comp_";
    
    auto label = prefix  + t2c::integral_label(integral) + "_" + integral.label();
    
    label = fstr::lowercase(label); 
        
    return {label.size() + 1, label};
}

std::pair<size_t, std::string>
auxilary_func_name(const I2CIntegral& integral)
{
    std::string prefix = "compAuxilary";
    
    const auto label = prefix  + t2c::integral_label(integral) + integral.label();
        
    return {label.size() + 1, label};
}

std::pair<size_t, std::string>
prim_compute_func_name(const I2CIntegral& integral,
                       const bool         sum_form)
{
    std::string prefix = (sum_form) ? "comp_sum_" : "comp_";
    
    auto label =  prefix + "prim_" + t2c::integral_label(integral) + "_" + integral.label();
    
    label = fstr::lowercase(label);
    
    return {label.size() + 1, label};
}

std::pair<size_t, std::string>
prim_compute_func_name(const TensorComponent& component,
                       const I2CIntegral&     integral,
                       const bool             sum_form,
                       const bool             bra_first)
{
    std::string prefix = (sum_form) ? "compSum" : "comp";
    
    auto label = prefix + "Primitive" + t2c::integral_label(integral) + integral.label();
    
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
                       const I2CIntegral&     integral,
                       const bool             sum_form)
{
    std::string prefix = (sum_form) ? "compSum" : "comp";
    
    auto label = prefix + "Primitive" + t2c::integral_label(integral) + integral.label();
    
    label += "_" + fstr::upcase(bra_component.label());
   
    label += "_" + fstr::upcase(ket_component.label());
    
    return {label.size() + 1, label};
}

std::string
auxilary_file_name(const I2CIntegral& integral)
{
    return t2c::integral_label(integral) + "Auxilary" + integral.label();
}

std::string
prim_file_name(const I2CIntegral& integral,
               const bool         sum_form)
{
    std::string prefix = (sum_form) ? "Sum" : "";
    
    return prefix + "Primitive" + t2c::integral_label(integral) + integral.label();
}

std::string
prim_file_name(const I2CIntegral& integral)
{
    return t2c::integral_label(integral) + "PrimRec" + integral.label();
}

std::string
prim_file_name(const TensorComponent& component,
               const I2CIntegral&     integral,
               const bool             sum_form,
               const bool             bra_first)
{
    std::string prefix = (sum_form) ? "Sum" : "";
    
    auto label = prefix + "Primitive" + t2c::integral_label(integral) + integral.label();
    
    if (bra_first)
    {
        label += "_" + fstr::upcase(component.label()) + "_T";
    }
    else
    {
        label += "_T_" + fstr::upcase(component.label());
    }
    
    return label;
}

std::string
prim_file_name(const TensorComponent& bra_component,
               const TensorComponent& ket_component,
               const I2CIntegral&     integral,
               const bool             sum_form)
{
    std::string prefix = (sum_form) ? "Sum" : "";
    
    auto label = prefix + "Primitive" + t2c::integral_label(integral) + integral.label();
    
    label += "_" + fstr::upcase(bra_component.label());
   
    label += "_" + fstr::upcase(ket_component.label());
    
    return label;
}

std::string
namespace_label(const I2CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    const auto iorder = std::to_string(integrand.shape().order());
    
    if (integrand.name() == "AG")
    {
        return "npotg0" + iorder + "0rec";
    }
    
    if (integrand.name() == "A")
    {
        return "npotrec";
    }
    
    if (integrand.name() == "T")
    {
        return "kinrec";
    }
    
    if (integrand.name() == "1")
    {
        return "ovlrec";
    }
    
    if (integrand.name() == "r")
    {
        if (iorder == "1") return "diprec";
        
        if (iorder == "2") return "quadrec";
        
        if (iorder == "3") return "octurec";
    }
    
    if (integrand.name() == "G(r)")
    {
        if (iorder == "0") return "t3ovlrec";
        
        if (iorder == "1") return "g3ovlrec";
    }
    
    return std::string();
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

bool
find_factor(const R2Group&     rgroup,
            const std::string& label)
{
    for (const auto& fact : rgroup.factors())
    {
        if (fact.label() == label) return true;
    }
    
    return false;
}

bool
find_factor(const R2Group&     rgroup,
            const std::string& label,
            const size_t       first,
            const size_t       last)
{
    for (size_t i = first; i < last; i++)
    {
        for (const auto& fact : rgroup[i].factors())
        {
            if (fact.label() == label) return true;
        }
    }
    
    return false;
}

std::string
get_factor_label(const R2CTerm& rterm,
                 const bool     first)
{
    const auto pre_fact = rterm.prefactor();
        
    auto plabel = pre_fact.label();
        
    if (plabel == "1.0")  plabel = "";
        
    if (plabel == "-1.0") plabel = "-";
        
    if (pre_fact.denominator() != 1)
    {
        if (pre_fact.numerator() < 0) plabel.erase(0, 1);
            
        plabel = "(" + plabel + ")";
            
        if (pre_fact.numerator() < 0) plabel = "-" + plabel;
    }
        
    const auto facts = rterm.factors();
        
    std::string flabel;
        
    for (const auto& fact : facts)
    {
        const auto norder = rterm.factor_order(fact);
            
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
    
    // fix for special case
    
    if (flabel == "-") flabel = "-1.0";
        
    return flabel;
}

int
boys_order(const I2CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    const auto order = integral[0] + integral[1];
    
    if (integrand.name() == "AG")
    {
        return order + integrand.shape().order();
    }
    
    if (integrand.name() == "A")
    {
        if (integral.is_simple())
        {
            return order;
        }
        else
        {
            auto morder = order;
            
            for (const auto& prefix : integral.prefixes())
            {
                morder += prefix.shape().order();
            }
            
            return morder;
        }
    }
    
    return -1;
}

bool
need_boys(const I2CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "AG")
    {
        return true;
    }
    
    if (integrand.name() == "A")
    {
        return true;
    }
    
    return false;
}

V4Auxilaries
get_unique_auxilaries(const R2Group& rgroup)
{
    V4Auxilaries auxs;
    
    for (size_t i = 0; i < rgroup.expansions(); i++)
    {
        for (size_t j = 0; j < rgroup[i].terms(); j++)
        {
            auxs.insert(t2c::get_auxilary(rgroup[i][j]));
        }
    }
   
    return auxs;
}

V4Auxilaries
get_unique_auxilaries(const R2Group& rgroup,
                      const size_t   first,
                      const size_t   last)
{
    V4Auxilaries auxs;
    
    for (size_t i = first; i < last; i++)
    {
        for (size_t j = 0; j < rgroup[i].terms(); j++)
        {
            auxs.insert(t2c::get_auxilary(rgroup[i][j]));
        }
    }
   
    return auxs;
}

size_t
get_auxilary_index(const V4Auxilaries& auxilaries,
                   const T4Index&      target)
{
    size_t index = 0;
    
    for (const auto& taux : auxilaries)
    {
        if (taux == target) return index;
        
        index++;
    }
   
    return -1;
}

T4Index
get_auxilary(const R2CTerm& rterm)
{
    const auto n = rterm.factor_order(Factor("N", "n"));
    
    const auto m = rterm.factor_order(Factor("M", "m"));
    
    const auto t = rterm.factor_order(Factor("T", "t"));
    
    const auto p = rterm.order();
    
    return {n, m, t, p};
}

T3Index
get_factor_decomposition(const T4Index& target)
{
    auto ft = std::min(target[0], std::min(target[1], target[2]));
    
    auto fm = std::min(target[0] - ft, target[2] - ft);
        
    auto fn = std::min(target[1] - ft, target[2] - ft - fm);
    
    return T3Index({ft, fm, fn});
}

T3Index
get_maximum_decomposition(const V4Auxilaries& auxilaries)
{
    auto mvals = T3Index({0, 0, 0});
    
    for (const auto& taux : auxilaries)
    {
        auto t3vals = t2c::get_factor_decomposition(taux);
        
        for (size_t i = 0; i < 3; i++)
        {
            mvals[i] = std::max(mvals[i], t3vals[i]);
        }
    }
    
    return mvals;
}

void
debug_info(const R2CDist& rdist)
{
    std::cout << "*** RECURSION FOR INTEGRAL COMPONENT: " << rdist.root().label() << std::endl;
    
    std::cout << " NUMBER OF TERMS:" << rdist.terms() << std::endl;
    
    for (size_t i = 0; i < rdist.terms(); i++)
    {
       // std::cout << " RECURSION TERM (" << i << "): " << rdist[i].label() << std::endl;
        
        std::cout << " RECURSION TERM (" << i << "): " << rdist[i].integral().bra().to_string() << " : "  << rdist[i].integral().ket().to_string() << std::endl;
    }
    
    std::cout << std::endl << std::endl;
}

SI2CIntegrals
get_integrals(const I2CIntegral& integral)
{
    SI2CIntegrals tints;
    
    if (integral.integrand().name() == "1")
    {
        V2IOverlapDriver ovl_drv;
        
        if (integral[0] > 0)
        {
            tints = ovl_drv.bra_vrr(integral);
        }
        else
        {
            tints = ovl_drv.ket_vrr(integral);
        }
    }
    
    return tints;
}

std::string
get_buffer_label(const I2CIntegral& integral,
                 const std::string& prefix)
{
    std::string label = prefix + "_buffer_";
    
    if (integral.integrand().name() == "1") label += "ovl_";
    
    label += fstr::lowercase(integral.label());

    return label;
}

} // t2c namespace
