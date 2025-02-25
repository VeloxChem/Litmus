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

#include "t3c_utils.hpp"

#include "string_formater.hpp"
#include "v3i_eri_driver.hpp"
#include "v3i_geom010_eri_driver.hpp"

namespace t3c { // t3c namespace

std::string
integral_label(const I3CIntegral& integral)
{
    std::string suffix = "Geom";
    
    const auto prefixes = integral.prefixes();
    
    if (const auto nterms = prefixes.size(); nterms == 3)
    {
        suffix += std::to_string(prefixes[0].shape().order());
        
        suffix += std::to_string(prefixes[1].shape().order());
        
        suffix += std::to_string(prefixes[2].shape().order());
    }
    
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "1/|r-r'|")
    {
        return (prefixes.empty()) ? "ThreeCenterElectronRepulsion" : "ThreeCenterElectronRepulsion" + suffix;
    }
    
    return std::string();
}

std::string
integral_split_label(const I3CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "1/|r-r'|")
    {
        return "Electron_Repulsion";
    }
    
    return std::string();
}

std::string
prim_file_name(const I3CIntegral& integral)
{
    return t3c::integral_label(integral) + "PrimRec" + integral.label();
}

std::string
hrr_file_name(const I3CIntegral& integral)
{
    const auto ket_one = Tensor(integral[1]);
    
    const auto ket_two = Tensor(integral[2]);
    
    return t3c::integral_label(integral) + "ContrRecX" + ket_one.label() + ket_two.label();
}

std::string
namespace_label(const I3CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "1/|r-r'|")
    {
        return "t3ceri";
    }
    
    return std::string();
}

std::string
integrand_label(const Operator& integrand)
{
    const auto iname = integrand.name();
    
    // TODO: Fix other stuff.
    
    return iname;
}

std::string
compute_func_name(const I3CIntegral& integral)
{
    std::string geom_label;
    
    auto tint_prefixes = integral.prefixes();
    
    if (!tint_prefixes.empty())
    {
        geom_label += "_geom";
        
        for (const auto& tint_prefix : tint_prefixes)
        {
            geom_label += std::to_string(tint_prefix.shape().order());
        }
    }
    
    auto label = "comp_" + t3c::integral_split_label(integral) +  geom_label + "_" + integral.label();
        
    return fstr::lowercase(label);
}

std::string
get_buffer_label(const I3CIntegral& integral,
                 const std::string& prefix)
{
    std::string label = prefix + "_buffer_";
    
    label += std::to_string(integral.order()) + "_";
    
    auto tint_prefixes = integral.prefixes();
    
    if (!tint_prefixes.empty())
    {
        label += "geom";
        
        for (const auto& tint_prefix : tint_prefixes)
        {
            label += std::to_string(tint_prefix.shape().order());
        }
        
       label += "_";
    }
    
    label += fstr::lowercase(integral.label());

    return label;
}

std::string
prim_compute_func_name(const I3CIntegral& integral)
{
    auto label =  "comp_prim_" + t3c::integral_split_label(integral) + "_" + integral.label();
    
    return fstr::lowercase(label);
}

SI3CIntegrals
get_vrr_integrals(const I3CIntegral& integral)
{
    SI3CIntegrals tints;
    
    if (integral.integrand().name() == "1/|r-r'|")
    {
        V3IElectronRepulsionDriver eri_drv;
        
        if (integral[0] > 0)
        {
            tints = eri_drv.bra_vrr(integral);
        }
        else
        {
            tints = eri_drv.ket_vrr(integral);
        }
    }
    
    return tints;
}

std::string
hrr_compute_func_name(const I3CIntegral& integral)
{
    const auto ket_one = Tensor(integral[1]);
    
    const auto ket_two = Tensor(integral[2]);
    
    auto label =  "comp_hrr_" + t3c::integral_split_label(integral) + "_x" + ket_one.label() + ket_two.label();
    
    return fstr::lowercase(label);
}

SI3CIntegrals
get_hrr_integrals(const I3CIntegral& integral)
{
    SI3CIntegrals tints;
    
    if (integral.integrand().name() == "1/|r-r'|")
    {
        V3IElectronRepulsionDriver eri_drv;
        
        if (integral[1] > 0)
        {
            tints = eri_drv.ket_hrr(integral);
        }
    }
    
    return tints;
}

SI3CIntegrals
get_geom_hrr_integrals(const I3CIntegral& integral)
{
    SI3CIntegrals tints;
    
    if (integral.integrand().name() == "1/|r-r'|")
    {
        V3IGeom010ElectronRepulsionDriver eri_drv;
        
        if (integral[1] == 0)
        {
            tints = eri_drv.ket_aux_hrr(integral);
        }
        else
        {
            tints = eri_drv.ket_hrr(integral);
        }
    }
    
    return tints;
}

std::string
get_index_label(const I3CIntegral& integral)
{
    const auto prefixes = integral.prefixes();
    
    std::string geom_label;
    
    if (!prefixes.empty())
    {
        geom_label = "geom_";
        
        for (const auto& prefix : prefixes)
        {
            geom_label += std::to_string(prefix.shape().order()) + "0";
        }
    }
    
    std::string label = "idx_";
    
    if (integral.integrand().name() == "1/|r-r'|")
    {
        label += "eri_" + std::to_string(integral.order()) + "_";
    }
    
    label += fstr::lowercase(integral.label());

    return label;
}

std::string
get_hrr_index(const I3CIntegral& integral)
{
    std::string label = "idx_";
    
    if (const auto geom_order = integral.prefixes_order(); !geom_order.empty())
    {
        label += "geom_" + std::to_string(geom_order[1]) + std::to_string(geom_order[2]) + "_";
    }
    
    const auto ket_one = Tensor(integral[1]);
        
    const auto ket_two = Tensor(integral[2]);
        
    label += "x" + ket_one.label() + ket_two.label();
    
    return fstr::lowercase(label);
}

std::string
get_hrr_buffer_label(const I3CIntegral& integral)
{
    std::string label = "contr_buffer_";
    
    const auto ket_one = Tensor(integral[1]);
        
    const auto ket_two = Tensor(integral[2]);
        
    label += "x" + ket_one.label() + ket_two.label();
    
    return fstr::lowercase(label);
}

std::string
prefixes_label(const I3CIntegral& integral)
{
    const auto prefixes = integral.prefixes();
    
    std::string label;
    
    if (const auto border = prefixes[0].shape().order(); border > 0)
    {
        label += "d^(" + std::to_string(border) + ")/dA^(" + std::to_string(border) + ")";
    }
    
    if (const auto border = prefixes[1].shape().order(); border > 0)
    {
        label += "d^(" + std::to_string(border) + ")/dC^(" + std::to_string(border) + ")";
    }
    
    if (const auto border = prefixes[2].shape().order(); border > 0)
    {
        label += "d^(" + std::to_string(border) + ")/dD^(" + std::to_string(border) + ")";
    }
    
    return label;
}

G3Term
prune_term(const G3Term& term)
{
    const auto tint = term.second;
    
    if (tint.prefixes_order() == std::vector<int>({1, 0, 0}) && (tint[0] == 0))
    {
        const auto cint = *tint.shift(1, 0);
        
        return G3Term({std::array<int, 3>({1, 0, 0}), cint.base()});
    }
    
    if (tint.prefixes_order() == std::vector<int>({0, 1, 0}) && (tint[1] == 0))
    {
        const auto cint = *tint.shift(1, 1);
        
        return G3Term({std::array<int, 3>({0, 1, 0}), cint.base()});
    }
    
    return term;
}

std::string
bra_geom_compute_func_name(const I3CIntegral& integral)
{
    const auto bra_one = Tensor(integral[0]);
    
    auto geom_orders = integral.prefixes_order();
    
    auto label =  "comp_bra_geom" + std::to_string(geom_orders[0]);
    
    if ((geom_orders[1] + geom_orders[2]) > 0)
    {
        label += std::to_string(geom_orders[1]) + std::to_string(geom_orders[2]);
    }
    
    label += "_" + t3c::integral_split_label(integral);
    
    label += "_" + bra_one.label() + "xx";
    
    return fstr::lowercase(label);
}

std::string
ket_geom_compute_func_name(const I3CIntegral& integral)
{
    const auto ket_one = Tensor(integral[1]);
    
    const auto ket_two = Tensor(integral[2]);
    
    auto geom_orders = integral.prefixes_order();
    
    auto label =  "comp_ket_geom" + std::to_string(geom_orders[0]);
    
    if ((geom_orders[1] + geom_orders[2]) > 0)
    {
        label += std::to_string(geom_orders[1]) + std::to_string(geom_orders[2]);
    }
    
    label += "_" + t3c::integral_split_label(integral);
    
    label += "_x" + ket_one.label() + ket_two.label();
    
    return fstr::lowercase(label);
}

SI3CIntegrals
get_bra_geom_integrals(const I3CIntegral& integral)
{
    const auto geom_order = integral.prefixes_order();
    
    SI3CIntegrals tints;
    
    if (geom_order == std::vector<int>({1, 0, 0}))
    {
        if (const auto tval = integral.shift(1, 0))
        {
            tints.insert(tval->base());
        }
        
        if (const auto tval = integral.shift(-1, 0))
        {
            tints.insert(tval->base());
        }
    }
    
    return tints;
}

std::string
bra_geom_file_name(const I3CIntegral& integral)
{
    const auto bra_one = Tensor(integral[0]);
    
    return t3c::integral_label(integral) + "ContrRec" + bra_one.label() + "XX";
}

std::string
ket_geom_file_name(const I3CIntegral& integral)
{
    const auto ket_one = Tensor(integral[1]);
    
    const auto ket_two = Tensor(integral[2]);
    
    return t3c::integral_label(integral) + "ContrRecX" + ket_one.label() + ket_two.label();
}

std::string
get_full_hrr_index(const I3CIntegral& integral,
                   const bool         use_ket)
{
    std::string label = "idx_";
    
    if (const auto geom_order = integral.prefixes_order(); !geom_order.empty())
    {
        label += "geom_" + std::to_string(geom_order[0]) + std::to_string(geom_order[1]);
      
        label += std::to_string(geom_order[2]) + "_";
    }
    
    if (use_ket)
    {
        const auto ket_one = Tensor(integral[1]);
        
        const auto ket_two = Tensor(integral[2]);
        
        label += "x" + ket_one.label() + ket_two.label();
    }
    else
    {
        const auto bra_one = Tensor(integral[0]);
        
        label += bra_one.label()  + "xx";
    }
    
    return fstr::lowercase(label);
}

std::string
get_hrr_buffer_label(const I3CIntegral& integral,
                     const bool         use_ket)
{
    std::string label = "contr_buffer_";
    
    if (use_ket)
    {
        const auto ket_one = Tensor(integral[1]);
        
        const auto ket_two = Tensor(integral[2]);
        
        label += "xx" + ket_one.label() + ket_two.label();
    }
    else
    {
        const auto bra_one = Tensor(integral[0]);
        
        label += bra_one.label() + "xx";
    }
    
    return fstr::lowercase(label);
}


} // t3c namespace
