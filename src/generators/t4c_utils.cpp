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

#include "string_formater.hpp"

#include "v4i_eri_driver.hpp"
#include "t4c_center_driver.hpp"
#include "v4i_geom10_eri_driver.hpp"

namespace t4c { // t4c namespace

std::string
integral_label(const I4CIntegral& integral)
{
    std::string suffix = "Geom";
    
    const auto prefixes = integral.prefixes();
    
    if (const auto nterms = prefixes.size(); nterms == 4)
    {
        suffix += std::to_string(prefixes[0].shape().order());
        
        suffix += std::to_string(prefixes[1].shape().order());
        
        suffix += std::to_string(prefixes[2].shape().order());
        
        suffix += std::to_string(prefixes[3].shape().order());
    }
    
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "1/|r-r'|")
    {
        return (prefixes.empty()) ? "ElectronRepulsion" : "ElectronRepulsion" + suffix;
    }
    
    return std::string(); 
}

std::string
integral_split_label(const I4CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "1/|r-r'|")
    {
        return "Electron_Repulsion";
    }
    
    return std::string();
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

std::string
geom_namespace_label()
{
    return std::string("t4c_geom");
}

std::string
integrand_label(const Operator& integrand)
{
    const auto iname = integrand.name();
    
    // TODO: Fix other stuff.
    
    return iname;
}

std::string
compute_func_name(const I4CIntegral& integral)
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
    
    auto label = "comp_" + t4c::integral_split_label(integral) +  geom_label + "_" + integral.label();
        
    return fstr::lowercase(label);
}

std::string
diag_compute_func_name(const I4CIntegral& integral)
{
    std::string geom_label;
    
    auto label = "comp_diag_" + t4c::integral_split_label(integral) + "_" + integral.label();
        
    return fstr::lowercase(label);
}

std::string
get_buffer_label(const I4CIntegral& integral,
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
get_geom_buffer_label(const I4CIntegral& integral)
{
    std::string label = "buffer_";
    
    auto tint_prefixes = integral.prefixes();
    
    if (!tint_prefixes.empty())
    {
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
get_hrr_buffer_label(const I4CIntegral& integral,
                     const bool         use_ket)
{
    std::string label = "contr_buffer_";
    
    if (use_ket)
    {
        const auto ket_one = Tensor(integral[2]);
        
        const auto ket_two = Tensor(integral[3]);
        
        label += "xx" + ket_one.label() + ket_two.label();
    }
    else
    {
        const auto bra_one = Tensor(integral[0]);
        
        const auto bra_two = Tensor(integral[1]);
        
        label += bra_one.label() + bra_two.label() + "xx"; 
    }
    
    return fstr::lowercase(label);
}

std::string
prim_compute_func_name(const I4CIntegral& integral)
{
    auto label =  "comp_prim_" + t4c::integral_split_label(integral) + "_" + integral.label();
    
    return fstr::lowercase(label);
}

std::string
geom_compute_func_name(const I4CIntegral& integral)
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
    
    auto label =  "comp" + geom_label + "_" + integral.label() + "_" +  std::to_string(integral.integrand().shape().order());
    
    return fstr::lowercase(label);
}

std::string
ket_hrr_compute_func_name(const I4CIntegral& integral)
{
    const auto ket_one = Tensor(integral[2]);
    
    const auto ket_two = Tensor(integral[3]);
    
    auto label =  "comp_ket_hrr_" + t4c::integral_split_label(integral) + "_xx" + ket_one.label() + ket_two.label();
    
    return fstr::lowercase(label);
}

std::string
bra_hrr_compute_func_name(const I4CIntegral& integral)
{
    const auto bra_one = Tensor(integral[0]);
    
    const auto bra_two = Tensor(integral[1]);
    
    auto label =  "comp_bra_hrr_" + t4c::integral_split_label(integral) + "_" + bra_one.label() + bra_two.label() + "xx";
    
    return fstr::lowercase(label);
}

std::string
bra_geom_hrr_compute_func_name(const I4CIntegral& integral)
{
    const auto bra_one = Tensor(integral[0]);
    
    const auto bra_two = Tensor(integral[1]);
    
    auto geom_orders = integral.prefixes_order();
    
    auto label =  "comp_bra_geom" + std::to_string(geom_orders[0]) + std::to_string(geom_orders[1]);
    
    label += "_hrr_" + t4c::integral_split_label(integral) + "_" + bra_one.label() + bra_two.label() + "xx";
    
    return fstr::lowercase(label);
}

SI4CIntegrals
get_vrr_integrals(const I4CIntegral& integral)
{
    SI4CIntegrals tints;
    
    if (integral.integrand().name() == "1/|r-r'|")
    {
        V4IElectronRepulsionDriver eri_drv;
        
        if (integral[1] > 0)
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

SI4CIntegrals
get_full_vrr_integrals(const I4CIntegral& integral)
{
    SI4CIntegrals tints;
    
    if (integral.integrand().name() == "1/|r-r'|")
    {
        V4IElectronRepulsionDriver eri_drv;
        
        if (integral[0] > 0)
        {
            tints = eri_drv.bra_vrr_a(integral);
        }
        
        if ((integral[1] > 0) && (integral[0] == 0))
        {
            tints = eri_drv.bra_vrr_b(integral);
        }
        
        if ((integral[2] > 0) && ((integral[0] + integral[1]) == 0))
        {
            tints = eri_drv.ket_vrr_c(integral);
        }
        
        if ((integral[3] > 0) && ((integral[0] + integral[1] + integral[2]) == 0))
        {
            tints = eri_drv.ket_vrr_d(integral);
        }
    }
    
    return tints;
}

SI4CIntegrals
get_ket_hrr_integrals(const I4CIntegral& integral)
{
    SI4CIntegrals tints;
    
    if (integral.integrand().name() == "1/|r-r'|")
    {
        V4IElectronRepulsionDriver eri_drv;
        
        if (integral[2] > 0)
        {
            tints = eri_drv.ket_hrr(integral);
        }
    }
    
    return tints;
}

SI4CIntegrals
get_bra_hrr_integrals(const I4CIntegral& integral)
{
    SI4CIntegrals tints;
    
    if (integral.integrand().name() == "1/|r-r'|")
    {
        V4IElectronRepulsionDriver eri_drv;
        
        if (integral[0] > 0)
        {
            tints = eri_drv.bra_hrr(integral);
        }
    }
    
    return tints;
}

SI4CIntegrals
get_bra_geom_hrr_integrals(const I4CIntegral& integral)
{
    const auto geom_order = integral.prefixes_order();
    
    SI4CIntegrals tints;
    
    if (geom_order == std::vector<int>({1, 0, 0, 0}))
    {
        V4IGeom10ElectronRepulsionDriver geom_drv;
        
        tints = geom_drv.bra_hrr(integral);
    }
    
    return tints;
}

SI4CIntegrals
get_geom_integrals(const I4CIntegral& integral)
{
    R4Group rgroup;
        
    T4CCenterDriver t4c_geom_drv;
        
    rgroup = t4c_geom_drv.create_recursion(integral.components<T2CPair, T2CPair>());
    
    SI4CIntegrals tints;
    
    for (size_t i = 0; i < rgroup.expansions(); i++)
    {
        for (size_t j = 0; j < rgroup[i].terms(); j++)
        {
            tints.insert(I4CIntegral(rgroup[i][j].integral().base()));
        }
    }
    
    return tints;
}

std::string
prim_file_name(const I4CIntegral& integral)
{
    return t4c::integral_label(integral) + "PrimRec" + integral.label();
}

std::string
geom_file_name(const I4CIntegral& integral)
{
    std::string label = "GeomDeriv";
    
    for (const auto& prefix : integral.prefixes())
    {
        label += std::to_string(prefix.shape().order());
    }
    
    if (integral.integrand().shape().order() == 0)
    {
        label += "OfScalar";
    }
    
    if (integral.integrand().shape().order() == 1)
    {
        label += "OfVector";
    }
    
    label += "For" + integral.label();
    
    return label;
}

std::string
ket_hrr_file_name(const I4CIntegral& integral)
{
    const auto ket_one = Tensor(integral[2]);
    
    const auto ket_two = Tensor(integral[3]);
    
    return t4c::integral_label(integral) + "ContrRecXX" + ket_one.label() + ket_two.label();
}

std::string
bra_hrr_file_name(const I4CIntegral& integral)
{
    const auto bra_one = Tensor(integral[0]);
    
    const auto bra_two = Tensor(integral[1]);
    
    return t4c::integral_label(integral) + "ContrRec" + bra_one.label() + bra_two.label() + "XX";
}

std::string
bra_geom_hrr_file_name(const I4CIntegral& integral)
{
    const auto bra_one = Tensor(integral[0]);
    
    const auto bra_two = Tensor(integral[1]);
        
    return t4c::integral_label(integral) + "ContrRec" + bra_one.label() + bra_two.label() + "XX";
}

// Only relevant for geometric derivatives
std::string
prefixes_label(const I4CIntegral& integral)
{
    const auto prefixes = integral.prefixes();
    
    std::string label;
    
    if (const auto border = prefixes[0].shape().order(); border > 0)
    {
        label += "d^(" + std::to_string(border) + ")/dA^(" + std::to_string(border) + ")";
    }
    
    if (const auto border = prefixes[1].shape().order(); border > 0)
    {
        label += "d^(" + std::to_string(border) + ")/dB^(" + std::to_string(border) + ")";
    }
    
    if (const auto border = prefixes[2].shape().order(); border > 0)
    {
        label += "d^(" + std::to_string(border) + ")/dC^(" + std::to_string(border) + ")";
    }
    
    if (const auto border = prefixes[3].shape().order(); border > 0)
    {
        label += "d^(" + std::to_string(border) + ")/dD^(" + std::to_string(border) + ")";
    }
    
    return label;
}

std::string
get_index_label(const I4CIntegral& integral)
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
get_hrr_index(const I4CIntegral& integral,
                     const bool         use_ket)
{
    std::string label = "idx_";
    
    if (use_ket)
    {
        const auto ket_one = Tensor(integral[2]);
        
        const auto ket_two = Tensor(integral[3]);
        
        label += "xx" + ket_one.label() + ket_two.label();
    }
    else
    {
        const auto bra_one = Tensor(integral[0]);
        
        const auto bra_two = Tensor(integral[1]);
        
        label += bra_one.label() + bra_two.label() + "xx";
    }
    
    return fstr::lowercase(label);
}

} // t4c namespace
