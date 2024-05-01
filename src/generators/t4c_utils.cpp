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

namespace t4c { // t4c namespace

std::string
integral_label(const I4CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "1/|r-r'|")
    {
        return "ElectronRepulsion";
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
integrand_label(const Operator& integrand)
{
    const auto iname = integrand.name();
    
    // TODO: Fix other stuff.
    
    return iname;
}

std::string
compute_func_name(const I4CIntegral& integral)
{
    auto label = "comp_" + t4c::integral_split_label(integral) + "_" + integral.label();
        
    return fstr::lowercase(label);
}

std::string
get_buffer_label(const I4CIntegral& integral,
                 const std::string& prefix)
{
    std::string label = prefix + "_buffer_";
    
    label += std::to_string(integral.order()) + "_";
    
    label += fstr::lowercase(integral.label());

    return label;
}

std::string
prim_compute_func_name(const I4CIntegral& integral)
{
    auto label =  "comp_prim_" + t4c::integral_split_label(integral) + "_" + integral.label();
    
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

std::string
prim_file_name(const I4CIntegral& integral)
{
    return t4c::integral_label(integral) + "PrimRec" + integral.label();
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

} // t4c namespace
