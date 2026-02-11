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

#include "t2c_prim_docs.hpp"

#include "file_stream.hpp"
#include "t2c_utils.hpp"

void
T2CPrimDocuDriver::write_doc_str(      std::ofstream& fstream,
                                 const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_compute_str(integral)});
    
    // TODO: Add special variables here
    
    for (const auto& label : _get_buffers_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_coordinates_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_recursion_variables_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CPrimDocuDriver::write_doc_str(      std::ofstream& fstream,
                                 const M2Integral&    integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_compute_str(integral.second)});
    
    // TODO: Add special variables here
    
    for (const auto& label : _get_buffers_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
//    for (const auto& label : _get_coordinates_str(integral))
//    {
//        lines.push_back({0, 0, 1, label});
//    }
//    
//    for (const auto& label : _get_recursion_variables_str(integral))
//    {
//        lines.push_back({0, 0, 1, label});
//    }
    
    ost::write_code_lines(fstream, lines);
}


std::string
T2CPrimDocuDriver::_get_compute_str(const I2CIntegral& integral) const
{
    const auto bra = Tensor(integral[0]);
        
    const auto ket = Tensor(integral[1]);
    
    const auto [bra_prefix, ket_prefix] = t2c::prefixes_label(integral);
    
    const auto integrand = integral.integrand();
    
    auto label = "/// @brief Computes primitive [" + bra_prefix + bra.label() + "|";
    
    if (integral.integrand().name() != "1")
    {
        label += t2c::integrand_label(integral.integrand()) + "|";
    }
    
    label += ket_prefix + ket.label() + "]";
    
    if (integral.integrand().name() == "U_l")
    {
        label += "_" + Tensor(integral.order()).label();
    }
    
    label += " integrals for set of data buffers.";
    
    return label;
}

std::vector<std::string>
T2CPrimDocuDriver::_get_buffers_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param pbuffer The primitive integrals buffer.");
    
    auto label = t2c::get_index_label(integral);
    
    vstr.push_back("/// @param " + label + " The index of integral in primitive integrals buffer.");
    
    for (const auto& tint : t2c::get_integrals(integral))
    {
        label = t2c::get_index_label(tint);
        
        vstr.push_back("/// @param " + label + " The index of integral in primitive integrals buffer.");
    }

    return vstr;
}

std::vector<std::string>
T2CPrimDocuDriver::_get_buffers_str(const M2Integral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param pbuffer The primitive integrals buffer.");
    
    auto label = t2c::get_index_label(integral);
    
    vstr.push_back("/// @param " + label + " The index of integral in primitive integrals buffer.");
    
    for (const auto& tint : t2c::get_common_integrals(integral))
    {
        label = t2c::get_index_label(tint);
        
        vstr.push_back("/// @param " + label + " The index of integral in primitive integrals buffer.");
    }

    return vstr;
}



std::vector<std::string>
T2CPrimDocuDriver::_get_coordinates_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param factors The primitive factors buffer.");
   
    if (integral.integrand().name() == "U_L") return vstr;
    
    if (
        (integral[0] > 0) &&
        (integral.integrand().name() != "GX(r)")  &&
        (integral.integrand().name() != "GR2(r)") &&
        (integral.integrand().name() != "GR.R2(r)")
        )
    {
        if (integral.integrand().name() == "G(r)")
        {
            vstr.push_back("/// @param idx_rga The vector of distances R(GA) = G - A.");
        }
        else
        {
            vstr.push_back("/// @param idx_rpa The vector of distances R(PA) = P - A.");
        }
        
    }
   
    if (
        (integral[0] == 0) && (integral[1] > 0)   &&
        (integral.integrand().name() != "GX(r)")  &&
        (integral.integrand().name() != "GR2(r)") &&
        (integral.integrand().name() != "GR.R2(r)")
        )
    {
        if (integral.integrand().name() == "G(r)")
        {
            vstr.push_back("/// @param idx_rgb The vector of distances R(GB) = G - B.");
        }
        else
        {
            vstr.push_back("/// @param idx_rpb The vector of distances R(PB) = P - B.");
        }
    }
    
    if (
        (integral.integrand().name() == "GX(r)")  ||
        (integral.integrand().name() == "GR2(r)") ||
        (integral.integrand().name() == "GR.R2(r)")
        )
    {
        vstr.push_back("/// @param idx_rgc The vector of distances R(GC) = G - C.");
    }
    
    if (_need_distances_pc(integral))
    {
        vstr.push_back("/// @param idx_rpc The vector of distances R(PC) = P - C.");
    }
    
    return vstr;
}

std::vector<std::string>
T2CPrimDocuDriver::_get_recursion_variables_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (integral.integrand().name() == "U_L") return vstr;
    
    if (_need_exponents(integral))
    {
        vstr.push_back("/// @param a_exp The primitive basis function exponent on center A.");
        
        if (
            (integral.integrand().name() == "G(r)")   ||
            (integral.integrand().name() == "GX(r)")  ||
            (integral.integrand().name() == "GR2(r)") ||
            (integral.integrand().name() == "GR.R2(r)")
            )
        {
            vstr.push_back("/// @param c_exp The primitive basis function exponent on center C.");
        }
    }
    
    return vstr;
}

bool
T2CPrimDocuDriver::_need_distances_pc(const I2CIntegral& integral) const
{
    if (integral.integrand().name() == "A") return true;
    
    if (integral.integrand().name() == "AG") return true;
    
    return false;
}

bool
T2CPrimDocuDriver::_need_exponents(const I2CIntegral& integral) const
{
    if (integral.integrand().name() == "T") return true;
    
    if (integral.integrand().name() == "GX(r)") return true;
    
    if (integral.integrand().name() == "GR2(r)") return true;
    
    if (integral.integrand().name() == "GR.R2(r)") return true;
    
    if (integral.integrand().name() == "r") return true;
    
    return (integral[0] + integral[1]) > 1;
}
