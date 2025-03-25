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

#include "t2c_prim_decl.hpp"

#include "file_stream.hpp"
#include "t2c_utils.hpp"

void
T2CPrimDeclDriver::write_func_decl(      std::ofstream&         fstream,
                                   const I2CIntegral&           integral,
                                   const bool                   terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_buffers_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    if (integral.is_simple())
    {
        for (const auto& label : _get_coordinates_str(integral, terminus))
        {
            lines.push_back({0, 0, 1, label});
        }
    }
    
    for (const auto& label : _get_recursion_variables_str(integral, terminus))
    {
        lines.push_back({0, 0, 1, label});
    }
        
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CPrimDeclDriver::_get_buffers_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t2c::prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    auto label = t2c::get_buffer_label(integral, "prim");
    
    vstr.push_back(name + "CSimdArray<double>& pbuffer, " );
    
    label = t2c::get_index_label(integral);
    
    vstr.push_back(spacer + "const size_t " + label + "," );
    
    for (const auto& tint : t2c::get_integrals(integral))
    {
        label = t2c::get_index_label(tint);
        
        vstr.push_back(spacer + "const size_t " + label + "," );
    }
    
    return vstr;
}

std::vector<std::string>
T2CPrimDeclDriver::_get_coordinates_str(const I2CIntegral& integral,
                                        const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    auto name = t2c::prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(spacer + "const CSimdArray<double>& factors,");
   
    if ((integral[0] > 0) && (integral.integrand().name() != "GX(r)"))
    {
        if (integral.integrand().name() == "G(r)")
        {
            vstr.push_back(spacer + "const size_t idx_rga,");
        }
        else
        {
            vstr.push_back(spacer + "const size_t idx_rpa,");
        }
    }
   
    if ((integral[0] == 0) && (integral[1] > 0) && (integral.integrand().name() != "GX(r)"))
    {
        if (integral.integrand().name() == "G(r)")
        {
            vstr.push_back(spacer + "const size_t idx_rgb,");
        }
        else
        {
            vstr.push_back(spacer + "const size_t idx_rpb,");
        }
    }
    
    if (integral.integrand().name() == "GX(r)")
    {
        vstr.push_back(spacer + "const size_t idx_rgc,");
    }
    
    if (_need_distances_pc(integral))
    {
        vstr.push_back(spacer + "const size_t idx_rpc,");
    }
        
    if (!_need_exponents(integral))
    {
        const auto idx = vstr.size() - 1;
        
        vstr[idx].pop_back();
        
        const auto tsymbol = (terminus) ? ";" : "";
        
        vstr[idx] += std::string(") -> void") + std::string(tsymbol);
    }
    
    return vstr;
}

std::vector<std::string>
T2CPrimDeclDriver::_get_recursion_variables_str(const I2CIntegral& integral,
                                                const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t2c::prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    if (!integral.is_simple())
    {
        vstr.push_back(spacer + "const double a_exp) -> void" + tsymbol);
    }
    else
    {
        if (((integral[0] + integral[1]) != 1) || (integral.integrand().name() == "T") || (integral.integrand().name() == "GX(r)"))
        {
            if ((integral.integrand().name() == "G(r)") || (integral.integrand().name() == "GX(r)"))
            {
                vstr.push_back(spacer + "const double a_exp,");
                
                vstr.push_back(spacer + "const double c_exp) -> void" + tsymbol);
            }
            else
            {
                vstr.push_back(spacer + "const double a_exp) -> void" + tsymbol);
            }
        }

        if ( ((integral[0] + integral[1]) == 1)  && (integral.integrand().name() == "r"))
        {
            vstr.push_back(spacer + "const double a_exp) -> void" + tsymbol);
        }
    }
    
    return vstr;
}

bool
T2CPrimDeclDriver::_need_exponents(const I2CIntegral& integral) const
{
    if (integral.integrand().name() == "T") return true;
    
    if (integral.integrand().name() == "GX(r)") return true;
    
    if (integral.integrand().name() == "r") return true;
    
    return (integral[0] + integral[1]) > 1;
}

bool
T2CPrimDeclDriver::_need_distances_pc(const I2CIntegral& integral) const
{
    if (integral.integrand().name() == "A") return true;
    
    if (integral.integrand().name() == "AG") return true;
    
    return false;
}
