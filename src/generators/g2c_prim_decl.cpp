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

#include "g2c_prim_decl.hpp"

#include "file_stream.hpp"
#include "t2c_utils.hpp"

void
G2CPrimDeclDriver::write_func_decl(      std::ofstream&         fstream,
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
G2CPrimDeclDriver::_get_buffers_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t2c::grid_prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(name + "CSubMatrix&  buffer," );
    
    auto label = t2c::get_index_label(integral);
    
    vstr.push_back(spacer + "const size_t " + label + "," );
    
    for (const auto& tint : t2c::get_integrals(integral))
    {
        label = t2c::get_index_label(tint);
        
        vstr.push_back(spacer + "const size_t " + label + "," );
    }
    
    return vstr;
}

std::vector<std::string>
G2CPrimDeclDriver::_get_coordinates_str(const I2CIntegral& integral,
                                        const bool         terminus) const
{
    const auto tsymbol = (terminus) ? ";" : "";
    
    std::vector<std::string> vstr;
    
    auto name = t2c::grid_prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
   
    if (integral[0] > 0)
    {
        vstr.push_back(spacer + "const double rpa_x,");
        
        vstr.push_back(spacer + "const double rpa_y,");
        
        if (!_need_exponents(integral))
        {
            vstr.push_back(spacer + "const double rpa_z) -> void" + std::string(tsymbol));
        }
        else
        {
            vstr.push_back(spacer + "const double rpa_z,");
        }
    }
   
    if ((integral[0] == 0) && (integral[1] > 0))
    {
        vstr.push_back(spacer + "const double rpb_x,");
        
        vstr.push_back(spacer + "const double rpb_y,");
        
        if (!_need_exponents(integral))
        {
            vstr.push_back(spacer + "const double rpb_z) -> void" + std::string(tsymbol));
        }
        else
        {
            vstr.push_back(spacer + "const double rpb_z,");
        }
    }
        
    return vstr;
}

std::vector<std::string>
G2CPrimDeclDriver::_get_recursion_variables_str(const I2CIntegral& integral,
                                                const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t2c::grid_prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    if (_need_exponents(integral))
    {
        vstr.push_back(spacer + "const double factor) -> void" + std::string(tsymbol));
    }
    
    return vstr;
}

bool
G2CPrimDeclDriver::_need_exponents(const I2CIntegral& integral) const
{
    return (integral[0] + integral[1]) > 1;
}

