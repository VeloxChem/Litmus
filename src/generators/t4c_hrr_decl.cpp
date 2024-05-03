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

#include "t4c_hrr_decl.hpp"

#include "file_stream.hpp"
#include "t4c_utils.hpp"

void
T4CHrrDeclDriver::write_ket_func_decl(      std::ofstream& fstream,
                                      const I4CIntegral&   integral,
                                      const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_ket_buffers_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_ket_coordinates_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_ket_recursion_variables_str(integral, terminus))
    {
        lines.push_back({0, 0, 1, label});
    }
        
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CHrrDeclDriver::_get_ket_buffers_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::ket_hrr_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    auto label = t4c::get_hrr_buffer_label(integral, true);
    
    vstr.push_back(name + "CSimdArray<double>& " + label + "," );
    
    for (const auto& tint : t4c::get_ket_hrr_integrals(integral))
    {
        auto label = t4c::get_hrr_buffer_label(tint, true);
        
        vstr.push_back(spacer + "const CSimdArray<double>& " + label + "," );
    }
    
    return vstr;
}

std::vector<std::string>
T4CHrrDeclDriver::_get_ket_coordinates_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::ket_hrr_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
   
    vstr.push_back(spacer + "const double* cd_x,");
        
    vstr.push_back(spacer + "const double* cd_y,");
        
    vstr.push_back(spacer + "const double* cd_z,");
        
    return vstr;
}

std::vector<std::string>
T4CHrrDeclDriver::_get_ket_recursion_variables_str(const I4CIntegral& integral,
                                                   const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t4c::ket_hrr_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(spacer + "const int a_angmom,");
        
    vstr.push_back(spacer + "const int b_angmom) -> void" + tsymbol);
 
    return vstr;
}
