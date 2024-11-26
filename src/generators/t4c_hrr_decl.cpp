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
    
    auto label = t4c::get_hrr_index(integral, true);
    
    vstr.push_back(name + "CSimdArray<double>& cbuffer," );
    
    vstr.push_back(spacer + "const size_t " + label + "," );
    
    if (integral[2] == 1)
    {
        vstr.push_back(spacer + "CSimdArray<double>& pbuffer," );
    }
    
    for (const auto& tint : t4c::get_ket_hrr_integrals(integral))
    {
        auto label = t4c::get_hrr_index(tint, true);
        
        vstr.push_back(spacer + "const size_t " + label + "," );
    }
    
    return vstr;
}

std::vector<std::string>
T4CHrrDeclDriver::_get_ket_coordinates_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::ket_hrr_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
   
    vstr.push_back(spacer + "const CSimdArray<double>& factors,");
        
    vstr.push_back(spacer + "const size_t idx_cd,");
        
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

void
T4CHrrDeclDriver::write_bra_func_decl(      std::ofstream& fstream,
                                      const I4CIntegral&   integral,
                                      const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_bra_buffers_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_bra_coordinates_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_bra_recursion_variables_str(integral, terminus))
    {
        lines.push_back({0, 0, 1, label});
    }
        
    ost::write_code_lines(fstream, lines);
}

void
T4CHrrDeclDriver::write_bra_geom_func_decl(      std::ofstream& fstream,
                                           const I4CIntegral&   integral,
                                           const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_bra_geom_buffers_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_bra_geom_coordinates_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_bra_geom_recursion_variables_str(integral, terminus))
    {
        lines.push_back({0, 0, 1, label});
    }
        
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CHrrDeclDriver::_get_bra_buffers_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::bra_hrr_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(name + "CSimdArray<double>& cbuffer," );
    
    auto label = t4c::get_hrr_index(integral, false);
    
    vstr.push_back(spacer + "const size_t " + label + "," );
   
    for (const auto& tint : t4c::get_bra_hrr_integrals(integral))
    {
        label = t4c::get_hrr_index(tint, false);
            
        vstr.push_back(spacer + "const size_t " + label + "," );
    }
   
    return vstr;
}

std::vector<std::string>
T4CHrrDeclDriver::_get_bra_geom_buffers_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::bra_geom_hrr_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(name + "CSimdArray<double>& cbuffer," );
    
    auto label = t4c::get_hrr_index(integral, false);
    
    vstr.push_back(spacer + "const size_t " + label + "," );
    
    if (integral[0] == 0)
    {
        for (const auto& tint : t4c::get_aux_geom_hrr_integrals(integral))
        {
            label = t4c::get_hrr_index(tint, false);
            
            vstr.push_back(spacer + "const size_t " + label + "," );
        }
    }
    else
    {
        for (const auto& tint : t4c::get_bra_geom_hrr_integrals(integral))
        {
            label = t4c::get_hrr_index(tint, false);
            
            vstr.push_back(spacer + "const size_t " + label + "," );
        }
    }
    
    return vstr;
}

std::vector<std::string>
T4CHrrDeclDriver::_get_bra_coordinates_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::ket_hrr_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
   
    vstr.push_back(spacer + "const TPoint<double>& r_ab,");
        
    return vstr;
}

std::vector<std::string>
T4CHrrDeclDriver::_get_bra_geom_coordinates_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::bra_geom_hrr_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    const bool no_rab = (integral.prefixes_order() == std::vector<int>({0, 1, 0, 0}))
    
                      && (integral[0] == 0);
    
    if (!no_rab)
    {
        vstr.push_back(spacer + "const TPoint<double>& r_ab,");
    }
   
    return vstr;
}

std::vector<std::string>
T4CHrrDeclDriver::_get_bra_recursion_variables_str(const I4CIntegral& integral,
                                                   const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t4c::ket_hrr_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(spacer + "const int c_angmom,");
        
    vstr.push_back(spacer + "const int d_angmom) -> void" + tsymbol);
 
    return vstr;
}

std::vector<std::string>
T4CHrrDeclDriver::_get_bra_geom_recursion_variables_str(const I4CIntegral& integral,
                                                        const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t4c::bra_geom_hrr_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(spacer + "const int c_angmom,");
        
    vstr.push_back(spacer + "const int d_angmom) -> void" + tsymbol);
 
    return vstr;
}
