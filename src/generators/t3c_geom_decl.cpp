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

#include "t3c_geom_decl.hpp"

#include "file_stream.hpp"
#include "t3c_utils.hpp"

void
T3CGeomDeclDriver::write_func_decl(      std::ofstream& fstream,
                                   const I3CIntegral&   integral,
                                   const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "template <class T>"});
    
    lines.push_back({0, 0, 1, "inline auto"});
    
    for (const auto& label : _get_matrices_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_gto_pair_blocks_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_indices_str(integral, terminus))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T3CGeomDeclDriver::_get_matrices_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t3c::compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(name + "T& distributor,");
    
    return vstr;
}

std::vector<std::string>
T3CGeomDeclDriver::_get_gto_pair_blocks_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t3c::compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(spacer + "const CGtoBlock& bra_gto_block,");
        
    vstr.push_back(spacer + "const CGtoPairBlock& ket_gto_pair_block,");
    
    return vstr;
}

std::vector<std::string>
T3CGeomDeclDriver::_get_indices_str(const I3CIntegral& integral,
                                    const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    auto name = t3c::compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    vstr.push_back(spacer + "const std::pair<size_t, size_t>& bra_range) -> void" + tsymbol);
    
    return vstr;
}

void
T3CGeomDeclDriver::write_bra_geom_func_decl(      std::ofstream& fstream,
                                            const I3CIntegral&   integral,
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
T3CGeomDeclDriver::_get_bra_geom_buffers_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    const auto gorders = integral.prefixes_order();
    
    auto name = t3c::bra_geom_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(name + "CSimdArray<double>& cbuffer," );
    
    std::string label;
    
    if (gorders == std::vector<int>({1, 0, 0}))
    {
        label = t3c::get_full_hrr_index(integral, false);
    }
    else
    {
        label = t3c::get_hrr_index(integral);
    }
  
    vstr.push_back(spacer + "const size_t " + label + "," );
    
    for (const auto& tint : t3c::get_bra_geom_integrals(integral))
    {
        if (gorders[0] > 0)
        {
            label = t3c::get_full_hrr_index(tint, false);
        }
        else
        {
            label = t3c::get_hrr_index(tint);
        }
            
        vstr.push_back(spacer + "const size_t " + label + "," );
    }
    
    return vstr;
}

std::vector<std::string>
T3CGeomDeclDriver::_get_bra_geom_coordinates_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
//    auto name = t3c::bra_geom_compute_func_name(integral) + "(";
//    
//    const auto spacer = std::string(name.size(), ' ');
//    
//    const bool no_rab = (integral.prefixes_order() == std::vector<int>({0, 1, 0, 0}))
//    
//                      && (integral[0] == 0);
//    
//    if (!no_rab)
//    {
//        vstr.push_back(spacer + "const TPoint<double>& r_ab,");
//    }
   
    return vstr;
}

std::vector<std::string>
T3CGeomDeclDriver::_get_bra_geom_recursion_variables_str(const I3CIntegral& integral,
                                                         const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t3c::bra_geom_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(spacer + "const int c_angmom,");
        
    vstr.push_back(spacer + "const int d_angmom) -> void" + tsymbol);
 
    return vstr;
}

void
T3CGeomDeclDriver::write_ket_geom_func_decl(      std::ofstream& fstream,
                                            const I3CIntegral&   integral,
                                            const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_ket_geom_buffers_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_ket_geom_coordinates_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_ket_geom_recursion_variables_str(integral, terminus))
    {
        lines.push_back({0, 0, 1, label});
    }
        
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T3CGeomDeclDriver::_get_ket_geom_buffers_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    const auto gorders = integral.prefixes_order();
    
    auto name = t3c::ket_geom_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(name + "CSimdArray<double>& cbuffer," );
    
    std::string label = t3c::get_hrr_index(integral);
    
    vstr.push_back(spacer + "const size_t " + label + "," );
    
    for (const auto& tint : t3c::get_geom_hrr_integrals(integral))
    {
        if (gorders[0] > 0)
        {
            label = t3c::get_full_hrr_index(tint, false);
        }
        else
        {
            label = t3c::get_hrr_index(tint);
        }
            
        vstr.push_back(spacer + "const size_t " + label + "," );
    }
    
    return vstr;
}

std::vector<std::string>
T3CGeomDeclDriver::_get_ket_geom_coordinates_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;

    auto name = t3c::ket_geom_compute_func_name(integral) + "(";
        
    const auto spacer = std::string(name.size(), ' ');
       
    vstr.push_back(spacer + "const CSimdArray<double>& factors,");
            
    vstr.push_back(spacer + "const size_t idx_cd,");
    
    return vstr;
}

std::vector<std::string>
T3CGeomDeclDriver::_get_ket_geom_recursion_variables_str(const I3CIntegral& integral,
                                                         const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t3c::ket_geom_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
        
    vstr.push_back(spacer + "const int a_angmom) -> void" + tsymbol);
 
    return vstr;
}

