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

#include "t4c_geom_decl.hpp"

#include "file_stream.hpp"
#include "t4c_utils.hpp"

void
T4CGeomDeclDriver::write_func_decl(      std::ofstream& fstream,
                                   const I4CIntegral&   integral,
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
T4CGeomDeclDriver::_get_matrices_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(name + "T& distributor,");
    
    return vstr;
}

std::vector<std::string>
T4CGeomDeclDriver::_get_gto_pair_blocks_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(spacer + "const CGtoPairBlock& bra_gto_pair_block,");
        
    vstr.push_back(spacer + "const CGtoPairBlock& ket_gto_pair_block,");
    
    return vstr;
}

std::vector<std::string>
T4CGeomDeclDriver::_get_indices_str(const I4CIntegral& integral,
                                const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    vstr.push_back(spacer + "const std::pair<size_t, size_t>& bra_indices,");
        
    vstr.push_back(spacer + "const std::pair<size_t, size_t>& ket_indices) -> void" + tsymbol);
    
    return vstr;
}

std::vector<std::string>
T4CGeomDeclDriver::_get_buffers_str(const SI4CIntegrals& geom_integrals,
                                    const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::geom_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    auto label = t4c::get_geom_buffer_label(integral);
    
    vstr.push_back(name + "CSimdArray<double>& " + label + "," );
    
    for (const auto& tint : geom_integrals)
    {
        auto label = t4c::get_geom_buffer_label(tint);
        
        vstr.push_back(spacer + "const CSimdArray<double>& " + label + "," );
    }
    
    return vstr;
}

std::vector<std::string>
T4CGeomDeclDriver::_get_recursion_variables_str(const I4CIntegral& integral,
                                                const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    int a_order = 0;
    
    int b_order = 0;
    
    int c_order = 0;
    
    int d_order = 0;
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        a_order = prefixes[0].shape().order();
        
        b_order = prefixes[1].shape().order();
        
        c_order = prefixes[2].shape().order();
        
        d_order = prefixes[3].shape().order();
    }
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t4c::geom_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    if (a_order > 0)
    {
        if ((b_order + c_order + d_order) != 0)
        {
            vstr.push_back(spacer + "const double a_exp,");
        }
        else
        {
            vstr.push_back(spacer + "const double a_exp) -> void" + tsymbol);
        }
    }
    
    if (b_order > 0)
    {
        if ((c_order + d_order) != 0)
        {
            vstr.push_back(spacer + "const double b_exp,");
        }
        else
        {
            vstr.push_back(spacer + "const double b_exp) -> void" + tsymbol);
        }
    }
    
    if (c_order > 0)
    {
        if (d_order != 0)
        {
            vstr.push_back(spacer + "const double* c_exps,");
        }
        else
        {
            vstr.push_back(spacer + "const double* c_exps) -> void" + tsymbol);
        }
    }
    
    if (d_order > 0)
    {
        vstr.push_back(spacer + "const double* d_exps) -> void" + tsymbol);
    }
    
    return vstr;
}
