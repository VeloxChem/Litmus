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

#include "t4c_prim_decl.hpp"

#include "file_stream.hpp"
#include "t4c_utils.hpp"

void
T4CPrimDeclDriver::write_func_decl(      std::ofstream&         fstream,
                                   const I4CIntegral&           integral,
                                   const bool                   terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_buffers_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_coordinates_str(integral, terminus))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_recursion_variables_str(integral, terminus))
    {
        lines.push_back({0, 0, 1, label});
    }
        
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CPrimDeclDriver::_get_buffers_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(name + "CSimdArray<double>& pbufer," );
    
    auto label = t4c::get_index_label(integral);
    
    vstr.push_back(spacer + "const size_t " + label + "," );
    
    for (const auto& tint : t4c::get_vrr_integrals(integral))
    {
        label = t4c::get_index_label(tint);
        
        vstr.push_back(spacer + "size_t " + label + "," );
    }
    
    return vstr;
}

std::vector<std::string>
T4CPrimDeclDriver::_get_coordinates_str(const I4CIntegral& integral,
                                        const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t4c::prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(spacer + "CSimdArray<double>& factors,");
   
    if (integral[1] > 0)
    {
        vstr.push_back(spacer + "const size_t idx_wp,");
        
        if ((integral[1] == 1) && (integral[3] == 0))
        {
            vstr.push_back(spacer + "const TPoint<double>& r_pb) -> void" + tsymbol);
        }
        else
        {
            vstr.push_back(spacer + "const TPoint<double>& r_pb,");
        }
    }
   
    if ((integral[1] == 0) && (integral[3] > 0))
    {
        vstr.push_back(spacer + "const size_t idx_qd,");
        
        if ((integral[1] == 0) && (integral[3] == 1))
        {
            vstr.push_back(spacer + "const size_t idx_wq) -> void" + tsymbol);
        }
        else
        {
            vstr.push_back(spacer + "const size_t idx_wq,");
        }
    }
    
    return vstr;
}

std::vector<std::string>
T4CPrimDeclDriver::_get_recursion_variables_str(const I4CIntegral& integral,
                                                const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t4c::prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    
//    if ((integral[1] > 1) && (integral[3] == 0))
//    {
//        vstr.push_back(spacer + "const double a_exp,");
//
//        vstr.push_back(spacer + "const double b_exp) -> void" + tsymbol);
//    }
//
//    if ((integral[1] == 0) && (integral[3] > 1))
//    {
//        vstr.push_back(spacer + "const double* c_exps,");
//
//        vstr.push_back(spacer + "const double* d_exps) -> void" + tsymbol);
//    }
    
    if ((integral[1] + integral[3]) > 1)
    {
        vstr.push_back(spacer + "const double a_exp,");
        
        vstr.push_back(spacer + "const double b_exp) -> void" + tsymbol);
    }
   
    if ((integral[1] + integral[3]) == 0)
    {
        vstr.push_back(spacer + "const size_t idx_ovl,");
        
        vstr.push_back(spacer + "const CSimdArray<double>& bf_data,");
        
        vstr.push_back(spacer + "const size_t idx_bvals) -> void" + tsymbol);
    }
    
    return vstr;
}
