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
// limitations under the License. by Zilvinas Rinkevicius on 2025-01-22.
//

#include "t3c_prim_decl.hpp"

#include "file_stream.hpp"
#include "t3c_utils.hpp"

void
T3CPrimDeclDriver::write_func_decl(      std::ofstream&         fstream,
                                   const I3CIntegral&           integral,
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
T3CPrimDeclDriver::_get_buffers_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t3c::prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(name + "CSimdArray<double>& pbuffer," );
    
    auto label = t3c::get_index_label(integral);
    
    vstr.push_back(spacer + "const size_t " + label + "," );
    
    for (const auto& tint : t3c::get_vrr_integrals(integral))
    {
        label = t3c::get_index_label(tint);
        
        vstr.push_back(spacer + "size_t " + label + "," );
    }
    
    return vstr;
}

std::vector<std::string>
T3CPrimDeclDriver::_get_coordinates_str(const I3CIntegral& integral,
                                        const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t3c::prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(spacer + "CSimdArray<double>& factors,");
   
    if (integral[0] > 0)
    {
        vstr.push_back(spacer + "const size_t idx_wa,");
    }
   
    if ((integral[0] == 0) && (integral[2] > 0))
    {
        vstr.push_back(spacer + "const size_t idx_qd,");
        
        if ((integral[0] == 0) && (integral[2] == 1))
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
T3CPrimDeclDriver::_get_recursion_variables_str(const I3CIntegral& integral,
                                                const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t3c::prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');

    if ((integral[0] + integral[2]) > 1)
    {
        vstr.push_back(spacer + "const double a_exp) -> void" + tsymbol);
    }
   
    return vstr;
}
