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

#include "t2c_geom_decl.hpp"

#include "file_stream.hpp"
#include "t2c_utils.hpp"

void
T2CGeomDeclDriver::write_func_decl(      std::ofstream& fstream,
                                   const SI2CIntegrals& geom_integrals,
                                   const I2CIntegral&   integral,
                                   const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_buffers_str(geom_integrals, integral))
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
T2CGeomDeclDriver::_get_buffers_str(const SI2CIntegrals& geom_integrals,
                                    const I2CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    auto name = t2c::prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(name + "CSimdArray<double>& pbuffer," );
    
    auto label = t2c::get_index_label(integral);
    
    vstr.push_back(spacer + "const size_t " + label + "," );
    
    for (const auto& tint : geom_integrals)
    {
        label = t2c::get_index_label(tint);
        
        vstr.push_back(spacer + "const size_t " + label + "," );
    }
    
    vstr.push_back(spacer + "const size_t op_comps," );
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        if (prefixes[1].shape().order() == 0)
        {
            vstr.push_back(spacer + "const size_t ket_comps," );
        }
    }
    
    return vstr;
}

std::vector<std::string>
T2CGeomDeclDriver::_get_recursion_variables_str(const I2CIntegral& integral,
                                                const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t2c::prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    const auto prefixes = integral.prefixes();
    
    if (prefixes[0].shape().order() > 0)
    {
        if (prefixes[1].shape().order() > 0)
        {
            vstr.push_back(spacer + "const CSimdArray<double>& factors,");
        }
        
        vstr.push_back(spacer + "const double a_exp) -> void" + tsymbol);
    }
    
    return vstr;
}
