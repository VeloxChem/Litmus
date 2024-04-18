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

#include "t2c_prim_body.hpp"

#include "t2c_utils.hpp"

void
T2CPrimFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                       const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    lines.push_back({1, 0, 2, "const auto ndims = " +
                              t2c::get_buffer_label(integral, "prim") +
                              ".number_of_columns();"});
    
    for (const auto& label : _get_buffers_str(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    const auto components = integral.components<T1CPair, T1CPair>();
    
    if ((integral[0] == 0) || (integral[1] == 0))
    {
        _add_recursion_loop(lines, integral, components, {0, static_cast<int>(components.size())});
    }
    else
    {
        // TODO : implement generic code
    }
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CPrimFuncBodyDriver::_get_buffers_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    for (const auto& tint : t2c::get_integrals(integral))
    {
        auto label = t2c::get_buffer_label(tint, "prim");
        
        vstr.push_back("/// Set up components of auxilary buffer : " + label);
        
        const auto tlabel = _get_tensor_label(tint);
        
        int index = 0;
        
        for (const auto& tcomp : tint.components<T1CPair, T1CPair>())
        {
            const auto line = "auto " + tlabel + "_" + tcomp.label() + " = " + label;
            
            vstr.push_back(line + "[" + std::to_string(index) + "];");
            
            index++;
        }
    }
    
    return vstr;
}

std::string
T2CPrimFuncBodyDriver::_get_tensor_label(const I2CIntegral& integral) const
{
    std::string label;
    
    if (integral.integrand().name() == "1") label = "ts";
    
    return label;
}

void
T2CPrimFuncBodyDriver::_add_recursion_loop(      VCodeLines&         lines,
                                           const I2CIntegral&        integral,
                                           const VT2CIntegrals&      components,
                                           const std::array<int, 2>& rec_range) const
{
    const bool full_rec = (rec_range[1] - rec_range[0]) == (static_cast<int>(components.size()));
    
    // set up recursion buffer
    
    auto label = t2c::get_buffer_label(integral, "prim");
    
    if (full_rec)
    {
        lines.push_back({1, 0, 2, "/// Set up components of targeted buffer : " + label});
    }
    else
    {
        lines.push_back({1, 0, 2, "/// Set up " +
                                   std::to_string(rec_range[0]) + "-" + std::to_string(rec_range[1]) +
                                   " components of targeted buffer : " + label});
    }
    
    const auto tlabel = _get_tensor_label(integral);
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        const auto line = "auto " + tlabel + "_" + components[i].label() + " = " + label;
        
        lines.push_back({1, 0, 2, line + "[" + std::to_string(i) + "];"});
    }
    
    // set up recursion loop
    
    lines.push_back({1, 0, 1, "for (std::size_t i = 0; i < ndims; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({1, 0, 1, "}"});
}
