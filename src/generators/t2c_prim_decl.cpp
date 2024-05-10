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
    
    vstr.push_back(name + "CSimdArray<double>& " + label + "," );
    
    for (const auto& tint : t2c::get_integrals(integral))
    {
        auto label = t2c::get_buffer_label(tint, "prim");
        
        vstr.push_back(spacer + "const CSimdArray<double>& " + label + "," );
    }
    
    return vstr;
}

std::vector<std::string>
T2CPrimDeclDriver::_get_coordinates_str(const I2CIntegral& integral,
                                        const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    auto name = t2c::prim_compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
   
    if (integral[0] > 0)
    {
        vstr.push_back(spacer + "const double* pa_x,");
        
        vstr.push_back(spacer + "const double* pa_y,");
      
        if ((integral[0] == 1)                   &&
            (integral[1] == 0)                   &&
            (integral.integrand().name() != "T") &&
            (integral.integrand().name() != "A"))
        {
            vstr.push_back(spacer + "const double* pa_z) -> void" + tsymbol);
        }
        else
        {
            vstr.push_back(spacer + "const double* pa_z,");
        }
    }
   
    if ((integral[0] == 0) && (integral[1] > 0))
    {
        vstr.push_back(spacer + "const double* pb_x,");
        
        vstr.push_back(spacer + "const double* pb_y,");
      
        if ((integral[0] == 0)                   &&
            (integral[1] == 1)                   &&
            (integral.integrand().name() != "T") &&
            (integral.integrand().name() != "A"))
        {
            vstr.push_back(spacer + "const double* pb_z) -> void" + tsymbol);
        }
        else
        {
            vstr.push_back(spacer + "const double* pb_z,");
        }
    }
    
    if ((integral.integrand().name() == "A") && ((integral[0] + integral[1]) != 0))
    {
        vstr.push_back(spacer + "const double* pc_x,");
        
        vstr.push_back(spacer + "const double* pc_y,");
        
        if ((integral[0] + integral[1]) > 1)
        {
            vstr.push_back(spacer + "const double* pc_z,");
        }
        else
        {
            vstr.push_back(spacer + "const double* pc_z) -> void" + tsymbol);
        }
    }
    
    if ((integral[0] + integral[1]) == 0)
    {
        if (integral.integrand().name() == "A")
        {
            vstr.push_back(spacer + "const double* bf_values,");
        }
        else
        {
            vstr.push_back(spacer + "const double* ab_x,");
            
            vstr.push_back(spacer + "const double* ab_y,");
        
            vstr.push_back(spacer + "const double* ab_z,");
        }
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
        vstr.push_back(spacer + "const double a_exp,");
        
        vstr.push_back(spacer + "const double* b_exps) -> void" + tsymbol);
        
    }
    else
    {
        if (((integral[0] + integral[1]) != 1) || (integral.integrand().name() == "T"))
        {
            vstr.push_back(spacer + "const double a_exp,");
            
            if (((integral[0] + integral[1]) == 0) && (integral.integrand().name() == "1"))
            {
                vstr.push_back(spacer + "const double* b_exps,");
            }
            else
            {
                vstr.push_back(spacer + "const double* b_exps) -> void" + tsymbol);
            }
        }
       
        if ((integral[0] + integral[1]) == 0)
        {
            if (integral.integrand().name() == "1")
            {
                vstr.push_back(spacer + "const double a_norm,");
                
                vstr.push_back(spacer + "const double* b_norms) -> void" + tsymbol);
            }
        }
    }
    

    return vstr;
}
