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

#include "t2c_decl.hpp"

#include "file_stream.hpp"
#include "t2c_utils.hpp"

void
T2CDeclDriver::write_func_decl(      std::ofstream&         fstream,
                               const I2CIntegral&           integral,
                               const std::pair<bool, bool>& rec_form,
                               const bool                   diagonal,
                               const bool                   terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "template <class T>"});
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_matrices_str(integral, rec_form))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_special_variables_str(integral, rec_form))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_gto_blocks_str(integral, rec_form, diagonal))
    {
        lines.push_back({0, 0, 1, label});
    }
    
//    for (const auto& label : _get_distributor_variables_str(integral, rec_form, diagonal))
//    {
//        lines.push_back({0, 0, 1, label});
//    }
    
    for (const auto& label : _get_indices_str(integral, rec_form, diagonal, terminus))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CDeclDriver::_get_matrices_str(const I2CIntegral&           integral,
                                 const std::pair<bool, bool>& rec_form) const
{
    std::vector<std::string> vstr;
    
    auto name = t2c::compute_func_name(integral, rec_form) + "(";
    
    vstr.push_back(name + "T* distributor,");
    
    return vstr;
}

// MR: This is the only place in this file where one should expect to have to make changes for new integral cases
std::vector<std::string>
T2CDeclDriver::_get_special_variables_str(const I2CIntegral& integral,
                                          const std::pair<bool, bool>& rec_form) const
{
    std::vector<std::string> vstr;
    
    const auto integrand = integral.integrand();
    
    auto name = t2c::compute_func_name(integral, rec_form) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    if (integrand.name() == "A")
    {
        if (rec_form.first)
        {
            vstr.push_back(spacer + "const std::vector<double>& charges,");
            
            vstr.push_back(spacer + "const std::vector<double>& coords_x,");
            
            vstr.push_back(spacer + "const std::vector<double>& coords_y,");
            
            vstr.push_back(spacer + "const std::vector<double>& coords_z,");
        }
        else
        {
            vstr.push_back(spacer + "const double charge,");
            
            vstr.push_back(spacer + "const double coord_x,");
            
            vstr.push_back(spacer + "const double coord_y,");
            
            vstr.push_back(spacer + "const double coord_z,");
        }
    }
    // MR: Here is additional declaration of variables in fn declaration for dipole integrals
    if (integrand.name() == "r")
    {
        vstr.push_back(spacer + "const double coord_x,");

        vstr.push_back(spacer + "const double coord_y,");

        vstr.push_back(spacer + "const double coord_z,");
    }
    
    return vstr;
}

std::vector<std::string>
T2CDeclDriver::_get_gto_blocks_str(const I2CIntegral&           integral,
                                   const std::pair<bool, bool>& rec_form,
                                   const bool                   diagonal) const
{
    std::vector<std::string> vstr;
    
    auto name = t2c::compute_func_name(integral, rec_form) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    if (diagonal)
    {
        vstr.push_back(spacer + "const CGtoBlock& gto_block,");
    }
    else
    {
        vstr.push_back(spacer + "const CGtoBlock& bra_gto_block,");
        
        vstr.push_back(spacer + "const CGtoBlock& ket_gto_block,");
    }
    
    return vstr;
}

std::vector<std::string>
T2CDeclDriver::_get_distributor_variables_str(const I2CIntegral&           integral,
                                              const std::pair<bool, bool>& rec_form,
                                              const bool                   diagonal) const
{
    std::vector<std::string> vstr;
    
    auto name = t2c::compute_func_name(integral, rec_form) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
 
    if (!diagonal)
    {
        if (integral[0] != integral[1])
        {
            vstr.push_back(spacer + "const bool ang_order,");
        }
        else
        {
            vstr.push_back(spacer + "const mat_t mat_type,");
        }
    }
    
    return vstr;
}

std::vector<std::string>
T2CDeclDriver::_get_indices_str(const I2CIntegral&           integral,
                                const std::pair<bool, bool>& rec_form,
                                const bool                   diagonal,
                                const bool                   terminus) const
{
    std::vector<std::string> vstr;
    
    auto name = t2c::compute_func_name(integral, rec_form) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    if (diagonal)
    {
        vstr.push_back(spacer + "const std::array<int, 2>& gto_range) -> void" + tsymbol);
    }
    else
    {
        vstr.push_back(spacer + "const std::array<int, 2>& bra_range,");
        
        vstr.push_back(spacer + "const std::array<int, 2>& ket_range) -> void" + tsymbol);
    }
   
    return vstr;
}
