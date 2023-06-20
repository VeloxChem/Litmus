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
T2CDeclDriver::write_func_decl(      std::ofstream& fstream,
                               const I2CIntegral&   integral,
                               const bool           diagonal,
                               const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_matrix_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_gto_blocks_str(integral, diagonal))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_indexes_str(integral,diagonal, terminus))
    {
        if  (label.find(";") == std::string::npos)
        {
            lines.push_back({0, 0, 1, label});
        }
        else
        {
            lines.push_back({0, 0, 2, label});
        }
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CDeclDriver::write_prim_func_decl(      std::ofstream& fstream,
                                    const I2CIntegral&   integral,
                                    const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_prim_buffer_str(integral, terminus))
    {
        if  (label.find(";") == std::string::npos)
        {
            lines.push_back({0, 0, 1, label});
        }
        else
        {
            lines.push_back({0, 0, 2, label});
        }
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CDeclDriver::write_prim_func_decl(      std::ofstream&   fstream,
                                    const TensorComponent& component,
                                    const I2CIntegral&     integral,
                                    const bool             bra_first,
                                    const bool             terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_prim_buffer_str(component, integral, bra_first, terminus))
    {
        if  (label.find(";") == std::string::npos)
        {
            lines.push_back({0, 0, 1, label});
        }
        else
        {
            lines.push_back({0, 0, 2, label});
        }
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CDeclDriver::write_prim_func_decl(      std::ofstream&   fstream,
                                    const TensorComponent& bra_component,
                                    const TensorComponent& ket_component,
                                    const I2CIntegral&     integral,
                                    const bool             terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_prim_buffer_str(bra_component, ket_component, integral, terminus))
    {
        if  (label.find(";") == std::string::npos)
        {
            lines.push_back({0, 0, 1, label});
        }
        else
        {
            lines.push_back({0, 0, 2, label});
        }
    }
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CDeclDriver::_get_matrix_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    const auto [nsize, name] = t2c::compute_func_name(integral);
    
    const auto labels = t2c::integrand_components(integral.integrand(), "matrix");
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        if (i == 0)
        {
            vstr.push_back(name + "(" + std::string(6, ' ') + "CSubMatrix* " + labels[i]);
        }
        else
        {
            vstr.push_back(std::string(nsize + 6, ' ') + "CSubMatrix* " + labels[i]);
        }
    }
    
    return vstr;
}

std::vector<std::string>
T2CDeclDriver::_get_gto_blocks_str(const I2CIntegral& integral,
                                   const bool         diagonal) const
{
    std::vector<std::string> vstr;
    
    const auto [nsize, name] = t2c::compute_func_name(integral);
    
    if (diagonal)
    {
        vstr.push_back(std::string(nsize, ' ') + "const CGtoBlock&  gto_block,");
    }
    else
    {
        vstr.push_back(std::string(nsize, ' ') + "const CGtoBlock&  bra_gto_block,");
        
        vstr.push_back(std::string(nsize, ' ') + "const CGtoBlock&  ket_gto_block,");
    }
    
    if (integral[0] != integral[1])
    {
        vstr.push_back(std::string(nsize, ' ') + "const bool        ang_order,");
    }
    
    return vstr;
}

std::vector<std::string>
T2CDeclDriver::_get_indexes_str(const I2CIntegral& integral,
                                const bool         diagonal,
                                const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto [nsize, name] = t2c::compute_func_name(integral);
    
    vstr.push_back(std::string(nsize, ' ') + "const int64_t     bra_first,");
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    if ((!diagonal) && (integral[0] == integral[1]))
    {
        vstr.push_back(std::string(nsize, ' ') + "const int64_t     bra_last,");
        
        vstr.push_back(std::string(nsize, ' ') + "const mat_t       mat_type) -> void" +  tsymbol);
    }
    else
    {
        vstr.push_back(std::string(nsize, ' ') + "const int64_t     bra_last) -> void" + tsymbol);
    }
    
    return vstr;
}

std::vector<std::string>
T2CDeclDriver::_get_prim_buffer_str(const I2CIntegral& integral,
                                    const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto [nsize, name] = t2c::prim_compute_func_name(integral);
    
    std::vector<std::string> labels({"buffer", });
    
    if (integral[0] > 0) labels = t2c::tensor_components(Tensor(integral[0]), "buffer");
    
    if (integral[1] > 0) labels = t2c::tensor_components(Tensor(integral[1]), "buffer");
    
    vstr.push_back(name + "(      TDoubleArray& " + labels[0] + ",");
   
    for (size_t i = 1; i < labels.size(); i++)
    {
        vstr.push_back(std::string(nsize + 6, ' ') + "TDoubleArray& " + labels[i] + ",");
    }
    
    _add_prim_variables(vstr, std::string(nsize, ' '), terminus); 
    
    return vstr;
}

std::vector<std::string>
T2CDeclDriver::_get_prim_buffer_str(const TensorComponent& component,
                                    const I2CIntegral&     integral,
                                    const bool             bra_first,
                                    const bool             terminus) const
{
    std::vector<std::string> vstr;
    
    const auto [nsize, name] = t2c::prim_compute_func_name(component, integral, bra_first);
    
    const auto labels = (bra_first) ? t2c::tensor_components(integral[1], "buffer")
                                    : t2c::tensor_components(integral[0], "buffer");
    
    vstr.push_back(name + "(      TDoubleArray& " + labels[0] + ",");
   
    for (size_t i = 1; i < labels.size(); i++)
    {
        vstr.push_back(std::string(nsize + 6, ' ') + "TDoubleArray& " + labels[i] + ",");
    }
    
    _add_prim_variables(vstr, std::string(nsize, ' '), terminus);
    
    return vstr;
}


std::vector<std::string>
T2CDeclDriver::_get_prim_buffer_str(const TensorComponent& bra_component,
                                    const TensorComponent& ket_component,
                                    const I2CIntegral&     integral,
                                    const bool             terminus) const
{
    std::vector<std::string> vstr;
    
    const auto [nsize, name] = t2c::prim_compute_func_name(bra_component, ket_component, integral);
    
    const auto labels = t2c::integrand_components(integral.integrand(), "buffer");
    
    vstr.push_back(name + "(      TDoubleArray& " + labels[0] + ",");
   
    for (size_t i = 1; i < labels.size(); i++)
    {
        vstr.push_back(std::string(nsize + 6, ' ') + "TDoubleArray& " + labels[i] + ",");
    }
    
    _add_prim_variables(vstr, std::string(nsize, ' '), terminus);
    
    return vstr;
}

void
T2CDeclDriver::_add_prim_variables(      std::vector<std::string>& vstrings,
                                   const std::string&              spacer,
                                   const bool                      terminus) const
{
    vstrings.push_back(spacer + "const double        bra_exp,");
    
    vstrings.push_back(spacer + "const double        bra_norm,");
        
    vstrings.push_back(spacer + "const TPoint3D&     bra_coord,");
        
    vstrings.push_back(spacer + "const TDoubleArray& ket_exps,");
        
    vstrings.push_back(spacer + "const TDoubleArray& ket_norms,");
        
    vstrings.push_back(spacer + "const TDoubleArray& ket_coords_x,");
        
    vstrings.push_back(spacer + "const TDoubleArray& ket_coords_y,");
        
    vstrings.push_back(spacer + "const TDoubleArray& ket_coords_z,");
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    vstrings.push_back(spacer + "const int64_t       ket_dim) -> void" + tsymbol);
}
