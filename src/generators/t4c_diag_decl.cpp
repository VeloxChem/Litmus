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

#include "t4c_diag_decl.hpp"

#include "file_stream.hpp"
#include "t4c_utils.hpp"

void
T4CDiagDeclDriver::write_func_decl(      std::ofstream& fstream,
                                   const I4CIntegral&   integral,
                                   const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_vars_str(integral, terminus))
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
T4CDiagDeclDriver::write_prim_func_decl(      std::ofstream& fstream,
                                        const T4CIntegral&   component,
                                        const I4CIntegral&   integral,
                                        const bool           diagonal,
                                        const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_prim_vars_str(component, integral, diagonal, terminus))
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
T4CDiagDeclDriver::_get_vars_str(const I4CIntegral& integral,
                                 const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto [nsize, name] = t4c::diag_compute_func_name(integral);
    
    auto label = name + "(const CGtoPairBlock& gto_pair_block) -> std::vector<double>";
    
    if (terminus) label += ";"; 
        
    vstr.push_back(label);
    
    return vstr;
}

std::vector<std::string>
T4CDiagDeclDriver::_get_prim_vars_str(const T4CIntegral& component,
                                      const I4CIntegral& integral,
                                      const bool         diagonal,
                                      const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto [nsize, name] = t4c::prim_diag_compute_func_name(component, integral);
    
    vstr.push_back(name + "(TDoubleArray& buffer,");
   
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_a_x,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_a_y,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_a_z,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_b_x,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_b_y,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_b_z,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& bra_exps_a,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& bra_exps_b,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& bra_norms,");
    
    if (!diagonal)
    {
        vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& ket_exps_c,");
        
        vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& ket_exps_d,");
        
        vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& ket_norms,");
    }
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    vstr.push_back(std::string(nsize, ' ') + "const int64_t       ndim) -> void" + tsymbol);
    
    return vstr;
}