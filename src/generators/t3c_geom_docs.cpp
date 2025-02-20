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

#include "t3c_geom_docs.hpp"

#include "file_stream.hpp"
#include "t3c_utils.hpp"

void
T3CGeomDocuDriver::write_doc_str(      std::ofstream& fstream,
                                 const I3CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_compute_str(integral)});
    
    for (const auto& label : _get_matrices_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }

    for (const auto& label : _get_gto_pair_blocks_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
        
    for (const auto& label : _get_indices_str())
    {
        lines.push_back({0, 0, 1, label});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T3CGeomDocuDriver::write_bra_geom_doc_str(      std::ofstream& fstream,
                                         const I3CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_bra_geom_compute_str(integral)});
    
   // TODO: Add special variables here
    
    for (const auto& label : _get_bra_geom_buffers_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_bra_coordinates_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_bra_recursion_variables_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    ost::write_code_lines(fstream, lines);
}

std::string
T3CGeomDocuDriver::_get_compute_str(const I3CIntegral& integral) const
{
    const auto bra_one = Tensor(integral[0]);
        
    const auto ket_one = Tensor(integral[1]);
    
    const auto ket_two = Tensor(integral[2]);
    
    const auto integrand = integral.integrand();
    
    std::string label = "/// @brief Computes ";

    label += t3c::prefixes_label(integral);
        
    label += "(" + bra_one.label() +  "|" + t3c::integrand_label(integral.integrand())  + "|";
   
    label += ket_one.label() + ket_two.label() + ")  integral derivatives.";
    
    return label;
}

std::vector<std::string>
T3CGeomDocuDriver::_get_matrices_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param distributor The pointer to Fock matrix/matrices distributor.");
             
    return vstr;
}

std::vector<std::string>
T3CGeomDocuDriver::_get_gto_pair_blocks_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param bra_gto_block The basis functions block on bra side.");
    
    vstr.push_back("/// @param ket_gto_pair_block The basis function pairs block on ket side.");
  
    return vstr;
}

std::vector<std::string>
T3CGeomDocuDriver::_get_indices_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param bra_range The range [bra_first, bra_last) of basis functions on bra side.");
    
    return vstr;
}

std::vector<std::string>
T3CGeomDocuDriver::_get_bra_coordinates_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
   
//    const bool no_rab = (integral.prefixes_order() == std::vector<int>({0, 1, 0, 0}))
//    
//                      && (integral[0] == 0);
//    
//    if (!no_rab)
//    {
//        vstr.push_back("/// @param r_ab The Cartesian distance R(AB) = A - B.");
//    }
            
    return vstr;
}

std::vector<std::string>
T3CGeomDocuDriver::_get_bra_recursion_variables_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param c_angmom The angular momentum on center C.");
    
    vstr.push_back("/// @param d_angmom The angular momentum on center D.");
                       
    return vstr;
}

void
T3CGeomDocuDriver::write_ket_geom_doc_str(      std::ofstream& fstream,
                                         const I3CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_ket_geom_compute_str(integral)});
    
   // TODO: Add special variables here
    
    for (const auto& label : _get_ket_geom_buffers_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_ket_coordinates_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_ket_recursion_variables_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    ost::write_code_lines(fstream, lines);
}

std::string
T3CGeomDocuDriver::_get_ket_geom_compute_str(const I3CIntegral& integral) const
{
    const auto ket_one = Tensor(integral[1]);
    
    const auto ket_two = Tensor(integral[2]);
    
    const auto integrand = integral.integrand();
    
    std::string label = "/// @brief Computes ";

    label += t3c::prefixes_label(integral);
        
    label += "(X|" + t3c::integrand_label(integral.integrand())  + "|";
   
    label += ket_one.label() + ket_two.label() + ")  integral derivatives.";
    
    return label;
}

std::vector<std::string>
T3CGeomDocuDriver::_get_ket_geom_buffers_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto label = t3c::get_hrr_index(integral);
    
    vstr.push_back("/// @param cbuffer The contracted integrals buffer.");
    
    vstr.push_back("/// @param " + label + " The contracted integrals buffer.");
    
    for (const auto& tint : t3c::get_geom_hrr_integrals(integral))
    {
        label = t3c::get_hrr_index(tint);
        
        vstr.push_back("/// @param " + label + " The contracted integrals buffer.");
    }
    
    
    return vstr;
}

std::vector<std::string>
T3CGeomDocuDriver::_get_ket_coordinates_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
   
    vstr.push_back("/// @param factors The factors buffer.");
                           
    vstr.push_back("/// @param idx_cd The vector of distances R(CD) = C - D.");
  
    return vstr;
}

std::vector<std::string>
T3CGeomDocuDriver::_get_ket_recursion_variables_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param a_angmom The angular momentum on center A.");
    
    return vstr;
}
