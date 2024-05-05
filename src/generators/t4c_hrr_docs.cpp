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

#include "t4c_hrr_docs.hpp"

#include "file_stream.hpp"
#include "t4c_utils.hpp"

void
T4CHrrDocuDriver::write_ket_doc_str(      std::ofstream& fstream,
                                    const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_ket_compute_str(integral)});
    
    // TODO: Add special variables here
    
    for (const auto& label : _get_ket_buffers_str(integral))
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

void
T4CHrrDocuDriver::write_bra_doc_str(      std::ofstream& fstream,
                                    const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_bra_compute_str(integral)});
    
    // TODO: Add special variables here
    
    for (const auto& label : _get_bra_buffers_str(integral))
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
T4CHrrDocuDriver::_get_ket_compute_str(const I4CIntegral& integral) const
{
    const auto ket_one = Tensor(integral[2]);
    
    const auto ket_two = Tensor(integral[3]);
    
    const auto integrand = integral.integrand();
    
    auto label = "/// Computes (XX|" + t4c::integrand_label(integral.integrand()) + "|";
   
    label += ket_one.label() + ket_two.label() + ")  integrals for set of data buffers.";
    
    return label;
}

std::vector<std::string>
T4CHrrDocuDriver::_get_ket_buffers_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto label = t4c::get_hrr_buffer_label(integral, true);
    
    vstr.push_back("/// - Parameter " + label + ": the contracted integrals buffer.");
    
    for (const auto& tint : t4c::get_ket_hrr_integrals(integral))
    {
        auto label = t4c::get_hrr_buffer_label(tint, true);
        
        vstr.push_back("/// - Parameter " + label + ": the contracted integrals buffer.");
    }

    return vstr;
}

std::vector<std::string>
T4CHrrDocuDriver::_get_ket_coordinates_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
   
    vstr.push_back("/// - Parameter cd_x: the vector of Cartesian X distances R(CD) = C - D.");
                       
    vstr.push_back("/// - Parameter cd_y: the vector of Cartesian Y distances R(CD) = C - D.");
                                      
    vstr.push_back("/// - Parameter cd_z: the vector of Cartesian Z distances R(CD) = C - D.");
        
    return vstr;
}

std::vector<std::string>
T4CHrrDocuDriver::_get_ket_recursion_variables_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// - Parameter a_angmom: the angular momentum on center A.");
    
    vstr.push_back("/// - Parameter b_angmom: the angular momentum on center B.");
                       
    return vstr;
}

std::string
T4CHrrDocuDriver::_get_bra_compute_str(const I4CIntegral& integral) const
{
    const auto bra_one = Tensor(integral[0]);
    
    const auto bra_two = Tensor(integral[1]);
    
    const auto integrand = integral.integrand();
    
    auto label = "/// Computes (" +  bra_one.label() + bra_two.label() + "|";
   
    label += t4c::integrand_label(integral.integrand()) + "XX)  integrals for set of data buffers.";
    
    return label;
}

std::vector<std::string>
T4CHrrDocuDriver::_get_bra_buffers_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto label = t4c::get_hrr_buffer_label(integral, false);
    
    vstr.push_back("/// - Parameter " + label + ": the contracted integrals buffer.");
    
    for (const auto& tint : t4c::get_bra_hrr_integrals(integral))
    {
        auto label = t4c::get_hrr_buffer_label(tint, false);
        
        vstr.push_back("/// - Parameter " + label + ": the contracted integrals buffer.");
    }

    return vstr;
}

std::vector<std::string>
T4CHrrDocuDriver::_get_bra_coordinates_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
   
    vstr.push_back("/// - Parameter ab_x: the Cartesian X distance R(AB) = A - B.");
    
    vstr.push_back("/// - Parameter ab_y: the Cartesian Y distance R(AB) = A - B.");
    
    vstr.push_back("/// - Parameter ab_z: the Cartesian Z distance R(AB) = A - B.");
        
    return vstr;
}

std::vector<std::string>
T4CHrrDocuDriver::_get_bra_recursion_variables_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// - Parameter c_angmom: the angular momentum on center C.");
    
    vstr.push_back("/// - Parameter d_angmom: the angular momentum on center D.");
                       
    return vstr;
}

