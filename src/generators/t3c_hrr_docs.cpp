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

#include "t3c_hrr_docs.hpp"

#include "file_stream.hpp"
#include "t3c_utils.hpp"

void
T3CHrrDocuDriver::write_doc_str(      std::ofstream& fstream,
                                const I3CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_compute_str(integral)});
    
    // TODO: Add special variables here
    
    for (const auto& label : _get_buffers_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_coordinates_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_recursion_variables_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    ost::write_code_lines(fstream, lines);
}

std::string
T3CHrrDocuDriver::_get_compute_str(const I3CIntegral& integral) const
{
    const auto ket_one = Tensor(integral[1]);
    
    const auto ket_two = Tensor(integral[2]);
    
    const auto integrand = integral.integrand();
    
    auto label = "/// Computes (X|" + t3c::integrand_label(integral.integrand()) + "|";
   
    label += ket_one.label() + ket_two.label() + ")  integrals for set of data buffers.";
    
    return label;
}

std::vector<std::string>
T3CHrrDocuDriver::_get_buffers_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto label = t3c::get_hrr_index(integral);
    
    vstr.push_back("/// @param cbuffer The contracted integrals buffer.");
    
    vstr.push_back("/// @param " + label + " The contracted integrals buffer.");
    
    for (const auto& tint : t3c::get_hrr_integrals(integral))
    {
        label = t3c::get_hrr_index(tint);
        
        vstr.push_back("/// @param " + label + " The contracted integrals buffer.");
    }

    return vstr;
}

std::vector<std::string>
T3CHrrDocuDriver::_get_coordinates_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
   
    vstr.push_back("/// @param factors The factors buffer.");
                       
    vstr.push_back("/// @param idx_cd The vector of distances R(CD) = C - D.");
                                              
    return vstr;
}

std::vector<std::string>
T3CHrrDocuDriver::_get_recursion_variables_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param a_angmom The angular momentum on center A.");
                       
    return vstr;
}
