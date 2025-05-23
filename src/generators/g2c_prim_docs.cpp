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

#include "g2c_prim_docs.hpp"

#include "file_stream.hpp"
#include "t2c_utils.hpp"

void
G2CPrimDocuDriver::write_doc_str(      std::ofstream& fstream,
                                 const I2CIntegral&   integral) const
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
G2CPrimDocuDriver::_get_compute_str(const I2CIntegral& integral) const
{
    const auto bra = Tensor(integral[0]);
        
    const auto ket = Tensor(integral[1]);
    
    const auto [bra_prefix, ket_prefix] = t2c::prefixes_label(integral);
    
    const auto integrand = integral.integrand();
    
    auto label = "/// @brief Computes primitive [" + bra_prefix + bra.label() + "|";
    
    if (integral.integrand().name() != "1")
    {
        label += t2c::integrand_label(integral.integrand()) + "|";
    }
    
    label += ket_prefix + ket.label() + "]  integrals for set of data buffers on given grid.";
    
    return label;
}

std::vector<std::string>
G2CPrimDocuDriver::_get_buffers_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param buffer The primitive integrals buffer.");
    
    auto label = t2c::get_index_label(integral);
    
    vstr.push_back("/// @param " + label + " The index of integral in primitive integrals buffer.");
    
    for (const auto& tint : t2c::get_integrals(integral))
    {
        label = t2c::get_index_label(tint);
        
        vstr.push_back("/// @param " + label + " The index of integral in primitive integrals buffer.");
    }

    return vstr;
}

std::vector<std::string>
G2CPrimDocuDriver::_get_coordinates_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (integral[0] > 0)
    {
        vstr.push_back("/// @param rpa_x The Cartesian X distance of R(PA) = P - A.");
        
        vstr.push_back("/// @param rpa_y The Cartesian Y distance of R(PA) = P - A.");
        
        vstr.push_back("/// @param rpa_z The Cartesian Z distance of R(PA) = P - A.");
    }
    
    if ((integral[0] ==  0)  && (integral[1] > 0))
    {
        vstr.push_back("/// @param rpb_x The Cartesian X distance of R(PB) = P - B.");
    
        vstr.push_back("/// @param rpb_y The Cartesian Y distance of R(PB) = P - B.");
        
        vstr.push_back("/// @param rpb_z The Cartesian Z distance of R(PB) = P - B.");
    }
   
    return vstr;
}

std::vector<std::string>
G2CPrimDocuDriver::_get_recursion_variables_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (_need_exponents(integral))
    {
        vstr.push_back("/// @param factor The combined exponential factor.");
    }
    
    return vstr;
}

bool
G2CPrimDocuDriver::_need_exponents(const I2CIntegral& integral) const
{
    return (integral[0] + integral[1]) > 1;
}
