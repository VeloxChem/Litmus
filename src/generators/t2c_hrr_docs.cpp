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

#include "t2c_hrr_docs.hpp"

#include "file_stream.hpp"
#include "t2c_utils.hpp"

void
T2CHRRDocuDriver::write_doc_str(      std::ofstream& fstream,
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
    
    ost::write_code_lines(fstream, lines);
}

std::string
T2CHRRDocuDriver::_get_compute_str(const I2CIntegral& integral) const
{
    const auto bra = Tensor(integral[0]);
        
    const auto ket = Tensor(integral[1]);
    
    auto label = "/// @brief Computes contracted [" + bra.label() + "|X|";
    
    label += ket.label() + "]  integrals for set of data buffers.";
    
    return label;
}

std::vector<std::string>
T2CHRRDocuDriver::_get_buffers_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param cbuffer The contracted integrals buffer.");
    
    auto label = t2c::get_index_label(integral);
    
    vstr.push_back("/// @param " + label + " The index of integral in contracted integrals buffer.");
    
    for (const auto& tint : t2c::get_hrr_integrals(integral, integral))
    {
        label = t2c::get_index_label(tint);
        
        vstr.push_back("/// @param " + label + " The index of integral in contracted integrals buffer.");
    }

    return vstr;
}

std::vector<std::string>
T2CHRRDocuDriver::_get_coordinates_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param factors The contracted factors buffer.");

    return vstr;
}
