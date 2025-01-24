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

#include "t3c_docs.hpp"

#include "file_stream.hpp"
#include "t3c_utils.hpp"
#include "string_formater.hpp"

void
T3CDocuDriver::write_doc_str(      std::ofstream& fstream,
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

std::string
T3CDocuDriver::_get_compute_str(const I3CIntegral& integral) const
{
    const auto bra_one = Tensor(integral[0]);
            
    const auto ket_one = Tensor(integral[1]);
    
    const auto ket_two = Tensor(integral[2]);
    
    const auto integrand = integral.integrand();
    
    auto label = "/// @brief Computes (" + bra_one.label();
    
    label +=  "|" + t3c::integrand_label(integral.integrand()) + "|";
   
    label += ket_one.label() + ket_two.label() + ")  integrals for ";
        
    label += "basis functions block and basis function pairs block.";
    
    return label;
}

std::vector<std::string>
T3CDocuDriver::_get_matrices_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param distributor The pointer to integrals distributor.");
             
    return vstr;
}

std::vector<std::string>
T3CDocuDriver::_get_gto_pair_blocks_str(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param bra_gto_block The basis functions block on bra side.");
    
    vstr.push_back("/// @param ket_gto_pair_block The basis function pairs block on ket side.");
        
    return vstr;
}

std::vector<std::string>
T3CDocuDriver::_get_indices_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param bra_indices The range [bra_first, bra_last) of basis functions on bra side.");
        
    return vstr;
}

