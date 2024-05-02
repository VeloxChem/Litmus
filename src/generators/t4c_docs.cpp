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

#include "t4c_docs.hpp"

#include "file_stream.hpp"
#include "t4c_utils.hpp"
#include "string_formater.hpp"


void
T4CDocuDriver::write_doc_str(      std::ofstream& fstream,
                             const I4CIntegral&   integral,
                             const bool           diagonal) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_compute_str(integral, diagonal)});
    
    for (const auto& label : _get_matrices_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }

    for (const auto& label : _get_gto_pair_blocks_str(integral, diagonal))
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
T4CDocuDriver::_get_compute_str(const I4CIntegral& integral,
                                const bool         diagonal) const
{
    const auto bra_one = Tensor(integral[0]);
    
    const auto bra_two = Tensor(integral[1]);
        
    const auto ket_one = Tensor(integral[2]);
    
    const auto ket_two = Tensor(integral[3]);
    
    const auto integrand = integral.integrand();
    
    auto label = "/// Computes (" + bra_one.label() + bra_two.label();
    
    label +=  "|" + t4c::integrand_label(integral.integrand()) + "|";
   
    label += ket_one.label() + ket_two.label() + ")  integrals for ";
        
    label += (diagonal) ? "GTOs pair block." : "two GTOs pair blocks.";
    
    return label;
}

std::vector<std::string>
T4CDocuDriver::_get_matrices_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// - Parameter distributor: the pointer to Fock matrix/matrices distributor.");
             
    return vstr;
}

std::vector<std::string>
T4CDocuDriver::_get_gto_pair_blocks_str(const I4CIntegral& integral,
                                        const bool         diagonal) const
{
    std::vector<std::string> vstr;
    
    if (diagonal)
    {
       vstr.push_back("/// - Parameter gto_pair_block: the GTOs pair block.");
    }
    else
    {
        vstr.push_back("/// - Parameter bra_gto_pair_block: the GTOs pair block on bra side.");
        
        vstr.push_back("/// - Parameter ket_gto_pair_block: the GTOs pair block on ket side.");
    }
        
    return vstr;
}

std::vector<std::string>
T4CDocuDriver::_get_indices_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// - Parameter bra_indices: the range [bra_first, bra_last) of GTOs on bra side.");
    
    vstr.push_back("/// - Parameter ket_indices: the range [ket_first, ket_last) of GTOs on ket side.");
    
    return vstr;
}
