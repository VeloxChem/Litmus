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

#include "t4c_diag_docs.hpp"

#include "file_stream.hpp"
#include "t4c_utils.hpp"
#include "string_formater.hpp"

void
T4CDiagDocuDriver::write_doc_str(      std::ofstream& fstream,
                                 const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
        
    lines.push_back({0, 0, 1, "/**"});
        
    lines.push_back({0, 0, 2, _get_compute_str(integral)});
    
    for (const auto& label : _get_vars_str())
    {
        lines.push_back({0, 1, 1, label});
    }
    
    lines.push_back({0, 0, 1, "*/"});
        
    ost::write_code_lines(fstream, lines);
}

std::string
T4CDiagDocuDriver::_get_compute_str(const I4CIntegral& integral) const
{
    const auto bra_a = Tensor(integral[0]);
    
    const auto bra_b = Tensor(integral[1]);
        
    const auto ket_a = Tensor(integral[2]);
    
    const auto ket_b = Tensor(integral[3]);
        
    const auto integrand = integral.integrand();
            
    auto label = " Evaluates <"  + bra_a.label() + bra_b.label() + "|";
        
    label += t4c::integrand_label(integral.integrand());
        
    label += "|" + ket_a.label()  + ket_b.label()+ ">  integrals for given ";
        
    label += "GTOs pair block.";
    
    return label;
    
}

std::vector<std::string>
T4CDiagDocuDriver::_get_vars_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("@param gto_pair_block the GTOs pair block for bra and ket sides.");
    
    vstr.push_back("@return the vector with largest Cartesian component of electron repulsion integrals.");
    
    return vstr;
}
