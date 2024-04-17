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

#include "t2c_docs.hpp"

#include "file_stream.hpp"
#include "t2c_utils.hpp"
#include "string_formater.hpp"

void
T2CDocuDriver::write_doc_str(      std::ofstream&         fstream,
                             const I2CIntegral&           integral,
                             const std::pair<bool, bool>& rec_form,
                             const bool                   diagonal) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_compute_str(integral, diagonal)});
    
    for (const auto& label : _get_matrices_str(integral, rec_form))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    // TODO: Add special variables here
    
    for (const auto& label : _get_gto_blocks_str(integral, diagonal))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_distributor_variables_str(integral, diagonal))
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
T2CDocuDriver::_get_compute_str(const I2CIntegral& integral,
                                const bool         diagonal) const
{
    const auto bra = Tensor(integral[0]);
        
    const auto ket = Tensor(integral[1]);
    
    const auto [bra_prefix, ket_prefix] = t2c::prefixes_label(integral);
    
    const auto integrand = integral.integrand();
    
    auto label = "/// Computes <" + bra_prefix + bra.label() + "|";
    
    if (integral.integrand().name() != "1")
    {
        label += t2c::integrand_label(integral.integrand()) + "|";
    }
    
    label += ket_prefix + ket.label() + ">  integrals for ";
        
    label += (diagonal) ? "GTOs block." : "pair of GTOs blocks.";
    
    return label;
}

std::vector<std::string>
T2CDocuDriver::_get_matrices_str(const I2CIntegral&           integral,
                                 const std::pair<bool, bool>& rec_form) const
{
    std::vector<std::string> vstr;
    
    if (integral.is_simple())
    {
        for (const auto& label : t2c::integrand_labels(integral, "matrix"))
        {
            vstr.push_back("/// - Parameter " + label + ": the pointer to matrix for storage of integrals.");
        }
    }
    else
    {
        // TODO: Add derrivatives
    }
             
    return vstr;
}

std::vector<std::string>
T2CDocuDriver::_get_gto_blocks_str(const I2CIntegral& integral,
                                   const bool         diagonal) const
{
    std::vector<std::string> vstr;
    
    if (diagonal)
    {
       vstr.push_back("/// - Parameter gto_block: the GTOs block.");
    }
    else
    {
        vstr.push_back("/// - Parameter bra_gto_block: the GTOs block on bra side.");
        
        vstr.push_back("/// - Parameter ket_gto_block: the GTOs block on ket side.");
    }
        
    return vstr;
}

std::vector<std::string>
T2CDocuDriver::_get_distributor_variables_str(const I2CIntegral& integral,
                                              const bool         diagonal) const
{
    std::vector<std::string> vstr;
    
    if (!diagonal)
    {
        if (integral[0] != integral[1])
        {
            vstr.push_back("/// - Parameter ang_order: the flag for matching angular order between matrix and pair of GTOs blocks.");
        }
        else
        {
            vstr.push_back("/// - Parameter mat_type: the matrix type.");
        }
    }
    
    return vstr;
}

std::vector<std::string>
T2CDocuDriver::_get_indices_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// - Parameter bra_indices: the range [bra_first, bra_last) of GTOs on bra side.");
    
    vstr.push_back("/// - Parameter ket_indices: the range [ket_first, ket_last) of GTOs on ket side.");
    
    return vstr;
}


