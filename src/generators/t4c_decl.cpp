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

#include "t4c_decl.hpp"

#include "file_stream.hpp"
#include "t4c_utils.hpp"

void
T4CDeclDriver::write_func_decl(      std::ofstream& fstream,
                               const I4CIntegral&   integral,
                               const bool           diagonal,
                               const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_matrices_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_gto_pair_blocks_str(integral, diagonal))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_indices_str(integral, diagonal, terminus))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CDeclDriver::_get_matrices_str(const I4CIntegral& integral) const
{
    if (!integral.is_simple())
    {
        // TODO: Add derrivatives
    }
    
    std::vector<std::string> vstr;
    
    auto name = t4c::compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(name + "CMatrix* fock,");
    
    vstr.push_back(spacer + "const CMatrix* density,");
    
    return vstr;
}

std::vector<std::string>
T4CDeclDriver::_get_gto_pair_blocks_str(const I4CIntegral& integral,
                                        const bool         diagonal) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    if (diagonal)
    {
        vstr.push_back(spacer + "const CGtoPairBlock& gto_pair_block,");
    }
    else
    {
        vstr.push_back(spacer + "const CGtoPairBlock& bra_gto_pair_block,");
        
        vstr.push_back(spacer + "const CGtoPairBlock& ket_gto_pair_block,");
    }
    
    return vstr;
}

std::vector<std::string>
T4CDeclDriver::_get_indices_str(const I4CIntegral& integral,
                                const bool         diagonal,
                                const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    auto name = t4c::compute_func_name(integral) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(spacer + "const std::array<int, 2>& bra_indices,");
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    vstr.push_back(spacer + "const std::array<int, 2>& ket_indices) -> void" + tsymbol);
    
    return vstr;
}
