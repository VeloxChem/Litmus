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

#include "t2c_decl.hpp"

#include "file_stream.hpp"
#include "t2c_utils.hpp"

void
T2CDeclDriver::write_func_decl(      std::ofstream&         fstream,
                               const I2CIntegral&           integral,
                               const std::pair<bool, bool>& rec_form,
                               const bool                   use_rs,
                               const bool                   terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "template <class T>"});
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_distributor_str(integral, rec_form, use_rs))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_gto_blocks_str(integral, rec_form, use_rs))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_indices_str(integral, rec_form, use_rs, terminus))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CDeclDriver::_get_distributor_str(const I2CIntegral&           integral,
                                    const std::pair<bool, bool>& rec_form,
                                    const bool                   use_rs) const
{
    std::vector<std::string> vstr;
    
    auto name = t2c::compute_func_name(integral, rec_form, use_rs) + "(";
    
    vstr.push_back(name + "T& distributor,");
    
    if (use_rs)
    {
        const auto spacer = std::string(name.size(), ' ');
        
        vstr.push_back(spacer + "const std::vector<double>& omegas,");
    }
    
    return vstr;
}

std::vector<std::string>
T2CDeclDriver::_get_gto_blocks_str(const I2CIntegral&           integral,
                                   const std::pair<bool, bool>& rec_form,
                                   const bool                   use_rs) const
{
    std::vector<std::string> vstr;
    
    auto name = t2c::compute_func_name(integral, rec_form, use_rs) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    vstr.push_back(spacer + "const CGtoBlock& bra_gto_block,");
        
    vstr.push_back(spacer + "const CGtoBlock& ket_gto_block,");
    
    return vstr;
}

std::vector<std::string>
T2CDeclDriver::_get_indices_str(const I2CIntegral&           integral,
                                const std::pair<bool, bool>& rec_form,
                                const bool                   use_rs,
                                const bool                   terminus) const
{
    std::vector<std::string> vstr;
    
    auto name = t2c::compute_func_name(integral, rec_form, use_rs) + "(";
    
    const auto spacer = std::string(name.size(), ' ');
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    vstr.push_back(spacer + "const std::pair<size_t, size_t>& bra_indices,");
        
    vstr.push_back(spacer + "const std::pair<size_t, size_t>& ket_indices," + tsymbol);
    
    vstr.push_back(spacer + "const bool bra_eq_ket) -> void" + tsymbol);
   
    return vstr;
}
