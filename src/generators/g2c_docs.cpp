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

#include "g2c_docs.hpp"

#include "file_stream.hpp"
#include "t2c_utils.hpp"
#include "string_formater.hpp"

void
G2CDocuDriver::write_doc_str(      std::ofstream&         fstream,
                             const I2CIntegral&           integral,
                             const bool                   use_rs) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_compute_str(integral, use_rs)});
    
    for (const auto& label : _get_distributor_str(use_rs))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_gto_blocks_str(integral))
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
G2CDocuDriver::_get_compute_str(const I2CIntegral& integral,
                                const bool         use_rs) const
{
    const auto bra = Tensor(integral[0]);
        
    const auto ket = Tensor(integral[1]);
    
    const auto [bra_prefix, ket_prefix] = t2c::prefixes_label(integral);
    
    const auto integrand = integral.integrand();
    
    auto label = "/// @brief Computes (" + bra_prefix + bra.label() + "|";
    
    if (integral.integrand().name() != "1")
    {
        if (use_rs)
        {
            label += "Erf(" + t2c::integrand_label(integral.integrand()) + ")|";
        }
        else
        {
            label += t2c::integrand_label(integral.integrand()) + "|";
        }
    }
    
    label += ket_prefix + ket.label() + ")  integrals for pair of basis functions on given grid.";
    
    return label;
}

std::vector<std::string>
G2CDocuDriver::_get_distributor_str(const bool use_rs) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param spher_buffer The spherical integrals buffer.");
    
    vstr.push_back("/// @param cart_buffer The Cartesian integrals buffer.");
    
    vstr.push_back("/// @param gcoords_x The Cartesian X coordinates of grid points.");
    
    vstr.push_back("/// @param gcoords_y The Cartesian Y coordinates of grid points.");
    
    vstr.push_back("/// @param gcoords_z The Cartesian Z coordinates of grid points.");
    
    vstr.push_back("/// @param gweights The weight of grid points.");
    
    return vstr;
}

std::vector<std::string>
G2CDocuDriver::_get_gto_blocks_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param bra_gto_block The basis functions block on bra side.");
        
    vstr.push_back("/// @param ket_gto_block The basis functions block on ket side.");
    
    return vstr;
}

std::vector<std::string>
G2CDocuDriver::_get_indices_str() const
{
    std::vector<std::string> vstr;

    vstr.push_back("/// @param bra_igto The index of basis function on bra side.");
        
    vstr.push_back("/// @param ket_igto The index of basis function on ket side.");
    
    return vstr;
}
