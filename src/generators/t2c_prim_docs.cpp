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

#include "t2c_prim_docs.hpp"

#include "file_stream.hpp"
#include "t2c_utils.hpp"

void
T2CPrimDocuDriver::write_doc_str(      std::ofstream& fstream,
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
T2CPrimDocuDriver::_get_compute_str(const I2CIntegral& integral) const
{
    const auto bra = Tensor(integral[0]);
        
    const auto ket = Tensor(integral[1]);
    
    const auto [bra_prefix, ket_prefix] = t2c::prefixes_label(integral);
    
    const auto integrand = integral.integrand();
    
    auto label = "/// Computes primitive [" + bra_prefix + bra.label() + "|";
    
    if (integral.integrand().name() != "1")
    {
        label += t2c::integrand_label(integral.integrand()) + "|";
    }
    
    label += ket_prefix + ket.label() + "]  integrals for set of data buffers.";
    
    return label;
}

std::vector<std::string>
T2CPrimDocuDriver::_get_buffers_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto label = t2c::get_buffer_label(integral, "prim");
    
    vstr.push_back("/// - Parameter prim_buffer : the primitive integrals buffer.");
    
    label = t2c::get_index_label(integral);
    
    vstr.push_back("/// - Parameter " + label + ": the index of integral in primitive integrals buffer.");
    
    for (const auto& tint : t2c::get_integrals(integral))
    {
        auto label = t2c::get_index_label(tint);
        
        vstr.push_back("/// - Parameter " + label + ": the index of integral in primitive integrals buffer.");
    }

    return vstr;
}

std::vector<std::string>
T2CPrimDocuDriver::_get_coordinates_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
   
    if (integral[0] > 0)
    {
        vstr.push_back("/// - Parameter pa_x: the vector of Cartesian X distances R(PA) = P - A.");
                       
        vstr.push_back("/// - Parameter pa_y: the vector of Cartesian Y distances R(PA) = P - A.");
                                      
        vstr.push_back("/// - Parameter pa_z: the vector of Cartesian Z distances R(PA) = P - A.");
    }
   
    if ((integral[0] == 0) && (integral[1] > 0))
    {
        vstr.push_back("/// - Parameter pb_x: the vector of Cartesian X distances R(PB) = P - B.");
                       
        vstr.push_back("/// - Parameter pb_y: the vector of Cartesian Y distances R(PB) = P - B.");
                                      
        vstr.push_back("/// - Parameter pb_z: the vector of Cartesian Z distances R(PB) = P - B.");
    }
    
    if ((integral[0] + integral[1]) > 0)
    {
        if (integral.integrand().name() == "A")
        {
            vstr.push_back("/// - Parameter pc_x: the vector of Cartesian X distances R(PC) = P - C.");
                           
            vstr.push_back("/// - Parameter pc_y: the vector of Cartesian Y distances R(PC) = P - C.");
                                          
            vstr.push_back("/// - Parameter pc_z: the vector of Cartesian Z distances R(PC) = P - C.");
        }
    }
    
    if ((integral[0] + integral[1]) == 0)
    {
        if (integral.integrand().name() == "A")
        {
            vstr.push_back("/// - Parameter bf_values: the vector of Boys function values.");
        }
        else
        {
            vstr.push_back("/// - Parameter ab_x: the vector of Cartesian X distances R(AB) = A - B.");
                           
            vstr.push_back("/// - Parameter ab_y: the vector of Cartesian Y distances R(AB) = A - B.");
                                          
            vstr.push_back("/// - Parameter ab_z: the vector of Cartesian Z distances R(AB) = A - B.");
        }
    }
   
    return vstr;
}

std::vector<std::string>
T2CPrimDocuDriver::_get_recursion_variables_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (((integral[0] + integral[1])  != 1) || (integral.integrand().name() == "T"))
    {
        vstr.push_back("/// - Parameter a_exp: the GTO exponent on center A.");
        
        vstr.push_back("/// - Parameter b_exps: the vector of GTOs exponents on center B.");
    }
   
    if ((integral[0] + integral[1]) == 0)
    {
        if (integral.integrand().name() == "1")
        {
            vstr.push_back("/// - Parameter a_norm: the GTO normalization factor on center A.");
            
            vstr.push_back("/// - Parameter b_norms: the vector of GTOs normalization factors on center B.");
        }
    }
    
    return vstr;
}
