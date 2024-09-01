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

#include "t4c_prim_docs.hpp"

#include "file_stream.hpp"
#include "t4c_utils.hpp"

void
T4CPrimDocuDriver::write_doc_str(      std::ofstream& fstream,
                                 const I4CIntegral&   integral) const
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
T4CPrimDocuDriver::_get_compute_str(const I4CIntegral& integral) const
{
    const auto bra_one = Tensor(integral[0]);
    
    const auto bra_two = Tensor(integral[1]);
        
    const auto ket_one = Tensor(integral[2]);
    
    const auto ket_two = Tensor(integral[3]);
    
    const auto integrand = integral.integrand();
    
    auto label = "/// Computes [" + bra_one.label() + bra_two.label();
    
    label +=  "|" + t4c::integrand_label(integral.integrand()) + "|";
   
    label += ket_one.label() + ket_two.label() + "]  integrals for set of data buffers.";
    
    return label;
}

std::vector<std::string>
T4CPrimDocuDriver::_get_buffers_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param pbuffer The primitive integrals buffer.");
    
    auto label = t4c::get_index_label(integral);
    
    vstr.push_back("/// @param " + label + " The index of integral in primitive integrals buffer.");
        
    for (const auto& tint : t4c::get_vrr_integrals(integral))
    {
        label = t4c::get_index_label(tint);
        
        vstr.push_back("/// @param " + label + " The primitive integrals buffer.");
    }

    return vstr;
}

std::vector<std::string>
T4CPrimDocuDriver::_get_coordinates_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// @param factors The primitive factors buffer.");
    
    if (integral[1] > 0)
    {
        vstr.push_back("/// @param idx_wp The vector of distances R(WP) = W - P.");
                       
        vstr.push_back("/// @param r_pb The Cartesiandistances R(PB) = P - B.");
    }
    
    if ((integral[1] == 0) && (integral[3] > 0))
    {
        vstr.push_back("/// @param idx_qd The vector of distances R(QD) = Q - D.");
        
        vstr.push_back("/// @param idx_wq The vector of distances R(WQ) = W - Q.");
    }
   
    return vstr;
}

std::vector<std::string>
T4CPrimDocuDriver::_get_recursion_variables_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
//    if ((integral[1] > 1) && (integral[3] == 0))
//    {
//        vstr.push_back("/// - Parameter a_exp: the exponent on center A.");
//
//        vstr.push_back("/// - Parameter b_exp: the exponent on center B.");
//    }
//
//    if ((integral[1] == 0) && (integral[3] > 1))
//    {
//        vstr.push_back("/// - Parameter c_exps: the vector of exponents on center C.");
//
//        vstr.push_back("/// - Parameter d_exps: the vector of exponents on center D.");
//    }
    
    if ((integral[1] + integral[3]) > 1)
    {
        vstr.push_back("/// @param a_exp The exponent on center A.");
                       
        vstr.push_back("/// @param b_exp The exponent on center B.");
    }
   
    if ((integral[1] + integral[3]) == 0)
    {
        vstr.push_back("/// @param idx_ovl The index of combined overlap factors.");
        
        vstr.push_back("/// @param bf_data The Boys function data.");
        
        vstr.push_back("/// @param idx_bvals The index of Boys function data.");
    }
    
    return vstr;
}
