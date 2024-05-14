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

#include "t4c_geom_docs.hpp"

#include "file_stream.hpp"
#include "t4c_utils.hpp"

void
T4CGeomDocuDriver::write_doc_str(      std::ofstream& fstream,
                                 const SI4CIntegrals& geom_integrals,
                                 const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_compute_str(integral)});
    
    for (const auto& label : _get_buffers_str(geom_integrals, integral))
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
T4CGeomDocuDriver::_get_compute_str(const I4CIntegral& integral) const
{
    const auto bra_one = Tensor(integral[0]);
    
    const auto bra_two = Tensor(integral[1]);
        
    const auto ket_one = Tensor(integral[2]);
    
    const auto ket_two = Tensor(integral[3]);
    
    const auto integrand = integral.integrand();
    
    std::string label = "/// Computes ";

    label += t4c::prefixes_label(integral);
        
    label += "[" + bra_one.label() + bra_two.label() + "|G|";
   
    label += ket_one.label() + ket_two.label() + "]  integrals for arbitrary scalar operator G.";
    
    return label;
}

std::vector<std::string>
T4CGeomDocuDriver::_get_buffers_str(const SI4CIntegrals& geom_integrals,
                                    const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    auto label = t4c::get_geom_buffer_label(integral);
    
    vstr.push_back("/// - Parameter " + label + ": the integral geometrical derivatives buffer.");
    
    for (const auto& tint : geom_integrals)
    {
        auto label = t4c::get_geom_buffer_label(tint);
        
        vstr.push_back("/// - Parameter " + label + ": the primitive integrals buffer.");
    }

    return vstr;
}

std::vector<std::string>
T4CGeomDocuDriver::_get_recursion_variables_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        if (prefixes[0].shape().order() > 0)
        {
            vstr.push_back("/// - Parameter a_exp: the exponent on center A.");
        }
        
        if (prefixes[1].shape().order() > 0)
        {
            vstr.push_back("/// - Parameter b_exp: the exponent on center B.");
        }
        
        if (prefixes[2].shape().order() > 0)
        {
            vstr.push_back("/// - Parameter c_exps: the vector of exponents on center C.");
        }
        
        if (prefixes[3].shape().order() > 0)
        {
            vstr.push_back("/// - Parameter d_exps: the vector of exponents on center D.");
        }
    }
    
    return vstr;
}
