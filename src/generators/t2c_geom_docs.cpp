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

#include "t2c_geom_docs.hpp"

#include "file_stream.hpp"
#include "t2c_utils.hpp"

void
T2CGeomDocuDriver::write_doc_str(      std::ofstream&      fstream,
                                 const SI2CIntegrals&      geom_integrals,
                                 const I2CIntegral&        integral,
                                 const std::array<int, 3>& geom_drvs) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_compute_str(integral, geom_drvs)});
    
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
T2CGeomDocuDriver::_get_compute_str(const I2CIntegral&        integral,
                                    const std::array<int, 3>& geom_drvs) const
{
    const auto bra = Tensor(integral[0]);
        
    const auto ket = Tensor(integral[2]);
    
    std::string label = "/// Computes ";

    const auto pref_labels = t2c::prefixes_label(integral);

    if (geom_drvs[0] > 0)
    {
        label += "[" + pref_labels.first + bra.label() + "|R|";
    }
    else
    {
        label += "[" + bra.label() + "|R|";
    }
    
    if (geom_drvs[2] > 0)
    {
        label += pref_labels.second + ket.label() + "]  integrals for arbitrary operator R.";
    }
    else
    {
        label += ket.label() + "]  integrals for arbitrary operator R.";
    }
    
    return label;
}

std::vector<std::string>
T2CGeomDocuDriver::_get_buffers_str(const SI2CIntegrals& geom_integrals,
                                    const I2CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// - Parameter prim_buffer : the primitive integrals buffer.");
    
    auto label = t2c::get_index_label(integral);
    
    vstr.push_back("/// - Parameter " + label + ": the index of integral in primitive integrals buffer.");
    
    for (const auto& tint : t2c::get_integrals(integral))
    {
        label = t2c::get_index_label(tint);
        
        vstr.push_back("/// - Parameter " + label + ": the index of integral in primitive integrals buffer.");
    }

    return vstr;
}

std::vector<std::string>
T2CGeomDocuDriver::_get_recursion_variables_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        if (prefixes[0].shape().order() > 0)
        {
            vstr.push_back("/// - Parameter a_exp: the exponent on center A.");
        }
        
        if (prefixes[2].shape().order() > 0)
        {
            vstr.push_back("/// - Parameter b_exps: the vector of exponents on center B.");
        }
    }
    
    return vstr;
}
