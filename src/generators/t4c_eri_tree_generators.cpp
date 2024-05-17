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

#include "t4c_eri_tree_generators.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "file_stream.hpp"

#include "t4c_utils.hpp"

void
T4CCallTreeGenerator::generate(const std::string& label,
                               const int          max_ang_mom) const
{
    std::string fname = "CallTreeFile.tmp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    SI4CIntegrals integrals;
    
    if (_is_available(label))
    {
        for (int i = 0; i <= max_ang_mom; i++)
        {
            for (int j = i; j <= max_ang_mom; j++)
            {
                for (int k = 0; k <= max_ang_mom; k++)
                {
                    for (int l = k; l <= max_ang_mom; l++)
                    {
                        integrals.insert(_get_integral(label, {i, j, k, l}));
                    }
                }
            }
        }
    }

    auto lines = VCodeLines();
    
    for (const auto& integral : integrals)
    {
        lines.push_back({0, 0, 1, "#include \"" + _file_name(integral) + ".hpp\""});
    }
    
    fstream << std::endl;
    
    for (const auto& integral : integrals)
    {
        lines.push_back({0, 0, 1, "if ((bra_angmoms == std::array<int, 2>({" + std::to_string(integral[0]) + ", " + std::to_string(integral[1]) + "})) &&"});
        
        lines.push_back({0, 0, 1, "    (ket_angmoms == std::array<int, 2>({" + std::to_string(integral[2]) + ", " + std::to_string(integral[3]) + "})))"});
        
        lines.push_back({0, 0, 1, "{"});
        
        lines.push_back({1, 0, 2, t4c::namespace_label(integral) + "::" + t4c::compute_func_name(integral) +  "(distributor, bra_gto_pair_block, ket_gto_pair_block, bra_range, ket_range);"});
            
        lines.push_back({1, 0, 1, "return;"});
        
        lines.push_back({0, 0, 2, "}"});
    }
    
    ost::write_code_lines(fstream, lines);
    
    fstream.close();
}

bool
T4CCallTreeGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "electron repulsion") return true;
    
    return false;
}

I4CIntegral
T4CCallTreeGenerator::_get_integral(const std::string&        label,
                                    const std::array<int, 4>& ang_moms) const
{
    // bra and ket sides

    const auto bpair = I2CPair("GA", ang_moms[0], "GB", ang_moms[1]);

    const auto kpair = I2CPair("GC", ang_moms[2], "GD", ang_moms[3]);

    // electron repulsion integrals

    if (fstr::lowercase(label) == "electron repulsion")
    {
        return I4CIntegral(bpair, kpair, Operator("1/|r-r'|"), 0, {});
    }
    
    return I4CIntegral();
}

std::string
T4CCallTreeGenerator::_file_name(const I4CIntegral& integral) const
{
    std::string label = "Rec" + integral.label();
    
    return t4c::integral_label(integral) + label;
}
