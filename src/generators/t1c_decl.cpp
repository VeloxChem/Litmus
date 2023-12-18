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

#include "t1c_decl.hpp"

#include "file_stream.hpp"

void
T1CDeclDriver::write_func_decl(      std::ofstream& fstream,
                               const int            angmom,
                               const int            gdrv,
                               const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_gto_str(angmom, gdrv))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_vars_str(angmom, gdrv, terminus))
    {
        if  (label.find(";") == std::string::npos)
        {
            lines.push_back({0, 0, 1, label});
        }
        else
        {
            lines.push_back({0, 0, 2, label});
        }
    }
    
    ost::write_code_lines(fstream, lines);
}

std::string
T1CDeclDriver::_func_name(const int angmom,
                          const int gdrv) const
{
    return "getGeom" + std::to_string(gdrv) + "ValuesRec" + Tensor(angmom).label();
}

std::vector<std::string>
T1CDeclDriver::_get_gto_str(const int angmom,
                            const int gdrv) const
{
    std::vector<std::string> vstr;
    
    const auto label = _func_name(angmom, gdrv) + "(const CGtoBlock&            gto_block," ;
    
    vstr.push_back(label); 
    
    return vstr;
}

std::vector<std::string>
T1CDeclDriver::_get_vars_str(const int  angmom,
                             const int  gdrv,
                             const bool terminus) const
{
    std::vector<std::string> vstr;
   
    const auto nsize = _func_name(angmom, gdrv).size() + 1;
   
    vstr.push_back(std::string(nsize, ' ') + "const std::vector<double>&  grid_coords_x,");
    
    vstr.push_back(std::string(nsize, ' ') + "const std::vector<double>&  grid_coords_y,");
    
    vstr.push_back(std::string(nsize, ' ') + "const std::vector<double>&  grid_coords_z,");
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    vstr.push_back(std::string(nsize, ' ') + "const std::vector<int64_t>& gtos_mask) -> CMatrix" + tsymbol);

    return vstr;
}
