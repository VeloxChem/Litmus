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

#include "t2c_geom_body.hpp"

void
T2CGeomFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                       const SI2CIntegrals& geom_integrals,
                                       const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
//    lines.push_back({1, 0, 2, "const auto ndims = " +
//                              t4c::get_geom_buffer_label(integral) +
//                              ".number_of_columns();"});
//
//    for (const auto& label : _get_buffers_str(geom_integrals, integral))
//    {
//        lines.push_back({1, 0, 2, label});
//    }
//
//    // generate integrals group
//
//    const auto components = integral.components<T2CPair, T2CPair>();
//
//    const auto rgroup = _generate_integral_group(components, integral);
//
//    // set up dimensions needed for splitting recursion into subblocks
//
//    const auto ncomps = static_cast<int>(components.size());
//
//    const auto acomps = t2c::number_of_cartesian_components(integral[0]);
//
//    const auto bcomps = t2c::number_of_cartesian_components(integral[1]);
//
//    const auto ccomps = t2c::number_of_cartesian_components(integral[2]);
//
//    const auto dcomps = t2c::number_of_cartesian_components(integral[3]);
//
//    const int ocomps = ncomps / (acomps * bcomps * ccomps * dcomps);
//
//    if ((acomps * bcomps * ccomps * dcomps) == 1)
//    {
//        _add_recursion_loop(lines, rgroup, integral, {0, ncomps});
//    }
//    else
//    {
//        int rcomps = dcomps;
//
//        if (rcomps == 1) rcomps = ccomps;
//
//        if (rcomps == 1) rcomps = bcomps;
//
//        if (rcomps == 1) rcomps = acomps;
//
//        int blkoff = 0;
//
//        for (int i = 0; i < ocomps; i++)
//        {
//            for (int j = 0; j < ncomps / (ocomps * rcomps); j++)
//            {
//                _add_recursion_loop(lines, rgroup, integral, {blkoff, blkoff + rcomps});
//
//                blkoff += rcomps;
//            }
//        }
//    }
        
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}
