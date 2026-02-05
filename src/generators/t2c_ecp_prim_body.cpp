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

#include "t2c_ecp_prim_body.hpp"

#include <algorithm>
#include <iostream>

#include "t2c_loc_ecp_driver.hpp"

void
T2CECPPrimFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                          const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    lines.push_back({1, 0, 2, "const auto nelems = pbuffer.number_of_active_elements();"});
    
    for (const auto& label : _get_factors_str(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    const auto components = integral.components<T1CPair, T1CPair>();
    
    const auto ncomps = static_cast<int>(components.size());
    
    std::vector<R2CDist> rec_dists;

    for (const auto& component : components)
    {
        rec_dists.push_back(_get_vrr_recursion(component));
    }

//    for (const auto& label : _get_buffers_str(rec_dists, integral))
//    {
//
//        lines.push_back({1, 0, 2, label});
//    }
//        
//    const auto kcomps = t2c::number_of_cartesian_components(integral[1]);
//    
//    if ((integral[0] == 0) || (integral[1] == 0))
//    {
//        const std::array<int, 2> rec_range({0, ncomps});
//
//        for (const auto& label : _get_buffers_str(integral, components, rec_range))
//        {
//
//            lines.push_back({1, 0, 2, label});
//        }
//        
//        _add_recursion_loop(lines, integral, components, rec_range);
//    }
//    else
//    {
//        const auto nblocks = ncomps / kcomps;
//        
//        for (int i = 0; i < nblocks; i++)
//        {
//            const std::array<int, 2> rec_range({i * kcomps, (i + 1) * kcomps});
//            
//            for (const auto& label : _get_buffers_str(integral, components, rec_range))
//            {
//                lines.push_back({1, 0, 2, label});
//            }
//            
//            _add_recursion_loop(lines, integral, components, {i * kcomps, (i + 1) * kcomps});
//            
//            if (i < (ncomps - 1))  lines.push_back({0, 0, 1, ""});;
//        }
//    }
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CECPPrimFuncBodyDriver::_get_factors_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (integral[0] > 0)
    {
        vstr.push_back("// Set up R(RA) distances");

        vstr.push_back("auto ra_x = factors.data(8);");
        
        vstr.push_back("auto ra_y = factors.data(9);");

        vstr.push_back("auto ra_z = factors.data(10);");
    }
    
    if (integral[1] > 0)
    {
        vstr.push_back("// Set up R(RB) distances");

        vstr.push_back("auto rb_x = factors.data(8);");
        
        vstr.push_back("auto rb_y = factors.data(9);");

        vstr.push_back("auto rb_z = factors.data(10);");
    }
    
    if ((integral[0] + integral[1]) > 1)
    {
        vstr.push_back("// Set up inverted 1/2 xi");

        vstr.push_back("auto fxi = factors.data(11);");
    }
            
    return vstr;
}

R2CDist
T2CECPPrimFuncBodyDriver::_get_vrr_recursion(const T2CIntegral& integral) const
{
    R2CDist rdist;

    if (integral.integrand().name() == "U_L")
    {
        T2CLocalECPDriver ecp_drv;
        
        if (integral[0].order() > 0)
        {
            rdist = ecp_drv.apply_bra_vrr(R2CTerm(integral));
        }
        else
        {
            rdist = ecp_drv.apply_ket_vrr(R2CTerm(integral));
        }
    }
   
    rdist.simplify();
    
    return rdist;
}
