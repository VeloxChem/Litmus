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
#include "t2c_utils.hpp"

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

    for (const auto& label : _get_buffers_str(rec_dists, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
       
    const std::array<int, 2> rec_range({0, ncomps});
    
    for (const auto& label : _get_buffers_str(integral, components, rec_range))
    {
        lines.push_back({1, 0, 2, label});
    }
        
    _add_recursion_loop(lines, integral, components, rec_range);
    
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
        vstr.push_back("// Set up inverted 1/2xi");

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

std::vector<std::string>
T2CECPPrimFuncBodyDriver::_get_buffers_str(const std::vector<R2CDist>& rec_dists,
                                           const I2CIntegral&          integral) const
{
    std::vector<std::string> vstr;
    
    for (const auto& tint : t2c::get_integrals(integral))
    {
        vstr.push_back("// Set up components of auxiliary buffer : " + tint.label());

        const auto tlabel = _get_tensor_label(tint);
        
        int index = 0;
        
        for (const auto& tcomp : tint.components<T1CPair, T1CPair>())
        {
            if (_find_integral(rec_dists, tcomp))
            {
                const auto line = "auto " + _get_component_label(tcomp) + " = pbuffer.data(";
                
                if (index > 0)
                {
                    vstr.push_back(line + t2c::get_index_label(tint) + " + " + std::to_string(index) + ");");
                }
                else
                {
                    vstr.push_back(line + t2c::get_index_label(tint) + ");");
                }
            }

            index++;
        }
    }
    
    return vstr;
}

std::string
T2CECPPrimFuncBodyDriver::_get_tensor_label(const I2CIntegral& integral) const
{
    std::string label;
    
    if (integral.integrand().name() == "U_L") label = "tg";

    return label;
}

std::string
T2CECPPrimFuncBodyDriver::_get_tensor_label(const T2CIntegral& integral) const
{
    std::string label;
    
    if (integral.integrand().name() == "U_L") label = "tg";

    return label;
}

bool
T2CECPPrimFuncBodyDriver::_find_integral(const std::vector<R2CDist>& rec_dists,
                                         const T2CIntegral&          integral) const
{
    for (const auto& rdist : rec_dists)
    {
        for (const auto& tint : rdist.unique_integrals())
        {
            if (integral == tint) return true;
        }
    }
    
    return false;
}

std::string
T2CECPPrimFuncBodyDriver::_get_component_label(const T2CIntegral& integral) const
{
    std::string label = _get_tensor_label(integral) + "_" + integral.label();
    
    return label;
}

std::vector<std::string>
T2CECPPrimFuncBodyDriver::_get_buffers_str(const I2CIntegral&        integral,
                                           const VT2CIntegrals&      components,
                                           const std::array<int, 2>& rec_range) const
{
    std::vector<std::string> vstr;
    
    auto label = t2c::get_buffer_label(integral, "prim");
    
    if ((rec_range[1] - rec_range[0]) == static_cast<int>(components.size()))
    {
        vstr.push_back("// Set up components of targeted buffer : " + integral.label());
    }
    else
    {
        vstr.push_back("// Set up " + std::to_string(rec_range[0]) + "-" + std::to_string(rec_range[1]) +
                       " components of targeted buffer : " + integral.label());
    }
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        const auto line = "auto " + _get_component_label(components[i]) + " = pbuffer.data(";
        
        if (i > 0)
        {
            vstr.push_back(line + t2c::get_index_label(integral) + " + " + std::to_string(i) + ");");
        }
        else
        {
            vstr.push_back(line + t2c::get_index_label(integral) + ");");
        }
    }
    
    return vstr;
}

void
T2CECPPrimFuncBodyDriver::_add_recursion_loop(      VCodeLines&         lines,
                                              const I2CIntegral&        integral,
                                              const VT2CIntegrals&      components,
                                              const std::array<int, 2>& rec_range) const
{
    std::vector<R2CDist> rec_dists;
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        rec_dists.push_back(_get_vrr_recursion(components[i]));
    }
    
    // set up recursion loop
    
    const auto var_str = _get_pragma_str(integral, rec_dists);
    
    lines.push_back({1, 0, 1, "#pragma omp simd aligned(" + var_str + " : 64)"});
    
    lines.push_back({1, 0, 1, "for (size_t i = 0; i < nelems; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    //_get_factor_lines(lines, rec_dists);
    
    for (size_t i = 0; i < rec_dists.size(); i++)
    {
        if (i < (rec_dists.size() - 1))
        {
            lines.push_back({2, 0, 2, _get_code_line(rec_dists[i])});
        }
        else
        {
            lines.push_back({2, 0, 1, _get_code_line(rec_dists[i])});
        }
    }
    
    lines.push_back({1, 0, 1, "}"});
}

std::string
T2CECPPrimFuncBodyDriver::_get_pragma_str(const I2CIntegral&          integral,
                                          const std::vector<R2CDist>& rec_distributions) const
{
    std::set<std::string> tlabels;
    
    for (const auto& rdist : rec_distributions)
    {
        auto tint = rdist.root().integral();
        
        tlabels.insert(_get_component_label(tint));
        
        for (size_t i = 0; i < rdist.terms(); i++)
        {
            auto tint = rdist[i].integral();
            
            tlabels.insert(_get_component_label(tint));
            
            for (const auto& fact : rdist[i].factors())
            {
                if (fact.order() > 0) tlabels.insert(fact.label());
            }
        }
    }
    
    if (integral[0] > 0)
    {
        tlabels.insert("ra_x");
        
        tlabels.insert("ra_y");
        
        tlabels.insert("ra_z");
    }
    
    if (integral[1] > 0)
    {
        tlabels.insert("rb_x");
        
        tlabels.insert("rb_y");
        
        tlabels.insert("rb_z");
    }
    
    if ((integral[0] + integral[1]) > 1)
    {
        tlabels.insert("fxi");
    }
    
    std::string label;
    
    for (const auto& tlabel : tlabels)
    {
        label += tlabel + ", ";
    }
    
    if (label[label.size() - 2] == ',') label.erase(label.end() - 2);
    
    return label;
}

std::string
T2CECPPrimFuncBodyDriver::_get_code_line(const R2CDist& rec_distribution) const
{
    auto tint = rec_distribution.root().integral();
    
    std::string line = _get_component_label(tint) + "[i] = ";

    for (size_t i = 0; i < rec_distribution.terms(); i++)
    {
        auto tint = rec_distribution[i].integral();
        
        line += _get_rterm_code(rec_distribution[i], i == 0);
    }

    return line + ";";
}

std::string
T2CECPPrimFuncBodyDriver::_get_rterm_code(const R2CTerm& rec_term,
                                          const bool     is_first) const
{
    const auto pre_fact = rec_term.prefactor();
        
    auto plabel = pre_fact.label();
        
    if (plabel == "1.0")  plabel = "";
        
    if (plabel == "-1.0") plabel = "-";
    
    if (plabel.size() > 1) plabel += " * ";
    
    auto tint = rec_term.integral();
    
    plabel += _get_component_label(tint) + "[i]";

    for (const auto& fact : rec_term.factors())
    {
        if (fact.label() == "fxi_0")
        {
            plabel += " * fxi[i]"; 
        }
        else
        {
            plabel += " * " + fact.label();
        }
        
        if (fact.order() > 0) plabel += "[i]";
    }
    
    if (!is_first)
    {
        if (plabel[0] == '-')
        {
            plabel.insert(plabel.begin() + 1, 1, ' ');
            
            plabel = " " + plabel;
        }
        else
        {
            plabel = " + " + plabel;
        }
    }
        
    return plabel;
}
