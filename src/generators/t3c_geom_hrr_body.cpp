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

#include "t3c_geom_hrr_body.hpp"

#include "t3c_utils.hpp"
#include "t2c_utils.hpp"
#include "t3c_geom_100_eri_driver.hpp"
#include "t3c_geom_010_eri_driver.hpp"
#include "string_formater.hpp"

void
T3CGeomHrrFuncBodyDriver::write_bra_func_body(      std::ofstream& fstream,
                                              const I3CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    lines.push_back({1, 0, 2, "const auto nelems = cbuffer.number_of_active_elements();"});
    
    lines.push_back({1, 0, 2, "const auto ccomps = tensor::number_of_cartesian_components(std::array<int, 1>{c_angmom,});"});
    
    lines.push_back({1, 0, 2, "const auto dcomps = tensor::number_of_cartesian_components(std::array<int, 1>{d_angmom,});"});
    
    lines.push_back({1, 0, 1, "for (int i = 0; i < ccomps; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 1, "for (int j = 0; j < dcomps; j++)"});
    
    lines.push_back({2, 0, 1, "{"});
    
    const auto components = integral.components<T1CPair, T2CPair>();
    
    std::vector<R3CDist> rec_dists;

    for (const auto& component : components)
    {
        rec_dists.push_back(_get_bra_hrr_recursion(component));
    }
    
    for (const auto& label : _get_bra_buffers_str(rec_dists, integral))
    {
        lines.push_back({3, 0, 2, label});
    }
   
    const auto bcomps = t2c::number_of_cartesian_components(integral[0]);

    const auto geom_orders = integral.prefixes_order();
    
    lines.push_back({3, 0, 2, "/// set up bra offset for " + t3c::get_hrr_buffer_label(integral, false)});
    
    if (geom_orders == std::vector<int>({1, 0, 0}))
    {
        lines.push_back({3, 0, 2, _get_full_bra_offset_def(integral)});
    }
    else
    {
        lines.push_back({3, 0, 2, _get_bra_offset_def(integral)});
    }

    if (geom_orders == std::vector<int>({1, 0, 0}))
    {
        for (int i = 0; i < 3; i++)
        {
            const std::array<int, 2> rec_range({i * bcomps, (i + 1) * bcomps});
            
            for (const auto& label : _get_bra_buffers_str(integral, components, rec_range))
            {
                lines.push_back({3, 0, 2, label});
            }
        
            _add_bra_recursion_loop(lines, integral, components, {i * bcomps, (i + 1) * bcomps});

            if (i < 2) lines.push_back({0, 0, 1, ""});;
        }
    }

    lines.push_back({2, 0, 1, "}"});
    
    lines.push_back({1, 0, 1, "}"});
   
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

R3CDist
T3CGeomHrrFuncBodyDriver::_get_bra_hrr_recursion(const T3CIntegral& integral) const
{
    R3CDist rdist;
    
    const auto geom_order = integral.prefixes_order();
    
    if (geom_order == std::vector<int>({1, 0, 0}))
    {
        if (integral.integrand().name() == "1/|r-r'|")
        {
            T3CGeom100ElectronRepulsionDriver eri_drv;
    
            if (integral[0].order() > 0)
            {
                rdist = eri_drv.apply_bra_hrr(R3CTerm(integral));
            }
        }
    }
    
    rdist.simplify();
    
    return rdist;
}

std::vector<std::string>
T3CGeomHrrFuncBodyDriver::_get_bra_buffers_str(const std::vector<R3CDist>& rec_dists,
                                               const I3CIntegral&          integral) const
{
    std::vector<std::string> vstr;
    
    const auto gorders = integral.prefixes_order();
    
    auto tints = t3c::get_bra_geom_integrals(integral);
    
    for (const auto& tint : tints)
    {
        auto label = "cbuffer.data(";
        
        vstr.push_back("/// Set up components of auxilary buffer : " + tint.label());
        
        if (gorders == std::vector<int>({1, 0, 0}))
        {
            vstr.push_back(_get_full_bra_offset_def(tint));
        }
        else
        {
            vstr.push_back(_get_bra_offset_def(tint));
        }
        
        int index = 0;
        
        for (const auto& tcomp : tint.components<T1CPair, T2CPair>())
        {
            std::string line;
            
            if (gorders == std::vector<int>({1, 0, 0}))
            {
                line += "auto " + _get_full_bra_component_label(tcomp) + " = " + label;
                
                line += _get_full_bra_offset_label(tint) + " + "  + std::to_string(index) + " * ccomps * dcomps);";
            }
            else
            {
                line += "auto " + _get_bra_component_label(tcomp) + " = " + label;
                
                line += _get_bra_offset_label(tint) + " + "  + std::to_string(index) + " * ccomps * dcomps);";
            }
            
            vstr.push_back(fstr::lowercase(line));
            
            index++;
        }
    }
    
    return vstr;
}

std::vector<std::string>
T3CGeomHrrFuncBodyDriver::_get_bra_buffers_str(const I3CIntegral&        integral,
                                               const VT3CIntegrals&      components,
                                               const std::array<int, 2>& rec_range) const
{
    std::vector<std::string> vstr;
    
    auto label = "cbuffer.data(";
    
    vstr.push_back("/// Set up " + std::to_string(rec_range[0]) + "-" + std::to_string(rec_range[1]) +
                       " components of targeted buffer : " + label);
    
    const auto gorders = integral.prefixes_order();
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        if (gorders == std::vector<int>({1, 0, 0}))
        {
            const auto line = "auto " + _get_full_bra_component_label(components[i]) + " = " + label;
            
            vstr.push_back(line  + _get_full_bra_offset_label(integral) + " + " + std::to_string(i) + " * ccomps * dcomps);");
        }
        else
        {
            const auto line = "auto " + _get_bra_component_label(components[i]) + " = " + label;
            
            vstr.push_back(line  + _get_bra_offset_label(integral) + " + " + std::to_string(i) + " * ccomps * dcomps);");
        }
    }
    
    return vstr;
}

std::string
T3CGeomHrrFuncBodyDriver::_get_tensor_label(const T3CIntegral& integral) const
{
    std::string label;
    
    if (integral.integrand().name() == "1/|r-r'|") label = "g";
    
    return label;
}

std::string
T3CGeomHrrFuncBodyDriver::_get_bra_component_label(const T3CIntegral& integral) const
{
    std::string label = _get_tensor_label(integral);
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        label += "_" + prefixes[0].label();
    }
    
    label += "_" + integral[0].label();
    
    return label;
}

std::string
T3CGeomHrrFuncBodyDriver::_get_full_bra_component_label(const T3CIntegral& integral) const
{
    std::string label = _get_tensor_label(integral);
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        label += "_" + prefixes[0].label() + "_" + prefixes[1].label();
        
        label += "_" + prefixes[2].label();
    }
    
    label += "_" + integral[0].label() + "_" +  integral[1].label();
    
    return label;
}

std::string
T3CGeomHrrFuncBodyDriver::_get_bra_offset_def(const I3CIntegral& integral) const
{
    const auto tlabel = std::to_string(integral.components<T1CPair, T2CPair>().size());
    
    auto label = "const auto " + _get_bra_offset_label(integral) + " = ";
    
    label += t3c::get_hrr_index(integral) +  " + i * dcomps + j;";
    
    return fstr::lowercase(label);
}

std::string
T3CGeomHrrFuncBodyDriver::_get_full_bra_offset_def(const I3CIntegral& integral) const
{
    const auto tlabel = std::to_string(integral.components<T1CPair, T2CPair>().size());
    
    auto label = "const auto " + _get_full_bra_offset_label(integral) + " = ";
    
    label += t3c::get_full_hrr_index(integral, false) +  " + i * dcomps + j;";
    
    return fstr::lowercase(label);
}

std::string
T3CGeomHrrFuncBodyDriver::_get_bra_offset_label(const I3CIntegral& integral) const
{
    const auto bra_one = Tensor(integral[0]);
    
    std::string label;
    
    if (const auto geom_orders = integral.prefixes_order(); !geom_orders.empty())
    {
        label = "_geom_" + std::to_string(geom_orders[0]) + std::to_string(geom_orders[1]);
    }
    
    label = bra_one.label() + label + "_off";
    
    return fstr::lowercase(label);
}

std::string
T3CGeomHrrFuncBodyDriver::_get_full_bra_offset_label(const I3CIntegral& integral) const
{
    const auto bra_one = Tensor(integral[0]);
    
    std::string label;
    
    if (const auto geom_orders = integral.prefixes_order(); !geom_orders.empty())
    {
        label += "_geom_" + std::to_string(geom_orders[0]) + std::to_string(geom_orders[1]);
        
        label += std::to_string(geom_orders[2]);
    }
    
    label = bra_one.label() + label + "_off";
    
    return fstr::lowercase(label);
}

void
T3CGeomHrrFuncBodyDriver::_add_bra_recursion_loop(      VCodeLines&         lines,
                                                  const I3CIntegral&        integral,
                                                  const VT3CIntegrals&      components,
                                                  const std::array<int, 2>& rec_range) const
{
    std::vector<R3CDist> rec_dists;
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        rec_dists.push_back(_get_bra_hrr_recursion(components[i]));
    }
    
    // set up recursion loop
    
    const auto var_str = _get_bra_pragma_str(integral, rec_dists);
    
    lines.push_back({3, 0, 1, "#pragma omp simd aligned(" + var_str + " : 64)"});
    
    lines.push_back({3, 0, 1, "for (size_t k = 0; k < nelems; k++)"});
    
    lines.push_back({3, 0, 1, "{"});
    
    for (size_t i = 0; i < rec_dists.size(); i++)
    {
        if (i < (rec_dists.size() - 1))
        {
            lines.push_back({4, 0, 2, _get_bra_code_line(rec_dists[i])});
        }
        else
        {
            lines.push_back({4, 0, 1, _get_bra_code_line(rec_dists[i])});
        }
    }
    
    lines.push_back({3, 0, 1, "}"});
}

std::string
T3CGeomHrrFuncBodyDriver::_get_bra_pragma_str(const I3CIntegral&          integral,
                                              const std::vector<R3CDist>& rec_distributions) const
{
    std::set<std::string> tlabels;
    
    const auto gorders = integral.prefixes_order();
    
    for (const auto& rdist : rec_distributions)
    {
        auto tint = rdist.root().integral();
        
        if (!gorders.empty())
        {
            tlabels.insert(_get_full_bra_component_label(tint));
        }
        else
        {
            tlabels.insert(_get_bra_component_label(tint));
        }
        
        for (size_t i = 0; i < rdist.terms(); i++)
        {
            auto tint = rdist[i].integral();
            
            if (!gorders.empty())
            {
                tlabels.insert(_get_full_bra_component_label(tint));
            }
            else
            {
                tlabels.insert(_get_bra_component_label(tint));
            }
        }
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
T3CGeomHrrFuncBodyDriver::_get_bra_code_line(const R3CDist& rec_distribution) const
{
    auto tint = rec_distribution.root().integral();
    
    std::string line;
    
    if (tint.prefixes_order() == std::vector<int>({1, 0, 0}))
    {
        line = _get_full_bra_component_label(tint) + "[k] = ";
    }
    else
    {
        line = _get_bra_component_label(tint) + "[k] = ";
    }
    
    for (size_t i = 0; i < rec_distribution.terms(); i++)
    {
        auto tint = rec_distribution[i].integral();
        
        line += _get_bra_rterm_code(rec_distribution[i], i == 0);
    }
    
    return line + ";";
}

std::string
T3CGeomHrrFuncBodyDriver::_get_bra_rterm_code(const R3CTerm& rec_term,
                                              const bool     is_first) const
{
    const auto pre_fact = rec_term.prefactor();
        
    auto plabel = pre_fact.label();
        
    if (plabel == "1.0")  plabel = "";
        
    if (plabel == "-1.0") plabel = "-";
    
    if (plabel.size() > 1) plabel += " * ";
    
    auto tint = rec_term.integral();
    
    const auto gorders = tint.prefixes_order();
    
    if (!gorders.empty())
    {
        if ((gorders[0] + gorders[1]) > 0)
        {
            plabel += _get_full_bra_component_label(tint) + "[k]";
        }
        else
        {
            plabel += _get_bra_component_label(tint) + "[k]";
        }
    }
    else
    {
        plabel += _get_full_bra_component_label(tint) + "[k]";
    }
    
    for (const auto& fact : rec_term.factors())
    {
        plabel+= " * " + fact.label();
            
        //if (fact.order() > 0) plabel += "[k]";
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

void
T3CGeomHrrFuncBodyDriver::write_ket_func_body(      std::ofstream& fstream,
                                              const I3CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    lines.push_back({1, 0, 2, "const auto nelems = cbuffer.number_of_active_elements();"});
    
    lines.push_back({1, 0, 2, "const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom,});"});
    
    lines.push_back({1, 0, 2, "// Set up R(CD) distances"});

    lines.push_back({1, 0, 2, "auto cd_x = factors.data(idx_cd);"});

    lines.push_back({1, 0, 2, "auto cd_y = factors.data(idx_cd + 1);"});

    lines.push_back({1, 0, 2, "auto cd_z = factors.data(idx_cd + 2);"});
    
    lines.push_back({1, 0, 1, "for (int i = 0; i < acomps; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    const auto components = integral.components<T1CPair, T2CPair>();
    
    std::vector<R3CDist> rec_dists;

    for (const auto& component : components)
    {
        rec_dists.push_back(_get_ket_hrr_recursion(component));
    }
    
    for (const auto& label : _get_ket_buffers_str(rec_dists, integral))
    {
        lines.push_back({2, 0, 2, label});
    }
    
    const auto bcomps = t2c::number_of_cartesian_components(integral[1]);
    
    const auto kcomps = t2c::number_of_cartesian_components(integral[2]);
    
    lines.push_back({2, 0, 2, "/// set up bra offset for " + t3c::get_hrr_buffer_label(integral, true)});
    
    lines.push_back({2, 0, 2, _get_ket_offset_def(integral)});
    
    const auto gorders = integral.prefixes_order();
    
    if (gorders == std::vector<int>({0, 1, 0}))
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < bcomps; j++)
            {
                const std::array<int, 2> rec_range({j * kcomps, (j + 1) * kcomps});
                
                for (const auto& label : _get_ket_geom_buffers_str(integral, components, rec_range, i, bcomps * kcomps))
                {
                    lines.push_back({2, 0, 2, label});
                }
                
                _add_ket_recursion_loop(lines, integral, components, {j * kcomps, (j + 1) * kcomps}, i, bcomps * kcomps);
                
                if (j < (bcomps - 1))  lines.push_back({0, 0, 1, ""});;
            }
        }
    }

    
    lines.push_back({1, 0, 1, "}"});
   
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

R3CDist
T3CGeomHrrFuncBodyDriver::_get_ket_hrr_recursion(const T3CIntegral& integral) const
{
    R3CDist rdist;
    
    const auto geom_order = integral.prefixes_order();
    
    if (geom_order == std::vector<int>({0, 1, 0}))
    {
        if (integral.integrand().name() == "1/|r-r'|")
        {
            T3CGeom010ElectronRepulsionDriver eri_drv;
    
            if (integral[1].order() > 0)
            {
                rdist = eri_drv.apply_ket_hrr(R3CTerm(integral));
            }
            else
            {
                rdist = eri_drv.apply_ket_aux_hrr(R3CTerm(integral));
            }
        }
    }
    
    rdist.simplify();
    
    return rdist;
}

std::vector<std::string>
T3CGeomHrrFuncBodyDriver::_get_ket_buffers_str(const std::vector<R3CDist>& rec_dists,
                                               const I3CIntegral&          integral) const
{
    std::vector<std::string> vstr;
    
    if (integral[1] == 0)
    {
        for (const auto& tint : t3c::get_geom_hrr_integrals(integral))
        {
            
            auto label = "cbuffer.data(";
            
            vstr.push_back("/// Set up components of auxilary buffer : " + tint.label());
            
            vstr.push_back(_get_ket_offset_def(tint));
            
            int index = 0;
            
            for (const auto& tcomp : tint.components<T1CPair, T2CPair>())
            {
                //if (_find_integral(rec_dists, tcomp))
                //{
                    auto line = "auto " + _get_ket_component_label(tcomp) + " = " + label;
                    
                    line +=  _get_ket_offset_label(tint) + " + "  + std::to_string(index) + ");";
                    
                    vstr.push_back(fstr::lowercase(line));
                //}
                
                index++;
            }
        }
    }
    else
    {
        for (const auto& tint : t3c::get_geom_hrr_integrals(integral))
        {
            const auto gorders = tint.prefixes_order();
            
            if (!gorders.empty())
            {
                std::string label = "cbuffer.data(";
                
                vstr.push_back("/// Set up components of auxilary buffer : " + tint.label());
                
                vstr.push_back(_get_ket_offset_def(tint));
                
                const auto bcomps = t2c::number_of_cartesian_components(tint[1]);
                
                const auto kcomps = t2c::number_of_cartesian_components(tint[2]);
                
                int index = 0;
                
                for (const auto& tcomp : tint.components<T1CPair, T2CPair>())
                {
                    const auto line = "auto " + _get_ket_component_label(tcomp) + " = " + label;
                    
                    const std::string glabel = std::to_string((index / (bcomps * kcomps)) * bcomps * kcomps) + " * acomps";
                    
                    vstr.push_back(line + _get_ket_offset_label(tint) + " + " + glabel + " + " + std::to_string(index % (bcomps * kcomps)) + ");");
                   
                    index++;
                }
            }
            else
            {
                std::string label = "cbuffer.data(";
                
                vstr.push_back("/// Set up components of auxilary buffer : " + tint.label());
                
                vstr.push_back(_get_ket_offset_def(tint));
                
                int index = 0;
                
                for (const auto& tcomp : tint.components<T1CPair, T2CPair>())
                {
                    //if (_find_integral(rec_dists, tcomp))
                    //{
                        auto line = "auto " + _get_ket_component_label(tcomp) + " = " + label;
                        
                        line +=  _get_ket_offset_label(tint) + " + "  + std::to_string(index) + ");";
                        
                        vstr.push_back(fstr::lowercase(line));
                    //}
                    
                    index++;
                }
            }
        }
    }

    return vstr;
}

std::string
T3CGeomHrrFuncBodyDriver::_get_ket_component_label(const T3CIntegral& integral) const
{
    std::string label = _get_tensor_label(integral);
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        label += "_" + prefixes[1].label() + "_" + prefixes[2].label();
    }
    
    label += "_" + integral[1].label() + "_" +  integral[2].label();
    
    return label;
}

std::string
T3CGeomHrrFuncBodyDriver::_get_ket_offset_def(const I3CIntegral& integral) const
{
    // const auto tlabel = std::to_string(integral.components<T2CPair, T2CPair>().size());
    
    const auto bcomps = t2c::number_of_cartesian_components(integral[1]);
    
    const auto kcomps = t2c::number_of_cartesian_components(integral[2]);
    
    const auto tlabel = std::to_string(bcomps * kcomps);
    
    auto label = "const auto " + _get_ket_offset_label(integral) + " = ";
    
    label += t3c::get_hrr_index(integral) + " + i * " + tlabel + ";";
    
    return fstr::lowercase(label);
}

std::string
T3CGeomHrrFuncBodyDriver::_get_ket_offset_label(const I3CIntegral& integral) const
{
    const auto ket_one = Tensor(integral[1]);
    
    const auto ket_two = Tensor(integral[2]);
    
    std::string label;
    
    if (const auto geom_orders = integral.prefixes_order(); !geom_orders.empty())
    {
        label = "_geom_" + std::to_string(geom_orders[1]) + std::to_string(geom_orders[2]);
    }
    
    label = ket_one.label() + ket_two.label() + label  + "_off";
    
    return fstr::lowercase(label);
}

std::vector<std::string>
T3CGeomHrrFuncBodyDriver::_get_ket_geom_buffers_str(const I3CIntegral&        integral,
                                                    const VT3CIntegrals&      components,
                                                    const std::array<int, 2>& rec_range,
                                                    const int                 ket_index,
                                                    const int                 ket_components) const
{
    std::vector<std::string> vstr;
    
    auto label = "cbuffer.data(";
    
    vstr.push_back("/// Set up " + std::to_string(rec_range[0]) + "-" + std::to_string(rec_range[1]) +
                       " components of targeted buffer : " + label);
    
    const int koff = ket_index * ket_components;
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        const auto line = "auto " + _get_ket_component_label(components[koff + i]) + " = " + label;
        
        const std::string glabel = std::to_string(koff) + " * acomps ";
        
        vstr.push_back(line + _get_ket_offset_label(integral) + " + " + glabel + " + " + std::to_string(i) + ");");
    }
    
    return vstr;
}

void
T3CGeomHrrFuncBodyDriver::_add_ket_recursion_loop(      VCodeLines&         lines,
                                                  const I3CIntegral&        integral,
                                                  const VT3CIntegrals&      components,
                                                  const std::array<int, 2>& rec_range,
                                                  const int                 ket_index,
                                                  const int                 ket_components) const
{
    std::vector<R3CDist> rec_dists;
    
    const auto koff = ket_index * ket_components;
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        rec_dists.push_back(_get_ket_hrr_recursion(components[i + koff]));
    }
    
    // set up recursion loop
    
    const auto var_str = _get_ket_pragma_str(integral, rec_dists);
    
    lines.push_back({2, 0, 1, "#pragma omp simd aligned(" + var_str + " : 64)"});
    
    lines.push_back({2, 0, 1, "for (size_t k = 0; k < nelems; k++)"});
    
    lines.push_back({2, 0, 1, "{"});
    
    for (size_t i = 0; i < rec_dists.size(); i++)
    {
        if (i < (rec_dists.size() - 1))
        {
            lines.push_back({3, 0, 2, _get_ket_code_line(rec_dists[i])});
        }
        else
        {
            lines.push_back({3, 0, 1, _get_ket_code_line(rec_dists[i])});
        }
    }
    
    lines.push_back({2, 0, 1, "}"});
}

std::string
T3CGeomHrrFuncBodyDriver::_get_ket_pragma_str(const I3CIntegral&          integral,
                                              const std::vector<R3CDist>& rec_distributions) const
{
    std::set<std::string> tlabels;
    
    for (const auto& rdist : rec_distributions)
    {
        auto tint = rdist.root().integral();
        
        tlabels.insert(_get_ket_component_label(tint));
        
        for (size_t i = 0; i < rdist.terms(); i++)
        {
            auto tint = rdist[i].integral();
            
            tlabels.insert(_get_ket_component_label(tint));
            
            for (const auto& fact : rdist[i].factors())
            {
                if (fact.order() > 0) tlabels.insert(fact.label());
            }
        }
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
T3CGeomHrrFuncBodyDriver::_get_ket_code_line(const R3CDist& rec_distribution) const
{
    auto tint = rec_distribution.root().integral();
    
    std::string line = _get_ket_component_label(tint) + "[k] = ";
    
    for (size_t i = 0; i < rec_distribution.terms(); i++)
    {
        auto tint = rec_distribution[i].integral();
        
        line += _get_ket_rterm_code(rec_distribution[i], i == 0);
    }
    
    return line + ";";
}

std::string
T3CGeomHrrFuncBodyDriver::_get_ket_rterm_code(const R3CTerm& rec_term,
                                              const bool     is_first) const
{
    const auto pre_fact = rec_term.prefactor();
        
    auto plabel = pre_fact.label();
        
    if (plabel == "1.0")  plabel = "";
        
    if (plabel == "-1.0") plabel = "-";
    
    if (plabel.size() > 1) plabel += " * ";
    
    auto tint = rec_term.integral();
    
    plabel += _get_ket_component_label(tint) + "[k]";
        
    for (const auto& fact : rec_term.factors())
    {
        plabel+= " * " + fact.label();
            
        if (fact.order() > 0) plabel += "[k]";
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
