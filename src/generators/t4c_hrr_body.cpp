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

#include "t4c_hrr_body.hpp"

#include <algorithm>

#include "t4c_utils.hpp"
#include "t2c_utils.hpp"
#include "t4c_vrr_eri_driver.hpp"
#include "t4c_hrr_eri_driver.hpp"
#include "string_formater.hpp"

void
T4CHrrFuncBodyDriver::write_ket_func_body(      std::ofstream& fstream,
                                          const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    lines.push_back({1, 0, 2, "const auto nelems = cbuffer.number_of_active_elements();"});
    
    lines.push_back({1, 0, 2, "const auto acomps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom,});"});
    
    lines.push_back({1, 0, 2, "const auto bcomps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom,});"});
    
    lines.push_back({1, 0, 2, "// Set up R(CD) distances"});

    lines.push_back({1, 0, 2, "auto cd_x = factors.data(idx_cd);"});

    lines.push_back({1, 0, 2, "auto cd_y = factors.data(idx_cd + 1);"});

    lines.push_back({1, 0, 2, "auto cd_z = factors.data(idx_cd + 2);"});
    
    lines.push_back({1, 0, 1, "for (int i = 0; i < acomps; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 1, "for (int j = 0; j < bcomps; j++)"});
    
    lines.push_back({2, 0, 1, "{"});
    
    const auto components = integral.components<T2CPair, T2CPair>();
    
    std::vector<R4CDist> rec_dists;
    
    for (const auto& component : components)
    {
        rec_dists.push_back(_get_ket_hrr_recursion(component));
    }
    
    for (const auto& label : _get_ket_buffers_str(rec_dists, integral))
    {
        lines.push_back({3, 0, 2, label});
    }
    
    const auto bcomps = t2c::number_of_cartesian_components(integral[2]);
    
    const auto kcomps = t2c::number_of_cartesian_components(integral[3]);
    
    lines.push_back({3, 0, 2, "/// set up bra offset for " + t4c::get_hrr_buffer_label(integral, true)});
    
    lines.push_back({3, 0, 2, _get_ket_offset_def(integral)});
    
    for (int i = 0; i < bcomps; i++)
    {
        const std::array<int, 2> rec_range({i * kcomps, (i + 1) * kcomps});
        
        for (const auto& label : _get_ket_buffers_str(integral, components, rec_range))
        {
            lines.push_back({3, 0, 2, label});
        }
        
        _add_ket_recursion_loop(lines, integral, components, {i * kcomps, (i + 1) * kcomps});
        
        if (i < (bcomps - 1))  lines.push_back({0, 0, 1, ""});;
    }
    
    lines.push_back({2, 0, 1, "}"});
    
    lines.push_back({1, 0, 1, "}"});
   
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

void
T4CHrrFuncBodyDriver::write_bra_func_body(      std::ofstream& fstream,
                                          const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    lines.push_back({1, 0, 2, "const auto nelems = cbuffer.number_of_active_elements();"});
    
    lines.push_back({1, 0, 2, "const auto ccomps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom,});"});
    
    lines.push_back({1, 0, 2, "const auto dcomps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom,});"});
    
    lines.push_back({1, 0, 2, "// set up R(AB) distances"});
        
    lines.push_back({1, 0, 2, "const auto xyz = r_ab.coordinates();"});
        
    lines.push_back({1, 0, 2, "const auto ab_x = xyz[0];"});
        
    lines.push_back({1, 0, 2, "const auto ab_y = xyz[1];"});
        
    lines.push_back({1, 0, 2, "const auto ab_z = xyz[2];"});
    
    lines.push_back({1, 0, 1, "for (int i = 0; i < ccomps; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 1, "for (int j = 0; j < dcomps; j++)"});
    
    lines.push_back({2, 0, 1, "{"});
    
    const auto components = integral.components<T2CPair, T2CPair>();
    
    std::vector<R4CDist> rec_dists;
    
    for (const auto& component : components)
    {
        rec_dists.push_back(_get_bra_hrr_recursion(component));
    }
    
    for (const auto& label : _get_bra_buffers_str(rec_dists, integral))
    {
        lines.push_back({3, 0, 2, label});
    }
    
    const auto bcomps = t2c::number_of_cartesian_components(integral[0]);
    
    const auto kcomps = t2c::number_of_cartesian_components(integral[1]);
    
    lines.push_back({3, 0, 2, "/// set up bra offset for " + t4c::get_hrr_buffer_label(integral, false)});
    
    lines.push_back({3, 0, 2, _get_bra_offset_def(integral)});
    
    for (int i = 0; i < bcomps; i++)
    {
        const std::array<int, 2> rec_range({i * kcomps, (i + 1) * kcomps});
        
        for (const auto& label : _get_bra_buffers_str(integral, components, rec_range))
        {
            lines.push_back({3, 0, 2, label});
        }
        
        _add_bra_recursion_loop(lines, integral, components, {i * kcomps, (i + 1) * kcomps});
        
        if (i < (bcomps - 1))  lines.push_back({0, 0, 1, ""});;
    }
    
    lines.push_back({2, 0, 1, "}"});
    
    lines.push_back({1, 0, 1, "}"});
   
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CHrrFuncBodyDriver::_get_ket_buffers_str(const std::vector<R4CDist>& rec_dists,
                                           const I4CIntegral&          integral) const
{
    std::vector<std::string> vstr;
    
    for (const auto& tint : t4c::get_ket_hrr_integrals(integral))
    {
        auto label = (integral[2] == 1) ? "pbuffer.data(" : "cbuffer.data(";
        
        vstr.push_back("/// Set up components of auxilary buffer : " + tint.label());
        
        vstr.push_back(_get_ket_offset_def(tint));
        
        int index = 0;
        
        for (const auto& tcomp : tint.components<T2CPair, T2CPair>())
        {
            if (_find_integral(rec_dists, tcomp))
            {
                auto line = "auto " + _get_ket_component_label(tcomp) + " = " + label;
                
                line +=  _get_ket_offset_label(tint) + " + "  + std::to_string(index) + ");";
                
                vstr.push_back(fstr::lowercase(line));
            }
            
            index++;
        }
    }
    
    return vstr;
}

std::vector<std::string>
T4CHrrFuncBodyDriver::_get_ket_buffers_str(const I4CIntegral&        integral,
                                           const VT4CIntegrals&      components,
                                           const std::array<int, 2>& rec_range) const
{
    std::vector<std::string> vstr;
    
    auto label = "cbuffer.data(";
    
    vstr.push_back("/// Set up " + std::to_string(rec_range[0]) + "-" + std::to_string(rec_range[1]) +
                       " components of targeted buffer : " + label);
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        const auto line = "auto " + _get_ket_component_label(components[i]) + " = " + label;
        
        vstr.push_back(line + _get_ket_offset_label(integral) + " + " + std::to_string(i) + ");");
    }
    
    return vstr;
}

std::vector<std::string>
T4CHrrFuncBodyDriver::_get_bra_buffers_str(const std::vector<R4CDist>& rec_dists,
                                           const I4CIntegral&          integral) const
{
    std::vector<std::string> vstr;
    
    for (const auto& tint : t4c::get_bra_hrr_integrals(integral))
    {
        auto label = "cbuffer.data(";
        
        vstr.push_back("/// Set up components of auxilary buffer : " + tint.label());
        
        vstr.push_back(_get_bra_offset_def(tint));
        
        int index = 0;
        
        for (const auto& tcomp : tint.components<T2CPair, T2CPair>())
        {
            if (_find_integral(rec_dists, tcomp))
            {
                auto line = "auto " + _get_bra_component_label(tcomp) + " = " + label;
                
                line += _get_bra_offset_label(tint) + " + "  + std::to_string(index) + " * ccomps * dcomps);";
                
                vstr.push_back(fstr::lowercase(line));
            }
            
            index++;
        }
    }
    
    return vstr;
}

std::vector<std::string>
T4CHrrFuncBodyDriver::_get_bra_buffers_str(const I4CIntegral&        integral,
                                           const VT4CIntegrals&      components,
                                           const std::array<int, 2>& rec_range) const
{
    std::vector<std::string> vstr;
    
    auto label = "cbuffer.data(";
    
    vstr.push_back("/// Set up " + std::to_string(rec_range[0]) + "-" + std::to_string(rec_range[1]) +
                       " components of targeted buffer : " + label);
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        const auto line = "auto " + _get_bra_component_label(components[i]) + " = " + label;
        
        vstr.push_back(line  + _get_bra_offset_label(integral) + " + " + std::to_string(i) + " * ccomps * dcomps);");
    }
    
    return vstr;
}

bool
T4CHrrFuncBodyDriver::_find_integral(const std::vector<R4CDist>& rec_dists,
                                     const T4CIntegral&          integral) const
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
T4CHrrFuncBodyDriver::_get_tensor_label(const I4CIntegral& integral) const
{
    std::string label;
    
    if (integral.integrand().name() == "1/|r-r'|") label = "g";
    
    return label;
}

std::string
T4CHrrFuncBodyDriver::_get_tensor_label(const T4CIntegral& integral) const
{
    std::string label;
    
    if (integral.integrand().name() == "1/|r-r'|") label = "g";
    
    return label;
}

void
T4CHrrFuncBodyDriver::_add_ket_recursion_loop(      VCodeLines&         lines,
                                              const I4CIntegral&        integral,
                                              const VT4CIntegrals&      components,
                                              const std::array<int, 2>& rec_range) const
{
    std::vector<R4CDist> rec_dists;
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        rec_dists.push_back(_get_ket_hrr_recursion(components[i]));
    }
    
    // set up recursion loop
    
    const auto var_str = _get_ket_pragma_str(integral, rec_dists);
    
    lines.push_back({3, 0, 1, "#pragma omp simd aligned(" + var_str + " : 64)"});
    
    lines.push_back({3, 0, 1, "for (size_t k = 0; k < nelems; k++)"});
    
    lines.push_back({3, 0, 1, "{"});
    
    for (size_t i = 0; i < rec_dists.size(); i++)
    {
        if (i < (rec_dists.size() - 1))
        {
            lines.push_back({4, 0, 2, _get_ket_code_line(rec_dists[i])});
        }
        else
        {
            lines.push_back({4, 0, 1, _get_ket_code_line(rec_dists[i])});
        }
    }
    
    lines.push_back({3, 0, 1, "}"});
}

std::string
T4CHrrFuncBodyDriver::_get_ket_pragma_str(const I4CIntegral&          integral,
                                          const std::vector<R4CDist>& rec_distributions) const
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

R4CDist
T4CHrrFuncBodyDriver::_get_ket_hrr_recursion(const T4CIntegral& integral) const
{
    R4CDist rdist;
    
    if (integral.integrand().name() == "1/|r-r'|")
    {
        T4CHrrElectronRepulsionDriver eri_drv;

        if (integral[2].order() > 0)
        {
            rdist = eri_drv.apply_ket_hrr(R4CTerm(integral));
        }
    }
    
    rdist.simplify();
    
    return rdist;
}

std::string
T4CHrrFuncBodyDriver::_get_ket_code_line(const R4CDist& rec_distribution) const
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
T4CHrrFuncBodyDriver::_get_ket_rterm_code(const R4CTerm& rec_term,
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

void
T4CHrrFuncBodyDriver::_add_bra_recursion_loop(      VCodeLines&         lines,
                                              const I4CIntegral&        integral,
                                              const VT4CIntegrals&      components,
                                              const std::array<int, 2>& rec_range) const
{
    std::vector<R4CDist> rec_dists;
    
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
T4CHrrFuncBodyDriver::_get_bra_pragma_str(const I4CIntegral&          integral,
                                          const std::vector<R4CDist>& rec_distributions) const
{
    std::set<std::string> tlabels;
    
    for (const auto& rdist : rec_distributions)
    {
        auto tint = rdist.root().integral();
        
        tlabels.insert(_get_bra_component_label(tint));
        
        for (size_t i = 0; i < rdist.terms(); i++)
        {
            auto tint = rdist[i].integral();
            
            tlabels.insert(_get_bra_component_label(tint));
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

R4CDist
T4CHrrFuncBodyDriver::_get_bra_hrr_recursion(const T4CIntegral& integral) const
{
    R4CDist rdist;
    
    if (integral.integrand().name() == "1/|r-r'|")
    {
        T4CHrrElectronRepulsionDriver eri_drv;

        if (integral[0].order() > 0)
        {
            rdist = eri_drv.apply_bra_hrr(R4CTerm(integral));
        }
    }
    
    rdist.simplify();
    
    return rdist;
}

std::string
T4CHrrFuncBodyDriver::_get_bra_code_line(const R4CDist& rec_distribution) const
{
    auto tint = rec_distribution.root().integral();
    
    std::string line = _get_bra_component_label(tint) + "[k] = ";
    
    for (size_t i = 0; i < rec_distribution.terms(); i++)
    {
        auto tint = rec_distribution[i].integral();
        
        line += _get_bra_rterm_code(rec_distribution[i], i == 0);
    }
    
    return line + ";";
}

std::string
T4CHrrFuncBodyDriver::_get_bra_rterm_code(const R4CTerm& rec_term,
                                          const bool     is_first) const
{
    const auto pre_fact = rec_term.prefactor();
        
    auto plabel = pre_fact.label();
        
    if (plabel == "1.0")  plabel = "";
        
    if (plabel == "-1.0") plabel = "-";
    
    if (plabel.size() > 1) plabel += " * ";
    
    auto tint = rec_term.integral();
    
    plabel += _get_bra_component_label(tint) + "[k]";
        
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

std::string
T4CHrrFuncBodyDriver::_get_component_label(const T4CIntegral& integral) const
{
    std::string label = _get_tensor_label(integral) + "_" + integral.label();
    
    label += "_" + std::to_string(integral.order());
    
    return label;
}

std::string
T4CHrrFuncBodyDriver::_get_ket_component_label(const T4CIntegral& integral) const
{
    std::string label = _get_tensor_label(integral) + "_" + integral[2].label() + "_" +  integral[3].label();
    
    return label;
}

std::string
T4CHrrFuncBodyDriver::_get_ket_offset_def(const I4CIntegral& integral) const
{
    const auto tlabel = std::to_string(integral.components<T2CPair, T2CPair>().size());
    
    auto label = "const auto " + _get_ket_offset_label(integral) + " = ";
    
    label += t4c::get_hrr_index(integral, true) + " + (i * bcomps + j) * " + tlabel + ";";
    
    return fstr::lowercase(label);
}

std::string
T4CHrrFuncBodyDriver::_get_ket_offset_label(const I4CIntegral& integral) const
{
    const auto ket_one = Tensor(integral[2]);
    
    const auto ket_two = Tensor(integral[3]);
    
    auto label = ket_one.label() + ket_two.label() + "_off";
    
    return fstr::lowercase(label);
}

std::string
T4CHrrFuncBodyDriver::_get_bra_component_label(const T4CIntegral& integral) const
{
    std::string label = _get_tensor_label(integral) + "_" + integral[0].label() + "_" +  integral[1].label();
    
    return label;
}

std::string
T4CHrrFuncBodyDriver::_get_bra_offset_def(const I4CIntegral& integral) const
{
    const auto tlabel = std::to_string(integral.components<T2CPair, T2CPair>().size());
    
    auto label = "const auto " + _get_bra_offset_label(integral) + " = ";
    
    label += t4c::get_hrr_index(integral, false) +  " + i * dcomps + j;";
    
    return fstr::lowercase(label);
}

std::string
T4CHrrFuncBodyDriver::_get_bra_offset_label(const I4CIntegral& integral) const
{
    const auto bra_one = Tensor(integral[0]);
    
    const auto bra_two = Tensor(integral[1]);
    
    auto label = bra_one.label() + bra_two.label() + "_off";
    
    return fstr::lowercase(label);
}

