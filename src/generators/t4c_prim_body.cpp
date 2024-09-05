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

#include "t4c_prim_body.hpp"

#include <algorithm>

#include "t4c_utils.hpp"
#include "t2c_utils.hpp"
#include "t4c_vrr_eri_driver.hpp"

void
T4CPrimFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                       const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    lines.push_back({1, 0, 2, "const auto nelems = pbuffer.number_of_active_elements();"});
    
    if ((integral[1] + integral[3]) > 1)
    {
        lines.push_back({1, 0, 2, "// Set up exponents"});

        lines.push_back({1, 0, 2, "auto c_exps = factors.data(0);"});
        
        lines.push_back({1, 0, 2, "auto d_exps = factors.data(1);"});
    }
    
    if (integral[1] > 0)
    {
        lines.push_back({1, 0, 2, "// Set up R(WP) distances"});

        lines.push_back({1, 0, 2, "auto wp_x = factors.data(idx_wp);"});

        lines.push_back({1, 0, 2, "auto wp_y = factors.data(idx_wp + 1);"});

        lines.push_back({1, 0, 2, "auto wp_z = factors.data(idx_wp + 2);"});
        
        lines.push_back({1, 0, 2, "// set up R(PB) distances"});

        lines.push_back({1, 0, 2, "const auto xyz = r_pb.coordinates();"});

        lines.push_back({1, 0, 2, "const auto pb_x = xyz[0];"});

        lines.push_back({1, 0, 2, "const auto pb_y = xyz[1];"});

        lines.push_back({1, 0, 2, "const auto pb_z = xyz[2];"});
    }
    
    if ((integral[1] == 0) && (integral[3] > 0))
    {
        lines.push_back({1, 0, 2, "// Set up R(QD) distances"});

        lines.push_back({1, 0, 2, "auto qd_x = factors.data(idx_qd);"});

        lines.push_back({1, 0, 2, "auto qd_y = factors.data(idx_qd + 1);"});

        lines.push_back({1, 0, 2, "auto qd_z = factors.data(idx_qd + 2);"});
        
        lines.push_back({1, 0, 2, "// Set up R(WQ) distances"});

        lines.push_back({1, 0, 2, "auto wq_x = factors.data(idx_wq);"});

        lines.push_back({1, 0, 2, "auto wq_y = factors.data(idx_wq + 1);"});

        lines.push_back({1, 0, 2, "auto wq_z = factors.data(idx_wq + 2);"});
    }
   
    const auto components = integral.components<T2CPair, T2CPair>();
    
    std::vector<R4CDist> rec_dists;
    
    for (const auto& component : components)
    {
        rec_dists.push_back(_get_vrr_recursion(component));
    }
    
    for (const auto& label : _get_buffers_str(rec_dists, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    if ((integral[1] == 0) || (integral[3] == 0))
    {
        const std::array<int, 2> rec_range({0, static_cast<int>(components.size())});
        
        for (const auto& label : _get_buffers_str(integral, components, rec_range))
        {
            lines.push_back({1, 0, 2, label});
        }
        
        _add_recursion_loop(lines, integral, components, rec_range);
    }
    else
    {
        const auto bcomps = t2c::number_of_cartesian_components(integral[1]);
        
        const auto kcomps = t2c::number_of_cartesian_components(integral[3]);
        
        for (int i = 0; i < bcomps; i++)
        {
            const std::array<int, 2> rec_range({i * kcomps, (i + 1) * kcomps});
            
            for (const auto& label : _get_buffers_str(integral, components, rec_range))
            {
                lines.push_back({1, 0, 2, label});
            }
            
            _add_recursion_loop(lines, integral, components, {i * kcomps, (i + 1) * kcomps});
            
            if (i < (bcomps - 1))  lines.push_back({0, 0, 1, ""});;
        }
    }
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CPrimFuncBodyDriver::_get_buffers_str(const std::vector<R4CDist>& rec_dists,
                                        const I4CIntegral&          integral) const
{
    std::vector<std::string> vstr;
    
    for (const auto& tint : t4c::get_vrr_integrals(integral))
    {
        auto label = "pbuffer.data(" + t4c::get_index_label(tint);
        
        vstr.push_back("/// Set up components of auxilary buffer : " + tint.label());
        
        const auto tlabel = _get_tensor_label(tint);
        
        int index = 0;
        
        for (const auto& tcomp : tint.components<T2CPair, T2CPair>())
        {
            if (_find_integral(rec_dists, tcomp))
            {
                const auto line = "auto " + _get_component_label(tcomp) + " = " + label;
                
                if (index > 0)
                {
                    vstr.push_back(line + " + " + std::to_string(index) + ");");
                }
                else
                {
                    vstr.push_back(line + ");");
                }
            }
            
            index++;
        }
    }
    
    return vstr;
}

bool
T4CPrimFuncBodyDriver::_find_integral(const std::vector<R4CDist>& rec_dists,
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

std::vector<std::string>
T4CPrimFuncBodyDriver::_get_buffers_str(const I4CIntegral&        integral,
                                        const VT4CIntegrals&      components,
                                        const std::array<int, 2>& rec_range) const
{
    std::vector<std::string> vstr;
    
    if ((rec_range[1] - rec_range[0]) == static_cast<int>(components.size()))
    {
        vstr.push_back("/// Set up components of targeted buffer : " + integral.label());
    }
    else
    {
        vstr.push_back("/// Set up " + std::to_string(rec_range[0]) + "-" + std::to_string(rec_range[1]) +
                       " components of targeted buffer : " + integral.label());
    }
    
    auto label = "pbuffer.data(" + t4c::get_index_label(integral);
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        const auto line = "auto " + _get_component_label(components[i]) + " = " + label;
        
        if (i > 0)
        {
            vstr.push_back(line + " + " + std::to_string(i) + ");");
        }
        else
        {
            vstr.push_back(line + ");");
        }
    }
    
    return vstr;
}

std::string
T4CPrimFuncBodyDriver::_get_tensor_label(const I4CIntegral& integral) const
{
    std::string label;
    
    if (integral.integrand().name() == "1/|r-r'|") label = "g";
    
    return label;
}

std::string
T4CPrimFuncBodyDriver::_get_tensor_label(const T4CIntegral& integral) const
{
    std::string label;
    
    if (integral.integrand().name() == "1/|r-r'|") label = "g";
    
    return label;
}

void
T4CPrimFuncBodyDriver::_add_recursion_loop(      VCodeLines&         lines,
                                           const I4CIntegral&        integral,
                                           const VT4CIntegrals&      components,
                                           const std::array<int, 2>& rec_range) const
{
    std::vector<R4CDist> rec_dists;
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        rec_dists.push_back(_get_vrr_recursion(components[i]));
    }
    
    // set up recursion loop
    
    const auto var_str = _get_pragma_str(integral, rec_dists);
    
    lines.push_back({1, 0, 1, "#pragma omp simd aligned(" + var_str + " : 64)"});
    
    lines.push_back({1, 0, 1, "for (size_t i = 0; i < nelems; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    _get_factor_lines(lines, rec_dists);
    
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
T4CPrimFuncBodyDriver::_get_pragma_str(const I4CIntegral&          integral,
                                       const std::vector<R4CDist>& rec_distributions) const
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
                if ((fact.order() > 0) &&
                    (fact.label() != "pb_x") &&
                    (fact.label() != "pb_y") &&
                    (fact.label() != "pb_z")) tlabels.insert(fact.label());
            }
        }
    }
    
    std::string label;
    
    for (const auto& tlabel : tlabels)
    {
        label += tlabel + ", ";
    }
    
    if ((integral[0] + integral[1]) > 1) label += "c_exps, d_exps ";
    
    if (label[label.size() - 2] == ',') label.erase(label.end() - 2);
    
    return label;
}

void
T4CPrimFuncBodyDriver::_get_factor_lines(                VCodeLines& lines,
                                         const std::vector<R4CDist>& rec_distributions) const
{
    std::set<std::string> tlabels;
    
    for (const auto& rdist : rec_distributions)
    {
        auto tint = rdist.root().integral();
        
        tlabels.insert(_get_tensor_label(tint) + "_" + tint.label());
        
        for (size_t i = 0; i < rdist.terms(); i++)
        {
            for (const auto& fact : rdist[i].factors())
            {
                if (fact.order() == 0) tlabels.insert(fact.label());
            }
        }
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "fi_ab_0") !=  tlabels.end())
    {
        lines.push_back({2, 0, 2, "const double fi_ab_0 = 0.5 / (a_exp + b_exp);"});
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "fi_cd_0") !=  tlabels.end())
    {
        lines.push_back({2, 0, 2, "const double fi_cd_0 = 0.5 / (c_exps[i] + d_exps[i]);"});
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "fi_abcd_0") !=  tlabels.end())
    {
        lines.push_back({2, 0, 2, "const double fi_abcd_0 = 0.5 / (a_exp + b_exp + c_exps[i] + d_exps[i]);"});
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "fti_ab_0") !=  tlabels.end())
    {
        if (std::find(tlabels.begin(), tlabels.end(), "fi_abcd_0") !=  tlabels.end())
        {
            lines.push_back({2, 0, 2, "const double fti_ab_0 = 2.0 * fi_abcd_0 * fi_ab_0 * (c_exps[i] + d_exps[i]);"});
        }
        else
        {
            lines.push_back({2, 0, 2, "const double fti_ab_0 =  fi_ab_0 * (c_exps[i] + d_exps[i]) / (a_exp + b_exp + c_exps[i] + d_exps[i]);"});
        }
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "fti_cd_0") !=  tlabels.end())
    {
        if (std::find(tlabels.begin(), tlabels.end(), "fi_abcd_0") !=  tlabels.end())
        {
            lines.push_back({2, 0, 2, "const double fti_cd_0 = 2.0 * fi_abcd_0 * fi_cd_0 * (a_exp + b_exp);"});
        }
        else
        {
            lines.push_back({2, 0, 2, "const double fti_cd_0 =  fi_cd_0 * (a_exp + b_exp) / (a_exp + b_exp + c_exps[i] + d_exps[i]);"});
        }
    }
}

R4CDist
T4CPrimFuncBodyDriver::_get_vrr_recursion(const T4CIntegral& integral) const
{
    R4CDist rdist;
    
    if (integral.integrand().name() == "1/|r-r'|")
    {
        T4CVrrElectronRepulsionDriver eri_drv;

        if (integral[1].order() > 0)
        {
            rdist = eri_drv.apply_bra_vrr(R4CTerm(integral));
        }
        else
        {
            rdist = eri_drv.apply_ket_vrr(R4CTerm(integral));
        }
    }
    
    rdist.simplify();
    
    return rdist;
}

std::string
T4CPrimFuncBodyDriver::_get_code_line(const R4CDist& rec_distribution) const
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
T4CPrimFuncBodyDriver::_get_rterm_code(const R4CTerm& rec_term,
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
        plabel+= " * " + fact.label();
            
        if ((fact.order() > 0)       &&
            (fact.label() != "pb_x") &&
            (fact.label() != "pb_y") &&
            (fact.label() != "pb_z")) plabel += "[i]";
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
T4CPrimFuncBodyDriver::_get_component_label(const T4CIntegral& integral) const
{
    std::string label = _get_tensor_label(integral) + "_" + integral.label();
    
    label += "_" + std::to_string(integral.order());
    
    return label;
}
