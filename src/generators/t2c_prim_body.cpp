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

#include "t2c_prim_body.hpp"

#include <algorithm>
#include <iostream>

#include "t2c_utils.hpp"
#include "t2c_ovl_driver.hpp"
#include "t2c_kin_driver.hpp"
#include "t2c_center_driver.hpp"
#include "t2c_npot_driver.hpp"
#include "t2c_dip_driver.hpp"
#include "t2c_linmom_driver.hpp"
#include "t2c_el_field_driver.hpp"
#include "t2c_eri_driver.hpp"
#include "t3c_ovl_driver.hpp"
#include "t3c_ovl_grad_driver.hpp"
#include "t3c_r2_driver.hpp"
#include "t3c_rr2_driver.hpp"

void
T2CPrimFuncBodyDriver::write_func_body(      std::ofstream& fstream,
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
        
    const auto kcomps = t2c::number_of_cartesian_components(integral[1]);
    
    if ((integral[0] == 0) || (integral[1] == 0))
    {
        const std::array<int, 2> rec_range({0, ncomps});

        for (const auto& label : _get_buffers_str(integral, components, rec_range))
        {

            lines.push_back({1, 0, 2, label});
        }
        
        _add_recursion_loop(lines, integral, components, rec_range);
    }
    else
    {
        const auto nblocks = ncomps / kcomps;
        
        for (int i = 0; i < nblocks; i++)
        {
            const std::array<int, 2> rec_range({i * kcomps, (i + 1) * kcomps});
            
            for (const auto& label : _get_buffers_str(integral, components, rec_range))
            {
                lines.push_back({1, 0, 2, label});
            }
            
            _add_recursion_loop(lines, integral, components, {i * kcomps, (i + 1) * kcomps});
            
            if (i < (ncomps - 1))  lines.push_back({0, 0, 1, ""});;
        }
    }
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CPrimFuncBodyDriver::_get_factors_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (_need_exponents(integral))
    {
        vstr.push_back("// Set up exponents");

        vstr.push_back("auto b_exps = factors.data(0);");
    }
    
    if (_need_distances_pa(integral))
    {
        if (integral.integrand().name() == "G(r)")
        {
            vstr.push_back("// Set up R(GA) distances");

            vstr.push_back("auto ga_x = factors.data(idx_rga);");
            
            vstr.push_back("auto ga_y = factors.data(idx_rga + 1);");

            vstr.push_back("auto ga_z = factors.data(idx_rga + 2);");
        }
        else
        {
            vstr.push_back("// Set up R(PA) distances");

            vstr.push_back("auto pa_x = factors.data(idx_rpa);");

            vstr.push_back("auto pa_y = factors.data(idx_rpa + 1);");

            vstr.push_back("auto pa_z = factors.data(idx_rpa + 2);");
        }
    }
    
    if (_need_distances_pb(integral))
    {
        if (integral.integrand().name() == "G(r)")
        {
            vstr.push_back("// Set up R(GB) distances");

            vstr.push_back("auto gb_x = factors.data(idx_rgb);");

            vstr.push_back("auto gb_y = factors.data(idx_rgb + 1);");

            vstr.push_back("auto gb_z = factors.data(idx_rgb + 2);");
        }
        else
        {
            vstr.push_back("// Set up R(PB) distances");

            vstr.push_back("auto pb_x = factors.data(idx_rpb);");

            vstr.push_back("auto pb_y = factors.data(idx_rpb + 1);");

            vstr.push_back("auto pb_z = factors.data(idx_rpb + 2);");
        }
    }
    
    if (_need_distances_pc(integral))
    {
        vstr.push_back("// Set up R(PC) distances");

        vstr.push_back("auto pc_x = factors.data(idx_rpc);");

        vstr.push_back("auto pc_y = factors.data(idx_rpc + 1);");

        vstr.push_back("auto pc_z = factors.data(idx_rpc + 2);");
    }
    
    if (
        (integral.integrand().name() == "GX(r)")  ||
        (integral.integrand().name() == "GR2(r)") ||
        (integral.integrand().name() == "GR.R2(r)")
        )
    {
        vstr.push_back("// Set up R(GC) distances");

        vstr.push_back("auto gc_x = factors.data(idx_rgc);");

        vstr.push_back("auto gc_y = factors.data(idx_rgc + 1);");

        vstr.push_back("auto gc_z = factors.data(idx_rgc + 2);");
    }
    
    return vstr;
}

std::vector<std::string>
T2CPrimFuncBodyDriver::_get_buffers_str(const std::vector<R2CDist>& rec_dists,
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

bool
T2CPrimFuncBodyDriver::_find_integral(const std::vector<R2CDist>& rec_dists,
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

std::vector<std::string>
T2CPrimFuncBodyDriver::_get_buffers_str(const I2CIntegral&        integral,
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

// MR: Change for new integral cases
std::string
T2CPrimFuncBodyDriver::_get_tensor_label(const I2CIntegral& integral) const
{
    std::string label;
    
    if (integral.integrand().name() == "1") label = "ts";
    
    if (integral.integrand().name() == "T") label = "tk";
    
    if (integral.integrand().name() == "A") label = "ta";

    if (integral.integrand().name() == "r") label = "tr";

    if (integral.integrand().name() == "p") label = "tp";

    if (integral.integrand().name() == "AG") label = "ta" + std::to_string(integral.integrand().shape().order());
    
    if (integral.integrand().name() == "1/|r-r'|") label = "g";

    return label;
}

// MR: Change for new integral cases
std::string
T2CPrimFuncBodyDriver::_get_tensor_label(const T2CIntegral& integral) const
{
    std::string label;
    
    if (integral.integrand().name() == "1") label = "ts";
    
    if (integral.integrand().name() == "T") label = "tk";
    
    if (integral.integrand().name() == "A") label = "ta";

    if (integral.integrand().name() == "r") label = "tr";

    if (integral.integrand().name() == "p") label = "tp";

    if (integral.integrand().name() == "AG") label = "ta" + std::to_string(integral.integrand().shape().order());
    
    if (integral.integrand().name() == "1/|r-r'|") label = "g";
    
    if (integral.integrand().name() == "G(r)") label = "ts";
    
    if (integral.integrand().name() == "GX(r)") label = "gs";
    
    if (integral.integrand().name() == "GR2(r)") label = "gr";
    
    if (integral.integrand().name() == "GR.R2(r)") label = "grr";

    return label;
}

void
T2CPrimFuncBodyDriver::_add_recursion_loop(      VCodeLines&         lines,
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
T2CPrimFuncBodyDriver::_get_pragma_str(const I2CIntegral&          integral,
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
    
    if (
        (integral.integrand().name() == "GR2(r)") ||
        (integral.integrand().name() == "GR.R2(r)")
        )
    {
        tlabels.insert("gc_x");
        
        tlabels.insert("gc_y");
        
        tlabels.insert("gc_z");
    }
    
    std::string label;
    
    for (const auto& tlabel : tlabels)
    {
        label += tlabel + ", ";
    }
    
    if (_need_exponents(integral)) label += "b_exps";
    
    if (label[label.size() - 2] == ',') label.erase(label.end() - 2);
    
    return label;
}

// MR: Possibly change for new integral cases
void
T2CPrimFuncBodyDriver::_get_factor_lines(                VCodeLines& lines,
                                         const std::vector<R2CDist>& rec_distributions) const
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
    
    if (std::find(tlabels.begin(), tlabels.end(), "fe_0") !=  tlabels.end())
    {
        lines.push_back({2, 0, 2, "const double fe_0 = 0.5 / (a_exp + b_exps[i]);"});
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "fz_0") !=  tlabels.end())
    {
        if (std::find(tlabels.begin(), tlabels.end(), "fe_0") !=  tlabels.end())
        {
            lines.push_back({2, 0, 2, "const double fz_0 = 2.0 * a_exp * b_exps[i] * fe_0;"});
        }
        else
        {
            lines.push_back({2, 0, 2, "const double fz_0 = a_exp * b_exps[i] / (a_exp + b_exps[i]);"});
        }
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "tbe_0") !=  tlabels.end())
    {
        lines.push_back({2, 0, 2, "const double tbe_0 = a_exp;"});
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "tce_0") !=  tlabels.end())
    {
        lines.push_back({2, 0, 2, "const double tce_0 = c_exp;"});
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "rgc2_0") !=  tlabels.end())
    {
        lines.push_back({2, 0, 2, "const double rgc2_0 = gc_x[i] * gc_x[i] + gc_y[i] * gc_y[i] + gc_z[i] * gc_z[i];"});
    }

    if (std::find(tlabels.begin(), tlabels.end(), "fbe_0") !=  tlabels.end())
    {
        lines.push_back({2, 0, 2, "const double fbe_0 = 0.5 / a_exp;"});
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "fke_0") !=  tlabels.end())
    {
        lines.push_back({2, 0, 2, "const double fke_0 = 0.5 / b_exps[i];"});
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "fz_be_0") !=  tlabels.end())
    {
        if (std::find(tlabels.begin(), tlabels.end(), "fe_0") !=  tlabels.end())
        {
            lines.push_back({2, 0, 2, "const double fz_be_0 =  2.0 * b_exps[i] * fe_0 * fbe_0;"});
        }
        else
        {
            lines.push_back({2, 0, 2, "const double fz_be_0 = b_exps[i] * fbe_0 / (a_exp + b_exps[i]);"});
        }
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "fz_ke_0") !=  tlabels.end())
    {
        if (std::find(tlabels.begin(), tlabels.end(), "fe_0") !=  tlabels.end())
        {
            lines.push_back({2, 0, 2, "const double fz_ke_0 =  2.0 * a_exp * fe_0 * fke_0;"});
        }
        else
        {
            lines.push_back({2, 0, 2, "const double fz_ke_0 = a_exp * fke_0 / (a_exp + b_exps[i]);"});
        }
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "gfe_0") !=  tlabels.end())
    {
        lines.push_back({2, 0, 2, "const double gfe_0 = 0.5 / (a_exp + b_exps[i] + c_exp);"});
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "gfe2_0") !=  tlabels.end())
    {
        lines.push_back({2, 0, 2, "const double gfe2_0 = gfe_0 * gfe_0;"});
    }
}

R2CDist
T2CPrimFuncBodyDriver::_get_vrr_recursion(const T2CIntegral& integral) const
{
    R2CDist rdist;

    if (!integral.prefixes().empty())
    {
        T2CCenterDriver geom_drv;

        const auto prefixes = integral.prefixes();

        if ((prefixes.size() == 2) || (prefixes.size() == 1))
        {

            if ((integral.prefixes()[0].shape().order() == 0) &&
                (integral.prefixes()[1].shape().order() >  0))
            {
                return geom_drv.apply_bra_ket_vrr(integral, 1);
            }
            else
            {
                return geom_drv.apply_bra_ket_vrr(integral, 0);
            }
        }
    }

    if (integral.integrand().name() == "1")
    {
        T2COverlapDriver ovl_drv;
        
        if (integral[0].order() > 0)
        {
            rdist = ovl_drv.apply_bra_vrr(R2CTerm(integral));
        }
        else
        {
            rdist = ovl_drv.apply_ket_vrr(R2CTerm(integral));
        }
    }
    
    if (integral.integrand().name() == "T")
    {
        T2CKineticEnergyDriver kin_drv;
        
        if (integral[0].order() > 0)
        {
            rdist = kin_drv.apply_bra_vrr(R2CTerm(integral));
        }
        else
        {
            rdist = kin_drv.apply_ket_vrr(R2CTerm(integral));
        }
    }


    if (integral.integrand().name() == "r")
    {
        T2CMultipoleDriver dip_drv;


        if (integral[0].order() > 0)
        {

            rdist = dip_drv.apply_bra_vrr(R2CTerm(integral));
        }
        else
        {
            rdist = dip_drv.apply_ket_vrr(R2CTerm(integral));
        }
    }

    if (integral.integrand().name() == "p")
    {
        T2CLinearMomentumDriver linmom_drv;

        rdist = linmom_drv.apply_op_vrr(R2CTerm(integral));
    }


    if (integral.integrand().name() == "A")
    {
        T2CNuclearPotentialDriver npot_drv;
        
        if (integral[0].order() > 0)
        {
            rdist = npot_drv.apply_bra_vrr(R2CTerm(integral));
        }
        else
        {
            rdist = npot_drv.apply_ket_vrr(R2CTerm(integral));
        }
    }

    if (integral.integrand().name() == "AG")
    {
        T2CElectricFieldDriver el_field_drv;

        if (integral[0].order() > 0)
        {
            rdist = el_field_drv.apply_bra_vrr(R2CTerm(integral));
        }
        else
        {
            rdist = el_field_drv.apply_ket_vrr(R2CTerm(integral));
        }
    }
    
    if (integral.integrand().name() == "G(r)")
    {
        T3COverlapDriver ovl_drv;

        if (integral[0].order() > 0)
        {
            rdist = ovl_drv.apply_bra_vrr(R2CTerm(integral));
        }
        else
        {
            rdist = ovl_drv.apply_ket_vrr(R2CTerm(integral));
        }
    }
    
    if (integral.integrand().name() == "GX(r)")
    {
        T3COverlapGradientDriver ovl_grad_drv;

        rdist = ovl_grad_drv.apply_aux_vrr(R2CTerm(integral));
    }
    
    if (integral.integrand().name() == "GR2(r)")
    {
        T3CR2Driver r2_drv;

        rdist = r2_drv.aux_vrr(R2CTerm(integral));
    }
    
    if (integral.integrand().name() == "GR.R2(r)")
    {
        T3CRR2Driver rr2_drv;

        rdist = rr2_drv.apply_aux_vrr(R2CTerm(integral));
    }
    
    if (integral.integrand().name() == "1/|r-r'|")
    {
        T2CElectronRepulsionDriver eri_drv;
        
        if (integral[0].order() > 0)
        {
            rdist = eri_drv.apply_bra_vrr(R2CTerm(integral));
        }
        else
        {
            rdist = eri_drv.apply_ket_vrr(R2CTerm(integral));
        }
    }
    
    rdist.simplify();
    
    return rdist;
}

std::string
T2CPrimFuncBodyDriver::_get_code_line(const R2CDist& rec_distribution) const
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
T2CPrimFuncBodyDriver::_get_rterm_code(const R2CTerm& rec_term,
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

// MR: May need to change for new integral cases
std::string
T2CPrimFuncBodyDriver::_get_component_label(const T2CIntegral& integral) const
{
    std::string label = _get_tensor_label(integral) + "_" + integral.label();
    
    if (integral.integrand().name() == "A")
    {
        label += "_" + std::to_string(integral.order());
    }

    if (integral.integrand().name() == "AG")
    {
        label += "_" + std::to_string(integral.order());
    }
    
    if (integral.integrand().name() == "1/|r-r'|")
    {
        label += "_" + std::to_string(integral.order());
    }
    
    return label;
}

bool
T2CPrimFuncBodyDriver::_need_distances_pa(const I2CIntegral& integral) const
{
    if (integral.integrand().name() == "GX(r)") return false;
    
    if (integral.integrand().name() == "GR2(r)") return false;
    
    if (integral.integrand().name() == "GR.R2(r)") return false;
    
    return integral[0] > 0;
}

bool
T2CPrimFuncBodyDriver::_need_distances_pb(const I2CIntegral& integral) const
{
    if (integral.integrand().name() == "GX(r)") return false;
    
    if (integral.integrand().name() == "GR2(r)") return false;
    
    if (integral.integrand().name() == "GR.R2(r)") return false;
    
    return (integral[0] == 0) && (integral[1] > 0);
}

bool
T2CPrimFuncBodyDriver::_need_distances_pc(const I2CIntegral& integral) const
{    
    if (integral.integrand().name() == "A") return true;
    
    if (integral.integrand().name() == "AG") return true;
    
    return false;
}

bool
T2CPrimFuncBodyDriver::_need_exponents(const I2CIntegral& integral) const
{
    if (integral.integrand().name() == "T") return true;
    
    if (integral.integrand().name() == "r") return true;
    
    if (integral.integrand().name() == "GX(r)") return true;
    
    if (integral.integrand().name() == "GR2(r)") return true;
    
    if (integral.integrand().name() == "GR.R2(r)") return true;
    
    return (integral[0] + integral[1]) > 1;
}


