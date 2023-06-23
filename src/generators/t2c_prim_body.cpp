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

#include "file_stream.hpp"
#include "t2c_utils.hpp"
#include "t2c_ovl_driver.hpp"
#include "t2c_kin_driver.hpp"

#include <iostream>

void
T2CPrimFuncBodyDriver::write_prim_func_body(      std::ofstream& fstream,
                                            const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_common_data_str())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_buffers_str(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_func_pragma(lines, integral);
    
    _add_loop_start(lines, integral);
    
     const auto tcomps = integral.components<T1CPair, T1CPair>();
   
    std::vector<std::string> labels({"fints", });
    
    if (integral[0] > 0) labels = t2c::tensor_components(Tensor(integral[0]), "fints");
    
    if (integral[1] > 0) labels = t2c::tensor_components(Tensor(integral[1]), "fints");
    
    _add_simd_code(lines, labels, tcomps, integral);
    
    _add_loop_end(lines);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CPrimFuncBodyDriver::write_prim_func_body(      std::ofstream&   fstream,
                                            const TensorComponent& component,
                                            const I2CIntegral&     integral,
                                            const bool             bra_first) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_common_data_str())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_buffers_str(component, integral, bra_first))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_func_pragma(lines, component, integral, bra_first);
    
    _add_loop_start(lines, integral);
    
    const auto tcomps = _select_integral_components(component, integral, bra_first);
    
    const auto labels = (bra_first) ? t2c::tensor_components(Tensor(integral[1]), "fints")
                                    : t2c::tensor_components(Tensor(integral[0]), "fints");
    
    _add_simd_code(lines, labels, tcomps, integral);
    
    _add_loop_end(lines);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CPrimFuncBodyDriver::write_prim_func_body(      std::ofstream&   fstream,
                                            const TensorComponent& bra_component,
                                            const TensorComponent& ket_component,
                                            const I2CIntegral&     integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_common_data_str())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_buffers_str(bra_component, ket_component, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_func_pragma(lines, bra_component, ket_component, integral);
    
    _add_loop_start(lines, integral);
    
    const auto tcomps = _select_integral_components(bra_component, ket_component, integral);
    
    const auto labels = t2c::integrand_components(integral.integrand(), "fints");
    
    _add_simd_code(lines, labels, tcomps, integral);
    
    _add_loop_end(lines);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CPrimFuncBodyDriver::_get_common_data_str() const
{
    std::vector<std::string> vstr;
        
    vstr.push_back("// set up math constants");
        
    vstr.push_back("const auto fpi = mathconst::getPiValue();");
        
    vstr.push_back("// set up coordinates for bra side");
        
    vstr.push_back("const auto bra_rx = bra_coord[0];");
        
    vstr.push_back("const auto bra_ry = bra_coord[1];");
        
    vstr.push_back("const auto bra_rz = bra_coord[2];");
        
    vstr.push_back("// set up coordinates for ket side");
        
    vstr.push_back("auto ket_rx = ket_coords_x.data();");
        
    vstr.push_back("auto ket_ry = ket_coords_y.data();");
        
    vstr.push_back("auto ket_rz = ket_coords_z.data();");
        
    vstr.push_back("// set exponents and normalization factors on ket side");
        
    vstr.push_back("auto ket_fe = ket_exps.data();");
        
    vstr.push_back("auto ket_fn = ket_norms.data();");
        
    return vstr;
}

std::vector<std::string>
T2CPrimFuncBodyDriver::_get_buffers_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// set up pointer to integrals buffer(s)");
    
    std::vector<std::string> labels({"buffer", });
    
    std::vector<std::string> flabels({"fints", });
    
    if (integral[0] > 0)
    {
        const auto bra = Tensor(integral[0]);
        
        labels = t2c::tensor_components(bra, "buffer");
        
        flabels = t2c::tensor_components(bra, "fints");
    }
       
    if (integral[1] > 0)
    {
        const auto ket = Tensor(integral[1]);
        
        labels = t2c::tensor_components(ket, "buffer");
        
        flabels = t2c::tensor_components(ket, "fints");
    }
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        vstr.push_back("auto " + flabels[i] + " = " + labels[i] + ".data();");
    }
    
    return vstr;
}

std::vector<std::string>
T2CPrimFuncBodyDriver::_get_buffers_str(const TensorComponent& component,
                                        const I2CIntegral&     integral,
                                        const bool             bra_first) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// set up pointer to integrals buffer(s)");
    
    std::vector<std::string> labels, flabels;
    
    if (bra_first)
    {
        const auto ket = Tensor(integral[1]);
        
        labels = t2c::tensor_components(ket, "buffer");
        
        flabels = t2c::tensor_components(ket, "fints");
    }
    else
    {
        const auto bra = Tensor(integral[0]);
        
        labels = t2c::tensor_components(bra, "buffer");
        
        flabels = t2c::tensor_components(bra, "fints");
    }
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        vstr.push_back("auto " + flabels[i] + " = " + labels[i] + ".data();");
    }
    
    return vstr;
}

std::vector<std::string>
T2CPrimFuncBodyDriver::_get_buffers_str(const TensorComponent& bra_component,
                                        const TensorComponent& ket_component,
                                        const I2CIntegral&     integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// set up pointer to integrals buffer(s)");
    
    const auto labels  = t2c::integrand_components(integral.integrand(), "buffer");
    
    const auto flabels = t2c::integrand_components(integral.integrand(), "fints");
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        vstr.push_back("auto " + flabels[i] + " = " + labels[i] + ".data();");
    }
    
    return vstr;
}

void
T2CPrimFuncBodyDriver::_add_func_pragma(      VCodeLines&  lines,
                                        const I2CIntegral& integral) const
{
    std::vector<std::string> labels({"fints", });
    
    if (integral[0] > 0) labels = t2c::tensor_components(Tensor(integral[0]), "fints");
    
    if (integral[1] > 0) labels = t2c::tensor_components(Tensor(integral[1]), "fints");
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        if (i == 0)
        {
            lines.push_back({1, 0, 1, "#pragma omp simd aligned(" + labels[i] + ",\\"});
        }
        else
        {
            lines.push_back({1, 25, 1, labels[i] + ",\\"});
        }
    }

    _add_common_pragma(lines);
}

void
T2CPrimFuncBodyDriver::_add_func_pragma(      VCodeLines&      lines,
                                        const TensorComponent& component,
                                        const I2CIntegral&     integral,
                                        const bool             bra_first) const
{
    const auto labels = (bra_first) ? t2c::tensor_components(Tensor(integral[1]), "fints")
                                    : t2c::tensor_components(Tensor(integral[0]), "fints");
        
    for (size_t i = 0; i < labels.size(); i++)
    {
        if (i == 0)
        {
            lines.push_back({1, 0, 1, "#pragma omp simd aligned(" + labels[i] + ",\\"});
        }
        else
        {
            lines.push_back({1, 25, 1, labels[i] + ",\\"});
        }
    }
    
    _add_common_pragma(lines);
}

void
T2CPrimFuncBodyDriver::_add_func_pragma(      VCodeLines&      lines,
                                        const TensorComponent& bra_component,
                                        const TensorComponent& ket_component,
                                        const I2CIntegral&     integral) const
{
    const auto labels = t2c::integrand_components(integral.integrand(), "fints");
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        if (i == 0)
        {
            lines.push_back({1, 0, 1, "#pragma omp simd aligned(" + labels[i] + ",\\"});
        }
        else
        {
            lines.push_back({1, 25, 1, labels[i] + ",\\"});
        }
    }
    
    _add_common_pragma(lines);
}

void
T2CPrimFuncBodyDriver::_add_common_pragma(VCodeLines& lines) const
{
    lines.push_back({1, 25, 1, "ket_fe,\\"});
        
    lines.push_back({1, 25, 1, "ket_fn,\\"});
        
    lines.push_back({1, 25, 1, "ket_rx,\\"});
        
    lines.push_back({1, 25, 1, "ket_ry,\\"});
        
    lines.push_back({1, 25, 1, "ket_rz : 64)"});
}

void
T2CPrimFuncBodyDriver::_add_loop_start(      VCodeLines&  lines,
                                       const I2CIntegral& integral) const
{
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < ket_dim; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "const auto ab_x = bra_rx - ket_rx[i];"});
    
    lines.push_back({2, 0, 2, "const auto ab_y = bra_ry - ket_ry[i];"});
    
    lines.push_back({2, 0, 2, "const auto ab_z = bra_rz - ket_rz[i];"});
    
    if (integral.integrand() == Operator("1"))
    {
        _add_overlap_vars(lines, integral);
    }
    
    if (integral.integrand() == Operator("T"))
    {
        _add_kinetic_energy_vars(lines, integral);
    }
}

void
T2CPrimFuncBodyDriver::_add_overlap_vars(      VCodeLines&  lines,
                                         const I2CIntegral& integral) const
{
    lines.push_back({2, 0, 2, "const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);"});
    
    lines.push_back({2, 0, 2, "auto fz_0 = bra_exp * ket_fe[i] * fe_0;"});
    
    lines.push_back({2, 0, 2, "fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);"});
    
    if ((integral[0] + integral[1]) == 0)
    {
        lines.push_back({2, 0, 1, "fints[i] += bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);"});
    }
    else
    {
        lines.push_back({2, 0, 2, "const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);"});
    }
}

void
T2CPrimFuncBodyDriver::_add_kinetic_energy_vars(      VCodeLines&  lines,
                                                const I2CIntegral& integral) const
{
    lines.push_back({2, 0, 2, "const auto r2ab = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;"});
    
    lines.push_back({2, 0, 2, "const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);"});
    
    lines.push_back({2, 0, 2, "const auto fz_0 = bra_exp * ket_fe[i] * fe_0;"});
    
    lines.push_back({2, 0, 2, "const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0 * r2ab);"});
    
    if ((integral[0] + integral[1]) == 0)
    {
        lines.push_back({2, 0, 1, "fints[i] += fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;"});
    }
    else
    {
        lines.push_back({2, 0, 2, "const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;"});
    }
}

void
T2CPrimFuncBodyDriver::_add_loop_end(VCodeLines& lines) const
{
    lines.push_back({1, 0, 1, "}"});
}

VT2CIntegrals
T2CPrimFuncBodyDriver::_select_integral_components(const TensorComponent& component,
                                                   const I2CIntegral&     integral,
                                                   const bool             bra_first) const
{
    VT2CIntegrals tcomps;
    
    for (const auto& tcomp : integral.components<T1CPair, T1CPair>())
    {
        if (bra_first)
        {
            if (tcomp.bra().shape() == component) tcomps.push_back(tcomp);
        }
        else
        {
            if (tcomp.ket().shape() == component) tcomps.push_back(tcomp);
        }
    }
        
    return tcomps;
}

VT2CIntegrals
T2CPrimFuncBodyDriver::_select_integral_components(const TensorComponent& bra_component,
                                                   const TensorComponent& ket_component,
                                                   const I2CIntegral&     integral) const
{
    VT2CIntegrals tcomps;
    
    for (const auto& tcomp : integral.components<T1CPair, T1CPair>())
    {
        if ((tcomp.bra().shape() == bra_component) &&
            (tcomp.ket().shape() == ket_component))
        {
            tcomps.push_back(tcomp);
        }
    }
        
    return tcomps;
}

void
T2CPrimFuncBodyDriver::_add_simd_code(      VCodeLines&               lines,
                                      const std::vector<std::string>& labels,
                                      const VT2CIntegrals&            components,
                                      const I2CIntegral&              integral) const
{
    const auto rgroup = _generate_integral_group(components, integral);
    
    _add_prefactors(lines, rgroup);
        
    _add_simd_lines(lines, labels, rgroup);
}

R2Group
T2CPrimFuncBodyDriver::_generate_integral_group(const VT2CIntegrals& components,
                                                const I2CIntegral&   integral) const
{
    R2Group rgroup;
    
    // Overlap inntegrals
    
    if (integral.integrand() == Operator("1"))
    {
        T2COverlapDriver t2c_ovl_drv;
        
        rgroup = t2c_ovl_drv.create_recursion(components);
    }
    
    // Kinetic energy inntegrals
    
    if (integral.integrand() == Operator("T"))
    {
        T2CKineticEnergyDriver t2c_kin_drv;
        
        rgroup = t2c_kin_drv.create_recursion(components);
    }
    
    // ... other integrals
    
    return rgroup;
}

void
T2CPrimFuncBodyDriver::_add_prefactors(      VCodeLines& lines,
                                       const R2Group&    rgroup) const
{
    if (t2c::find_factor(rgroup, "rpa_x"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_x = -ket_fe[i] * ab_x * fe_0;"});
    }
    
    if (t2c::find_factor(rgroup, "rpa_y"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_y = -ket_fe[i] * ab_y * fe_0;"});
    }
    
    if (t2c::find_factor(rgroup, "rpa_z"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_z = -ket_fe[i] * ab_z * fe_0;"});
    }
    
    if (t2c::find_factor(rgroup, "rpb_x"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_x = bra_exp * ab_x * fe_0;"});
    }
    
    if (t2c::find_factor(rgroup, "rpb_y"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_y = bra_exp * ab_y * fe_0;"});
    }
    
    if (t2c::find_factor(rgroup, "rpb_z"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_z = bra_exp * ab_z * fe_0;"});
    }
    
    if (t2c::find_factor(rgroup, "fke_0"))
    {
        lines.push_back({2, 0, 2, "const auto fke_0 = 1.0 / ket_fe[i];"});
    }
    
    if (t2c::find_factor(rgroup, "fbe_0"))
    {
        lines.push_back({2, 0, 2, "const auto fbe_0 = 1.0 / bra_exp;"});
    }
}

void
T2CPrimFuncBodyDriver::_add_simd_lines(      VCodeLines&               lines,
                                       const std::vector<std::string>& labels,
                                       const R2Group&                  rgroup) const
{
    if (rgroup.empty()) return;
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        const auto rdist = rgroup[i];
        
        for (const auto& tint : rdist.unique_integrals())
        {
            const auto tdist = rdist.split(tint);
            
            _add_simd_lines_block(lines, labels[i], tint, tdist);
        }
    }
}

void
T2CPrimFuncBodyDriver::_add_simd_lines_block(      VCodeLines&  lines,
                                             const std::string& label,
                                             const T2CIntegral& integral,
                                             const R2CDist&     rdist) const
{
    const auto tlabel = _get_aux_label(integral);
    
    const auto nterms = rdist.terms();

    auto nbatches = nterms / 5;

    if ((nterms % 5) != 0) nbatches++;

    for (size_t i = 0; i < nbatches; i++)
    {
        const auto sterm = 5 * i;
    
        const auto eterm = ((sterm + 5) > nterms) ? nterms : sterm + 5;
    
        std::string simd_str;
    
        for (size_t j = sterm; j < eterm; j++)
        {
            simd_str += t2c::get_factor_label(rdist[j], j == sterm);
        }
    
        if ((eterm - sterm) > 1) simd_str = "(" + simd_str + ")";

        lines.push_back({2, 0, 2, label + "[i] += " + tlabel + " * " + simd_str + ";"});
    }
}

std::string
T2CPrimFuncBodyDriver::_get_aux_label(const T2CIntegral& integral) const
{
    if (integral.integrand() == OperatorComponent("1"))
    {
        return std::string("fss");
    }
    
    if (integral.integrand() == OperatorComponent("T"))
    {
        return std::string("ftt");
    }
    
    return std::string();
}
