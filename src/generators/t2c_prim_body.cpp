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
#include "t2c_npot_driver.hpp"
#include "t2c_npot_geom_driver.hpp"
#include "t2c_mpol_driver.hpp"
#include "t2c_center_driver.hpp"

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
    
    for (const auto& label : _get_special_vars_str(integral, true))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_buffers_str(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_boys_vars_str(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_boys_compute_lines(lines, integral);
    
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
    
    for (const auto& label : _get_special_vars_str(integral, true))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_buffers_str(component, integral, bra_first))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_boys_vars_str(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_boys_compute_lines(lines, integral);
    
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
    
    for (const auto& label : _get_special_vars_str(integral, true))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_buffers_str(bra_component, ket_component, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_boys_vars_str(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_boys_compute_lines(lines, integral);
    
    _add_func_pragma(lines, bra_component, ket_component, integral);
    
    _add_loop_start(lines, integral);
    
    const auto tcomps = _select_integral_components(bra_component, ket_component, integral);
    
    std::vector<std::string> labels;
    
    const auto prefixes = integral.prefixes();
    
    if (prefixes.empty())
    {
        labels = t2c::integrand_components(integral.integrand(), "fints");
    }
    
    if (prefixes.size() == 1)
    {
        labels = t2c::integrand_components(prefixes[0].shape(), integral.integrand(), "fints");
    }
    
    if (prefixes.size() == 2)
    {
        labels = t2c::integrand_components(prefixes[0].shape(),prefixes[1].shape(), integral.integrand(), "fints");
    }
    
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
    
    std::vector<std::string> labels, flabels;
    
    const auto prefixes = integral.prefixes();
    
    if (prefixes.empty())
    {
        labels  = t2c::integrand_components(integral.integrand(), "buffer");
        
        flabels = t2c::integrand_components(integral.integrand(), "fints");
    }
    
    if (prefixes.size() == 1)
    {
        labels  = t2c::integrand_components(prefixes[0].shape(), integral.integrand(), "buffer");
        
        flabels = t2c::integrand_components(prefixes[0].shape(), integral.integrand(), "fints");
    }
    
    if (prefixes.size() == 2)
    {
        labels = t2c::integrand_components(prefixes[0].shape(), prefixes[1].shape(), integral.integrand(), "buffer");
        
        flabels = t2c::integrand_components(prefixes[0].shape(),prefixes[1].shape(), integral.integrand(), "fints");
    }
    
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
    std::vector<std::string> labels;
    
    const auto prefixes = integral.prefixes();
    
    if (prefixes.empty())
    {
        labels = t2c::integrand_components(integral.integrand(), "fints");
    }
    
    if (prefixes.size() == 1)
    {
        labels = t2c::integrand_components(prefixes[0].shape(), integral.integrand(), "fints");
    }
    
    if (prefixes.size() == 2)
    {
        labels = t2c::integrand_components(prefixes[0].shape(),prefixes[1].shape(), integral.integrand(), "fints");
    }
    
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
    
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "1")
    {
        _add_overlap_vars(lines, integral);
    }
    
    if (integrand.name() == "T")
    {
        _add_kinetic_energy_vars(lines, integral);
    }
    
    if (integrand.name() == "A")
    {
        _add_nuclear_potential_vars(lines, integral);
    }
    
    if (integrand.name() == "AG")
    {
        _add_nuclear_potential_geom_vars(lines, integral);
    }
    
    if (integrand.name() == "r")
    {
        _add_multipole_vars(lines, integral);
    }
}

void
T2CPrimFuncBodyDriver::_add_overlap_vars(      VCodeLines&  lines,
                                         const I2CIntegral& integral) const
{
    lines.push_back({2, 0, 2, "const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);"});
    
    lines.push_back({2, 0, 2, "auto fz_0 = bra_exp * ket_fe[i] * fe_0;"});
    
    lines.push_back({2, 0, 2, "fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);"});
    
    if (((integral[0] + integral[1]) == 0) && integral.is_simple())
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
    
    if (((integral[0] + integral[1]) == 0) && (integral.is_simple()))
    {
        lines.push_back({2, 0, 1, "fints[i] += fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;"});
    }
    else
    {
        lines.push_back({2, 0, 2, "const auto ftt = fz_0 * (3.0 - 2.0 * fz_0 * r2ab) * fss;"});
    }
}

void
T2CPrimFuncBodyDriver::_add_nuclear_potential_vars(      VCodeLines&  lines,
                                                   const I2CIntegral& integral) const
{
    lines.push_back({2, 0, 2, "const auto fxi_0 = bra_exp + ket_fe[i];"});
    
    lines.push_back({2, 0, 2, "const auto fe_0 = 1.0 / fxi_0;"});
    
    lines.push_back({2, 0, 2, "const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);"});
    
    if ((integral[0] + integral[1]) == 0)
    {
        lines.push_back({2, 0, 1, "fints[i] += b0_vals[i] * 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);"});
    }
    else
    {
        lines.push_back({2, 0, 2, "const auto fss = 2.0 * charge * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);"});
    }
}

void
T2CPrimFuncBodyDriver::_add_nuclear_potential_geom_vars(      VCodeLines&  lines,
                                                        const I2CIntegral& integral) const
{
    lines.push_back({2, 0, 2, "const auto fxi_0 = bra_exp + ket_fe[i];"});
    
    lines.push_back({2, 0, 2, "const auto fe_0 = 1.0 / fxi_0;"});
    
    lines.push_back({2, 0, 2, "const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);"});
    
    lines.push_back({2, 0, 2, "const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);"});

    lines.push_back({2, 0, 2, "const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);"});

    lines.push_back({2, 0, 2, "const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);"});
    
    lines.push_back({2, 0, 2, "const auto fss = 2.0 * std::sqrt(fxi_0 / fpi) * bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);"});
    
    const auto op_gdrv = integral.integrand().shape().order();

    if ((integral[0] + integral[1]) == 0)
    {
        if (op_gdrv == 1)
        {
            lines.push_back({2, 0, 2, "fints_x[i] += dip_x * b1_vals[i] * 2.0 * fxi_0 * rpc_x * fss;"});
            
            lines.push_back({2, 0, 2, "fints_y[i] += dip_y * b1_vals[i] * 2.0 * fxi_0 * rpc_y * fss;"});
            
            lines.push_back({2, 0, 2, "fints_z[i] += dip_z * b1_vals[i] * 2.0 * fxi_0 * rpc_z * fss;"});
        }
        
        if (op_gdrv == 2)
        {
            lines.push_back({2, 0, 2, "fints_xx[i] += qpol_xx * fss * (4.0 * fxi_0 * fxi_0 * rpc_x * rpc_x * b2_vals[i] - 2.0 * fxi_0 * b1_vals[i]);"});
            
            lines.push_back({2, 0, 2, "fints_xy[i] += qpol_xy * fss * 4.0 * fxi_0 * fxi_0 * rpc_x * rpc_y * b2_vals[i];"});
            
            lines.push_back({2, 0, 2, "fints_xz[i] += qpol_xz * fss * 4.0 * fxi_0 * fxi_0 * rpc_x * rpc_z * b2_vals[i];"});
            
            lines.push_back({2, 0, 2, "fints_yy[i] += qpol_yy * fss * (4.0 * fxi_0 * fxi_0 * rpc_y * rpc_y * b2_vals[i] - 2.0 * fxi_0 * b1_vals[i]);"});
            
            lines.push_back({2, 0, 2, "fints_yz[i] += qpol_yz * fss * 4.0 * fxi_0 * fxi_0 * rpc_y * rpc_z * b2_vals[i];"});
            
            lines.push_back({2, 0, 2, "fints_zz[i] += qpol_zz * fss * (4.0 * fxi_0 * fxi_0 * rpc_z * rpc_z * b2_vals[i] - 2.0 * fxi_0 * b1_vals[i]);"});
        }
    }
    else
    {
        if (op_gdrv == 1)
        {
            lines.push_back({2, 0, 2, "const auto faa_x = dip_x * 2.0 * fxi_0 * rpc_x * fss;"});
            
            lines.push_back({2, 0, 2, "const auto faa_y = dip_y * 2.0 * fxi_0 * rpc_y * fss;"});
            
            lines.push_back({2, 0, 2, "const auto faa_z = dip_z * 2.0 * fxi_0 * rpc_z * fss;"});
        }
        
        if (op_gdrv == 2)
        {
            lines.push_back({2, 0, 2, "const auto faa_xx = qpol_xx * fss * 4.0 * fxi_0 * fxi_0 * rpc_x * rpc_x;"});
            
            lines.push_back({2, 0, 2, "const auto faa_xy = qpol_xy * fss * 4.0 * fxi_0 * fxi_0 * rpc_x * rpc_y;"});
            
            lines.push_back({2, 0, 2, "const auto faa_xz = qpol_xz * fss * 4.0 * fxi_0 * fxi_0 * rpc_x * rpc_z;"});
            
            lines.push_back({2, 0, 2, "const auto faa_yy = qpol_yy * fss * 4.0 * fxi_0 * fxi_0 * rpc_y * rpc_y;"});
            
            lines.push_back({2, 0, 2, "const auto faa_yz = qpol_yz * fss * 4.0 * fxi_0 * fxi_0 * rpc_y * rpc_z;"});
            
            lines.push_back({2, 0, 2, "const auto faa_zz = qpol_zz * fss * 4.0 * fxi_0 * fxi_0 * rpc_z * rpc_z;"});
            
            lines.push_back({2, 0, 2, "const auto faa_x = 2.0 * fxi_0 * rpc_x * fss;"});
            
            lines.push_back({2, 0, 2, "const auto faa_y = 2.0 * fxi_0 * rpc_y * fss;"});
            
            lines.push_back({2, 0, 2, "const auto faa_z = 2.0 * fxi_0 * rpc_z * fss;"});
            
            lines.push_back({2, 0, 2, "const auto faa = -2.0 * fxi_0 * fss;"});
        }
    }
}

void
T2CPrimFuncBodyDriver::_add_multipole_vars(      VCodeLines&  lines,
                                           const I2CIntegral& integral) const
{
    lines.push_back({2, 0, 2, "const auto fxi_0 = bra_exp + ket_fe[i];"});
        
    lines.push_back({2, 0, 2, "const auto fe_0 = 1.0 / fxi_0;"});
        
    lines.push_back({2, 0, 2, "const auto fz_0 = bra_exp * ket_fe[i] * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);"});
        
    lines.push_back({2, 0, 2, "const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);"});

    lines.push_back({2, 0, 2, "const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);"});

    lines.push_back({2, 0, 2, "const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);"});
        
    lines.push_back({2, 0, 2, "const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);"});
        
    const auto op_gdrv = integral.integrand().shape().order();

    if ((integral[0] + integral[1]) == 0)
    {
        if (op_gdrv == 1)
        {
            lines.push_back({2, 0, 2, "fints_x[i] += rpc_x * fss;"});
                
            lines.push_back({2, 0, 2, "fints_y[i] += rpc_y * fss;"});
                
            lines.push_back({2, 0, 2, "fints_z[i] += rpc_z * fss;"});
        }
            
        if (op_gdrv == 2)
        {
            lines.push_back({2, 0, 2, "fints_xx[i] += fss * (rpc_x * rpc_x + 0.5 * fe_0);"});
                
            lines.push_back({2, 0, 2, "fints_xy[i] += fss * rpc_x * rpc_y;"});
                
            lines.push_back({2, 0, 2, "fints_xz[i] += fss * rpc_x * rpc_z;"});
                
            lines.push_back({2, 0, 2, "fints_yy[i] += fss * (rpc_y * rpc_y + 0.5 * fe_0);"});
                
            lines.push_back({2, 0, 2, "fints_yz[i] += fss * rpc_y * rpc_z;"});
                
            lines.push_back({2, 0, 2, "fints_zz[i] += fss * (rpc_z * rpc_z + 0.5 * fe_0);"});
        }
        
        if (op_gdrv == 3)
        {
            lines.push_back({2, 0, 2, "fints_xxx[i] += fss * (rpc_x * rpc_x * rpc_x + 1.5 * fe_0 * rpc_x);"});
                
            lines.push_back({2, 0, 2, "fints_xxy[i] += fss * (rpc_x * rpc_x * rpc_y + 0.5 * fe_0 * rpc_y);"});
                
            lines.push_back({2, 0, 2, "fints_xxz[i] += fss * (rpc_x * rpc_x * rpc_z + 0.5 * fe_0 * rpc_z);"});
                
            lines.push_back({2, 0, 2, "fints_xyy[i] += fss * (rpc_x * rpc_y * rpc_y + 0.5 * fe_0 * rpc_x);"});
                
            lines.push_back({2, 0, 2, "fints_xyz[i] += fss * rpc_x * rpc_y * rpc_z;"});
                
            lines.push_back({2, 0, 2, "fints_xzz[i] += fss * (rpc_x * rpc_z * rpc_z + 0.5 * fe_0 * rpc_x);"});
            
            lines.push_back({2, 0, 2, "fints_yyy[i] += fss * (rpc_y * rpc_y * rpc_y + 1.5 * fe_0 * rpc_y);"});
                
            lines.push_back({2, 0, 2, "fints_yyz[i] += fss * (rpc_y * rpc_y * rpc_z + 0.5 * fe_0 * rpc_z);"});
                
            lines.push_back({2, 0, 2, "fints_yzz[i] += fss * (rpc_y * rpc_z * rpc_z + 0.5 * fe_0 * rpc_y);"});
            
            lines.push_back({2, 0, 2, "fints_zzz[i] += fss * (rpc_z * rpc_z * rpc_z + 1.5 * fe_0 * rpc_z);"});
        }
    }
    else
    {
        if (op_gdrv >= 1)
        {
            lines.push_back({2, 0, 2, "const auto faa_x = rpc_x * fss;"});
                
            lines.push_back({2, 0, 2, "const auto faa_y = rpc_y * fss;"});
                
            lines.push_back({2, 0, 2, "const auto faa_z = rpc_z * fss;"});
        }
            
        if (op_gdrv >= 2)
        {
            lines.push_back({2, 0, 2, "const auto faa_xx = fss * (rpc_x * rpc_x + 0.5 * fe_0);"});
                    
            lines.push_back({2, 0, 2, "const auto faa_xy = fss * rpc_x * rpc_y;"});
                    
            lines.push_back({2, 0, 2, "const auto faa_xz = fss * rpc_x * rpc_z;"});
                    
            lines.push_back({2, 0, 2, "const auto faa_yy = fss * (rpc_y * rpc_y + 0.5 * fe_0);"});
                    
            lines.push_back({2, 0, 2, "const auto faa_yz = fss * rpc_y * rpc_z;"});
                    
            lines.push_back({2, 0, 2, "const auto faa_zz = fss * (rpc_z * rpc_z + 0.5 * fe_0);"});
        }
        
        if (op_gdrv >= 3)
        {
            lines.push_back({2, 0, 2, "const auto faa_xxx = fss * (rpc_x * rpc_x * rpc_x + 1.5 * fe_0 * rpc_x);"});
                
            lines.push_back({2, 0, 2, "const auto faa_xxy = fss * (rpc_x * rpc_x * rpc_y + 0.5 * fe_0 * rpc_y);"});
                
            lines.push_back({2, 0, 2, "const auto faa_xxz = fss * (rpc_x * rpc_x * rpc_z + 0.5 * fe_0 * rpc_z);"});
                
            lines.push_back({2, 0, 2, "const auto faa_xyy = fss * (rpc_x * rpc_y * rpc_y + 0.5 * fe_0 * rpc_x);"});
                
            lines.push_back({2, 0, 2, "const auto faa_xyz = fss * rpc_x * rpc_y * rpc_z;"});
                
            lines.push_back({2, 0, 2, "const auto faa_xzz = fss * (rpc_x * rpc_z * rpc_z + 0.5 * fe_0 * rpc_x);"});
            
            lines.push_back({2, 0, 2, "const auto faa_yyy = fss * (rpc_y * rpc_y * rpc_y + 1.5 * fe_0 * rpc_y);"});
                
            lines.push_back({2, 0, 2, "const auto faa_yyz = fss * (rpc_y * rpc_y * rpc_z + 0.5 * fe_0 * rpc_z);"});
                
            lines.push_back({2, 0, 2, "const auto faa_yzz = fss * (rpc_y * rpc_z * rpc_z + 0.5 * fe_0 * rpc_y);"});
            
            lines.push_back({2, 0, 2, "const auto faa_zzz = fss * (rpc_z * rpc_z * rpc_z + 1.5 * fe_0 * rpc_z);"});
        }
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
    
    _add_prefactors(lines, rgroup, integral);
        
    _add_simd_lines(lines, labels, rgroup);
}

R2Group
T2CPrimFuncBodyDriver::_generate_integral_group(const VT2CIntegrals& components,
                                                const I2CIntegral&   integral) const
{
    R2Group rgroup;
    
    if (!integral.is_simple())
    {
        T2CCenterDriver t2c_geom_drv;
        
        rgroup = t2c_geom_drv.create_recursion(components);
    }
    
    // Overlap inntegrals
    
    if (integral.integrand() == Operator("1"))
    {
        T2COverlapDriver t2c_ovl_drv;
        
        if (integral.is_simple())
        {
            rgroup = t2c_ovl_drv.create_recursion(components);
        }
        else
        {
            t2c_ovl_drv.apply_recursion(rgroup);
        }
    }
    
    // Kinetic energy inntegrals
    
    if (integral.integrand() == Operator("T"))
    {
        T2CKineticEnergyDriver t2c_kin_drv;
        
        if (integral.is_simple())
        {
            rgroup = t2c_kin_drv.create_recursion(components);
        }
        else
        {
            t2c_kin_drv.apply_recursion(rgroup);
        }
    }
    
    // Nuclear potential inntegrals
    
    if (const auto integrand = integral.integrand();
        (integrand.name() == "A") && (integrand.shape() == Tensor(0)))
    {
        T2CNuclearPotentialDriver t2c_npot_drv;
        
        rgroup = t2c_npot_drv.create_recursion(components);
    }
    
    // Nuclear potential geometrical derivative inntegrals
    
    if (const auto integrand = integral.integrand();
        (integrand.name() == "AG") && (integrand.shape() != Tensor(0)))
    {
        T2CNuclearPotentialGeometryDriver t2c_npot_geom_drv;
        
        rgroup = t2c_npot_geom_drv.create_recursion(components);
    }
    
    // Multipole inntegrals
    
    if (const auto integrand = integral.integrand();
        (integrand.name() == "r") && (integrand.shape() != Tensor(0)))
    {
        T2CMultipoleDriver t2c_mpol_drv;
        
        rgroup = t2c_mpol_drv.create_recursion(components);
    }
    
    // ... other integrals
    
    return rgroup;
}

void
T2CPrimFuncBodyDriver::_add_prefactors(      VCodeLines&  lines,
                                       const R2Group&     rgroup,
                                       const I2CIntegral& integral) const
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
    
    if (t2c::find_factor(rgroup, "tke_0"))
    {
        lines.push_back({2, 0, 2, "const auto tke_0 = ket_fe[i];"});
    }
    
    if (t2c::find_factor(rgroup, "tbe_0"))
    {
        lines.push_back({2, 0, 2, "const auto tbe_0 = bra_exp;"});
    }
    
    // special variables, excluded for specific integrals
    
    const auto integrand = integral.integrand();
    
    if ((integrand.name() != "AG") && (integrand.name() != "r"))
    {
        if (t2c::find_factor(rgroup, "rpc_x"))
        {
            lines.push_back({2, 0, 2, "const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);"});
        }
        
        if (t2c::find_factor(rgroup, "rpc_y"))
        {
            lines.push_back({2, 0, 2, "const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);"});
        }

        if (t2c::find_factor(rgroup, "rpc_z"))
        {
            lines.push_back({2, 0, 2, "const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);"});
        }
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
    const auto tlabel = _get_aux_label(integral, rdist.root().integral());
    
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
    
        if (((eterm - sterm) > 1) || (simd_str[0] == '-')) simd_str = "(" + simd_str + ")";

        if (simd_str.empty())
        {
            lines.push_back({2, 0, 2, label + "[i] += " + tlabel + ";"});
        }
        else
        {
            lines.push_back({2, 0, 2, label + "[i] += " + tlabel + " * " + simd_str + ";"});
        }
    }
}

std::string
T2CPrimFuncBodyDriver::_get_aux_label(const T2CIntegral& integral,
                                      const T2CIntegral& base) const
{
    const auto bname = base.integrand().name();
    
    const auto iname = integral.integrand().name();
    
    const auto border = base.integrand().shape().order();
    
    const auto iorder = integral.integrand().shape().order();
    
    if (bname == "1")
    {
        return std::string("fss");
    }
    
    if (bname == "T")
    {
        if (iname == "1")
        {
            return std::string("fss");
        }
        else
        {
            return std::string("ftt");
        }
    }
    
    if (bname == "A")
    {
        return "fss * b" + std::to_string(integral.order()) + "_vals[i]";
    }
    
    if ((bname == "AG") && (border == 1))
    {
        const auto torder = integral.order();
        
        if (iname == "A")
        {
            const auto dlabel = "dip_" + base.integrand().shape().label();
            
            return dlabel + " * fss * b" + std::to_string(torder) + "_vals[i]";
        }
        
        if (iname == "AG")
        {
            const auto flabel = "faa_" + integral.integrand().shape().label();
            
            return flabel + " * b" + std::to_string(torder + 1) + "_vals[i]";
        }
    }
    
    if ((bname == "AG") && (border == 2))
    {
        const auto torder = integral.order();
        
        const auto qlabel = "qpol_" + base.integrand().shape().label();
        
        if (iname == "A")
        {
            return qlabel + " * fss * b" + std::to_string(torder) + "_vals[i]";
        }
        
        if (iname == "AG")
        {
            const auto flabel = "faa_" + integral.integrand().shape().label();
            
            if (iorder == 1)
            {
                return qlabel + " * " + flabel + " * b" + std::to_string(torder + 1) + "_vals[i]";
            }
            
            if (iorder == 2)
            {
                if ((flabel == "faa_xx") || (flabel == "faa_yy") || (flabel == "faa_zz"))
                {
                    return "(" + flabel + " * b" + std::to_string(torder + 2) + "_vals[i] + faa * b" + std::to_string(torder + 1) + "_vals[i])"; 
                }
                else
                {
                    return flabel + " * b" + std::to_string(torder + 2) + "_vals[i]";
                }
            }
        }
    }
    
    if (bname == "r")
    {
        if (iname == "1")
        {
            return std::string("fss");
        }
        else
        {
            return "faa_" + integral.integrand().shape().label();
        }
    }

    return std::string();
}

std::vector<std::string>
T2CPrimFuncBodyDriver::_get_special_vars_str(const I2CIntegral& integral,
                                              const bool        geom_form) const
{
    std::vector<std::string> vstr;
    
    if (t2c::need_boys(integral) || (integral.integrand().name() == "r"))
    {
        if (geom_form)
        {
            vstr.push_back("// set up coordinates for C center");
                
            vstr.push_back("const auto c_rx = point[0];");
                
            vstr.push_back("const auto c_ry = point[1];");
                
            vstr.push_back("const auto c_rz = point[2];");
        }
    }
    
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "AG")
    {
        if (integrand.shape().order() == 1)
        {
            vstr.push_back("// set up dipole components");
                
            vstr.push_back("const auto dip_x = dipole[0];");
                
            vstr.push_back("const auto dip_y = dipole[1];");
                
            vstr.push_back("const auto dip_z = dipole[2];");
        }
        
        if (integrand.shape().order() == 2)
        {
            vstr.push_back("// set up quadrupole components");
                
            vstr.push_back("const auto qpol_xx = quadrupole[0];");
            
            vstr.push_back("const auto qpol_xy = quadrupole[1];");
            
            vstr.push_back("const auto qpol_xz = quadrupole[2];");
            
            vstr.push_back("const auto qpol_yy = quadrupole[3];");
            
            vstr.push_back("const auto qpol_yz = quadrupole[4];");
            
            vstr.push_back("const auto qpol_zz = quadrupole[5];");
        }
    }
    
    return vstr;
}

std::vector<std::string>
T2CPrimFuncBodyDriver::_get_boys_vars_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (t2c::need_boys(integral))
    {
        const auto order = t2c::boys_order(integral);
        
        vstr.push_back("// set up Boys function variables");
                
        vstr.push_back("const CBoysFunc<" + std::to_string(order) + "> bf_table;");
                
        vstr.push_back("alignas(64) TDoubleArray bf_args;");
                
        vstr.push_back("TDoubleArray2D<" + std::to_string(order + 1) + "> bf_values;");
        
        size_t istart = 0;
        
        if (integral.integrand().name() == "AG") istart = 1;
        
        for (size_t i = istart; i < (order + 1); i++)
        {
            vstr.push_back("auto b" + std::to_string(i) + "_vals = bf_values[" + std::to_string(i) + "].data();");
        }
        
        vstr.push_back("auto targs = bf_args.data();"); 
    }
    
    return vstr;
}

void
T2CPrimFuncBodyDriver::_add_boys_compute_lines(      VCodeLines&  lines,
                                               const I2CIntegral& integral) const
{
    // nuclear potential integrals
    
    if (t2c::need_boys(integral))
    {
        lines.push_back({1, 0, 2, "// compute Boys function values"});
        
        lines.push_back({1, 0, 1, "#pragma omp simd aligned(targs, ket_fe, ket_rx, ket_ry, ket_rz : 64)"});
        
        lines.push_back({1, 0, 1, "for (int64_t i = 0; i < ket_dim; i++)"});
        
        lines.push_back({1, 0, 1, "{"});
        
        lines.push_back({2, 0, 2, "const auto fxi_0 = bra_exp + ket_fe[i];"});
        
        lines.push_back({2, 0, 2, "const auto fe_0 = 1.0 / fxi_0;"});
        
        lines.push_back({2, 0, 2, "const auto rpc_x = fe_0 * (bra_exp * bra_rx + ket_fe[i] * ket_rx[i] - fxi_0 * c_rx);"});

        lines.push_back({2, 0, 2, "const auto rpc_y = fe_0 * (bra_exp * bra_ry + ket_fe[i] * ket_ry[i] - fxi_0 * c_ry);"});

        lines.push_back({2, 0, 2, "const auto rpc_z = fe_0 * (bra_exp * bra_rz + ket_fe[i] * ket_rz[i] - fxi_0 * c_rz);"});
        
        lines.push_back({2, 0, 1, "targs[i] = fxi_0 * (rpc_x * rpc_x + rpc_y * rpc_y + rpc_z * rpc_z);"});
        
        lines.push_back({1, 0, 2, "}"});
        
        const auto order = t2c::boys_order(integral);
        
        lines.push_back({1, 0, 2, "bf_table.compute<" + std::to_string(order + 1) + ">(bf_values, bf_args, ket_dim);"});
    }
}
