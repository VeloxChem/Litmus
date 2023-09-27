#include "t4c_diag_prim_body.hpp"

#include "t4c_diag_eri_driver.hpp"
#include "t4c_full_eri_driver.hpp"
#include "t4c_utils.hpp"

void
T4CDiagPrimFuncBodyDriver::write_prim_func_body(      std::ofstream& fstream,
                                                const T4CIntegral&   component,
                                                const I4CIntegral&   integral,
                                                const bool           diagonal) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_common_data_str(diagonal))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    if (!diagonal)
    {
        for (const auto& label : _get_boys_vars_str(integral))
        {
            lines.push_back({1, 0, 2, label});
        }
        
        _add_boys_compute_lines(lines, integral); 
    }
    
    lines.push_back({1, 0, 2, "// set up pointer to integrals buffer"});
    
    lines.push_back({1, 0, 2, "auto fints = buffer.data();"});
    
    _add_func_pragma(lines, integral, diagonal);
    
    _add_loop_start(lines, integral, diagonal);
    
    _add_simd_code(lines, component, integral, diagonal);
    
    _add_loop_end(lines, diagonal);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CDiagPrimFuncBodyDriver::_get_common_data_str(const bool diagonal) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// set up math constants");
        
    vstr.push_back("const auto fpi = mathconst::getPiValue();");
    
    vstr.push_back("const auto invfpi = 1.0 / mathconst::getPiValue();");
    
    vstr.push_back("// set up coordinates for bra center A");

    vstr.push_back("const auto ra_x = coords_a_x.data();");

    vstr.push_back("const auto ra_y = coords_a_y.data();");

    vstr.push_back("const auto ra_z = coords_a_z.data();");
    
    vstr.push_back("// set up coordinates for bra center B");

    vstr.push_back("const auto rb_x = coords_b_x.data();");

    vstr.push_back("const auto rb_y = coords_b_y.data();");

    vstr.push_back("const auto rb_z = coords_b_z.data();");
    
    vstr.push_back("// set up bra side data");
    
    vstr.push_back("const auto fexps_a = bra_exps_a.data();");
    
    vstr.push_back("const auto fexps_b = bra_exps_b.data();");
    
    vstr.push_back("const auto fnorms = bra_norms.data();");
    
    if (!diagonal)
    {
        vstr.push_back("// set up ket side data");
        
        vstr.push_back("const auto fexps_c = ket_exps_c.data();");
        
        vstr.push_back("const auto fexps_d = ket_exps_d.data();");
        
        vstr.push_back("const auto knorms = ket_norms.data();");
    }
    
    return vstr;
}

std::vector<std::string>
T4CDiagPrimFuncBodyDriver::_get_boys_vars_str(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    const auto order = t4c::boys_order(integral);
    
    vstr.push_back("// set up Boys function variables");
            
    vstr.push_back("const CBoysFunc<" + std::to_string(order) + "> bf_table;");
            
    vstr.push_back("alignas(64) TDoubleArray bf_args;");
            
    vstr.push_back("TDoubleArray2D<" + std::to_string(order + 1) + "> bf_values;");
        
    for (size_t i = 0; i < (order + 1); i++)
    {
        vstr.push_back("auto b" + std::to_string(i) + "_vals = bf_values[" + std::to_string(i) + "].data();");
    }
    
    vstr.push_back("auto targs = bf_args.data();");
    
    return vstr;
}

void
T4CDiagPrimFuncBodyDriver::_add_boys_compute_lines(      VCodeLines&  lines,
                                                   const I4CIntegral& integral) const
{
    lines.push_back({1, 0, 2, "// compute Boys function values"});
    
    lines.push_back({1, 0, 1, "#pragma omp simd aligned(targs, ra_x, ra_y, ra_z, rb_x, rb_y, rb_z, fexps_a, fexps_b, fexps_c, fexps_d : 64)"});
    
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < ndim; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "const auto fe_ab = fexps_a[i] + fexps_b[i];"});
    
    lines.push_back({2, 0, 2, "const auto fe_cd = fexps_c[i] + fexps_d[i];"});
    
    lines.push_back({2, 0, 2, "const auto fi_ab = 1.0 / fe_ab;"});

    lines.push_back({2, 0, 2, "const auto fi_cd = 1.0 / fe_cd;"});

    lines.push_back({2, 0, 2, "const auto fm_ac = fi_ab * fexps_a[i] - fi_cd * fexps_c[i];"});
    
    lines.push_back({2, 0, 2, "const auto fm_bd = fi_ab * fexps_b[i] - fi_cd * fexps_d[i];"});
    
    lines.push_back({2, 0, 2, "const auto rpq_x = ra_x[i] * fm_ac + rb_x[i] * fm_bd;"});
    
    lines.push_back({2, 0, 2, "const auto rpq_y = ra_y[i] * fm_ac + rb_y[i] * fm_bd;"});
    
    lines.push_back({2, 0, 2, "const auto rpq_z = ra_z[i] * fm_ac + rb_z[i] * fm_bd;"});
    
    lines.push_back({2, 0, 2, "targs[i] = fe_ab * fe_cd * (rpq_x * rpq_x + rpq_y * rpq_y + rpq_z * rpq_z) / (fe_ab + fe_cd);"});
    
    lines.push_back({1, 0, 2, "}"});
    
    const auto order = t4c::boys_order(integral);
    
    lines.push_back({1, 0, 2, "bf_table.compute<" + std::to_string(order + 1) + ">(bf_values, bf_args, ket_dim);"});
}

void
T4CDiagPrimFuncBodyDriver::_add_func_pragma(      VCodeLines&  lines,
                                            const I4CIntegral& integral,
                                            const bool         diagonal) const
{
    std::vector<std::string> labels({"fints", "ra_x", "ra_y", "ra_z", "rb_x", "rb_y", "rb_z",
                                     "fexps_a", "fexps_b", "fnorms"});
    
    if (!diagonal)
    {
        labels.push_back("fexps_c");
        
        labels.push_back("fexps_d");
        
        labels.push_back("knorms");
        
        const auto order = t4c::boys_order(integral);
        
        for (size_t i = 0; i < (order + 1); i++)
        {
            labels.push_back("b" + std::to_string(i) + "_vals");
        }
        
    }
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        if (i == 0)
        {
            lines.push_back({1, 0, 1, "#pragma omp simd aligned(" + labels[i] + ",\\"});
        }
        else
        {
            if ((i + 1) == labels.size())
            {
                lines.push_back({1, 25, 1, labels[i] + " : 64)"});
            }
            else
            {
                lines.push_back({1, 25, 1, labels[i] + ",\\"});
            }
        }
    }
}

void
T4CDiagPrimFuncBodyDriver::_add_loop_start(      VCodeLines&  lines,
                                           const I4CIntegral& integral,
                                           const bool         diagonal) const
{
    
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < ndim; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "const auto ab_x = ra_x[i] - rb_x[i];"});
    
    lines.push_back({2, 0, 2, "const auto ab_y = ra_y[i] - rb_y[i];"});
    
    lines.push_back({2, 0, 2, "const auto ab_z = ra_z[i] - rb_z[i];"});
    
    if (diagonal)
    {
        lines.push_back({2, 0, 2, "const auto fe_ab = fexps_a[i] + fexps_b[i];"});

        lines.push_back({2, 0, 2, "const auto fi_ab = 1.0 / fe_ab;"});

        lines.push_back({2, 0, 2, "const auto fz_ab = fexps_a[i] * fexps_b[i] * fi_ab;"});

        lines.push_back({2, 0, 2, "const auto fss_ab = fnorms[i] * std::pow(fi_ab * fpi, 1.50) * std::exp(-fz_ab * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z));"});
    }
    else
    {
       
        lines.push_back({2, 0, 2, "const auto fe_ab = fexps_a[i] + fexps_b[i];"});
        
        lines.push_back({2, 0, 2, "const auto fe_cd = fexps_c[i] + fexps_d[i];"});

        lines.push_back({2, 0, 2, "const auto fi_ab = 1.0 / fe_ab;"});
        
        lines.push_back({2, 0, 2, "const auto fi_cd = 1.0 / fe_cd;"});

        lines.push_back({2, 0, 2, "const auto fz_ab = fexps_a[i] * fexps_b[i] * fi_ab;"});
        
        lines.push_back({2, 0, 2, "const auto fz_cd = fexps_c[i] * fexps_d[i] * fi_cd;"});
        
        lines.push_back({2, 0, 2, "const auto fss_abcd = bnorms[i] * knorms[i] * std::pow(fi_ab * fi_cd * fpi * fpi, 1.50)"});
        
        lines.push_back({2, 0, 2, "                    * std::exp(-(fz_ab + fz_cd) * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z));"});
    }
}


void
T4CDiagPrimFuncBodyDriver::_add_simd_code(      VCodeLines&  lines,
                                          const T4CIntegral& component,
                                          const I4CIntegral& integral,
                                          const bool         diagonal) const
{
    auto rgroup = _generate_integral_group(component, integral, diagonal);
        
    //t4c::debug_info(rgroup[0]);
    
}

void
T4CDiagPrimFuncBodyDriver::_add_loop_end(      VCodeLines& lines,
                                         const bool        diagonal) const
{
    if (diagonal)
    {
        lines.push_back({2, 0, 1, "fints[i] += 2.0 * fss_ab * fss_ab * std::sqrt(0.5 * fe_ab * invfpi) * fact;"});
    }
    else
    {
        lines.push_back({2, 0, 1, "fints[i] += 4.0 * fss_abcd * std::sqrt(invfpi * fe_ab * fe_cd / (fe_ab + fe_cd)) * fact;"});
    }
    
    lines.push_back({1, 0, 1, "}"});
}



R4Group
T4CDiagPrimFuncBodyDriver::_generate_integral_group(const T4CIntegral& component,
                                                    const I4CIntegral& integral,
                                                    const bool         diagonal) const
{
    R4Group rgroup;
    
    if (diagonal)
    {
        T4CDiagElectronRepulsionDriver t4c_eri_drv;
        
        rgroup = t4c_eri_drv.create_recursion({component, });
    }
    else
    {
        T4CFullElectronRepulsionDriver t4c_eri_drv;
        
        rgroup = t4c_eri_drv.create_recursion({component, });
    }
    
    return rgroup;
}
