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
    
    vstr.push_back("const auto bnorms = bra_norms.data();");
    
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
    
    lines.push_back({1, 0, 2, "bf_table.compute<" + std::to_string(order + 1) + ">(bf_values, bf_args, ndim);"});
}

void
T4CDiagPrimFuncBodyDriver::_add_func_pragma(      VCodeLines&  lines,
                                            const I4CIntegral& integral,
                                            const bool         diagonal) const
{
    std::vector<std::string> labels({"fints", "ra_x", "ra_y", "ra_z", "rb_x", "rb_y", "rb_z",
                                     "fexps_a", "fexps_b", "bnorms"});
    
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
        lines.push_back({2, 0, 2, "const auto fe_ab_0 = fexps_a[i] + fexps_b[i];"});

        lines.push_back({2, 0, 2, "const auto fi_ab_0 = 1.0 / fe_ab_0;"});

        lines.push_back({2, 0, 2, "const auto fz_ab_0 = fexps_a[i] * fexps_b[i] * fi_ab_0;"});

        lines.push_back({2, 0, 2, "const auto fss_ab = bnorms[i] * std::pow(fi_ab_0 * fpi, 1.50) * std::exp(-fz_ab_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z));"});
    }
    else
    {
       
        lines.push_back({2, 0, 2, "const auto fe_ab_0 = fexps_a[i] + fexps_b[i];"});
        
        lines.push_back({2, 0, 2, "const auto fe_cd_0 = fexps_c[i] + fexps_d[i];"});

        lines.push_back({2, 0, 2, "const auto fi_ab_0 = 1.0 / fe_ab_0;"});
        
        lines.push_back({2, 0, 2, "const auto fi_cd_0 = 1.0 / fe_cd_0;"});

        lines.push_back({2, 0, 2, "const auto fz_ab_0 = fexps_a[i] * fexps_b[i] * fi_ab_0;"});
        
        lines.push_back({2, 0, 2, "const auto fz_cd_0 = fexps_c[i] * fexps_d[i] * fi_cd_0;"});
        
        lines.push_back({2, 0, 2, "const auto fss_abcd = bnorms[i] * knorms[i] * std::pow(fi_ab_0 * fi_cd_0 * fpi * fpi, 1.50)"});
        
        lines.push_back({2, 0, 2, "                    * std::exp(-(fz_ab_0 + fz_cd_0) * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z));"});
    }
}


void
T4CDiagPrimFuncBodyDriver::_add_simd_code(      VCodeLines&  lines,
                                          const T4CIntegral& component,
                                          const I4CIntegral& integral,
                                          const bool         diagonal) const
{
    const auto rdist = (_generate_integral_group(component, integral, diagonal))[0];
    
    _add_prefactors(lines, rdist);
    
    size_t index = 0;
    
    for (const auto& tint : rdist.unique_integrals())
    {
        const auto tdist = rdist.split(tint);
        
        _add_simd_lines_block(lines, tint, tdist, index, diagonal);
        
        index++;
    }
}

void
T4CDiagPrimFuncBodyDriver::_add_loop_end(      VCodeLines& lines,
                                         const bool        diagonal) const
{
    if (diagonal)
    {
        lines.push_back({2, 0, 1, "fints[i] += 2.0 * fss_ab * fss_ab * std::sqrt(0.5 * fe_ab_0 * invfpi) * fact;"});
    }
    else
    {
        lines.push_back({2, 0, 1, "fints[i] += 4.0 * fss_abcd * std::sqrt(invfpi * fe_ab_0 * fe_cd_0 / (fe_ab_0 + fe_cd_0)) * fact;"});
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

void
T4CDiagPrimFuncBodyDriver::_add_simd_lines_block(      VCodeLines&  lines,
                                                 const T4CIntegral& integral,
                                                 const R4CDist&     rdist,
                                                 const size_t       index,
                                                 const bool         diagonal) const
{
    const auto tlabel = _get_aux_label(integral, rdist.root().integral(), diagonal);
    
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
            simd_str += t4c::get_factor_label(rdist[j], j == sterm);
        }
    
        if (((eterm - sterm) > 1) || (simd_str[0] == '-')) simd_str = "(" + simd_str + ")";
        
        std::string var_str = ((index == 0) && (i == 0)) ? "auto fact = " : "fact += ";

        if (simd_str.empty())
        {
            if (tlabel.empty())
            {
                lines.push_back({2, 0, 2, var_str + "1.0;"});
            }
            else
            {
                lines.push_back({2, 0, 2, var_str + tlabel + ";"});
            }
        }
        else
        {
            if (tlabel.empty())
            {
                lines.push_back({2, 0, 2, var_str  + simd_str + ";"});
            }
            else
            {
                lines.push_back({2, 0, 2, var_str + tlabel + " * " + simd_str + ";"});
            }
        }
    }
}

std::string
T4CDiagPrimFuncBodyDriver::_get_aux_label(const T4CIntegral& integral,
                                          const T4CIntegral& base,
                                          const bool         diagonal) const
{
    const auto order = integral.order();
    
    if (diagonal)
    {
        if (order > 0)
        {
            return "(1.0 / " + std::to_string(2 * order +1) + ".0)";
        }
        
        return std::string();
    }
    else
    {
        return  "b" + std::to_string(order) + "_vals[i]";
    }
}

void
T4CDiagPrimFuncBodyDriver::_add_prefactors(      VCodeLines&  lines,
                                           const R4CDist&     rdist) const
{
    if (t4c::find_factor(rdist, "rpa_x"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_x = -fexps_b[i] * ab_x * fi_ab_0;"});
    }
    
    if (t4c::find_factor(rdist, "rpa_y"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_y = -fexps_b[i] * ab_y * fi_ab_0;"});
    }
    
    if (t4c::find_factor(rdist, "rpa_z"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_z = -fexps_b[i] * ab_z * fi_ab_0;"});
    }
    
    if (t4c::find_factor(rdist, "rpb_x"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_x = fexps_a[i] * ab_x * fi_ab_0;"});
    }
    
    if (t4c::find_factor(rdist, "rpb_y"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_y = fexps_a[i] * ab_y * fi_ab_0;"});
    }
    
    if (t4c::find_factor(rdist, "rpb_z"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_z = fexps_a[i] * ab_z * fi_ab_0;"});
    }
    
    if (t4c::find_factor(rdist, "rqc_x"))
    {
        lines.push_back({2, 0, 2, "const auto rqc_x = -fexps_d[i] * ab_x * fi_cd_0;"});
    }
    
    if (t4c::find_factor(rdist, "rqc_y"))
    {
        lines.push_back({2, 0, 2, "const auto rqc_y = -fexps_d[i] * ab_y * fi_cd_0;"});
    }
    
    if (t4c::find_factor(rdist, "rqc_z"))
    {
        lines.push_back({2, 0, 2, "const auto rqc_z = -fexps_d[i] * ab_z * fi_cd_0;"});
    }
    
    if (t4c::find_factor(rdist, "rqd_x"))
    {
        lines.push_back({2, 0, 2, "const auto rqd_x = fexps_c[i] * ab_x * fi_cd_0;"});
    }
    
    if (t4c::find_factor(rdist, "rqd_y"))
    {
        lines.push_back({2, 0, 2, "const auto rqd_y = fexps_c[i] * ab_y * fi_cd_0;"});
    }
    
    if (t4c::find_factor(rdist, "rqd_z"))
    {
        lines.push_back({2, 0, 2, "const auto rqd_z = fexps_c[i] * ab_z * fi_cd_0;"});
    }
    
    if (t4c::find_factor(rdist, "fi_abcd_0") ||
        t4c::find_factor(rdist, "fti_ab_0")  ||
        t4c::find_factor(rdist, "fti_cd_0")  ||
        t4c::find_factor(rdist, "rwp_x")     ||
        t4c::find_factor(rdist, "rwp_y")     ||
        t4c::find_factor(rdist, "rwp_z")     ||
        t4c::find_factor(rdist, "rwq_x")     ||
        t4c::find_factor(rdist, "rwq_y")     ||
        t4c::find_factor(rdist, "rwq_z"))
    {
        lines.push_back({2, 0, 2, "const auto fi_abcd_0 = 1.0 / (fe_ab_0 + fe_cd_0);"});
    }
    
    if (t4c::find_factor(rdist, "fti_ab_0"))
    {
        lines.push_back({2, 0, 2, "const auto fti_ab_0 = fe_cd_0 * fi_ab_0 * fi_abcd_0;"});
    }
    
    if (t4c::find_factor(rdist, "fti_cd_0"))
    {
        lines.push_back({2, 0, 2, "const auto fti_cd_0 = fe_ab_0 * fi_cd_0 * fi_abcd_0;"});
    }
    
    if (t4c::find_factor(rdist, "rwp_x")     ||
        t4c::find_factor(rdist, "rwp_y")     ||
        t4c::find_factor(rdist, "rwp_z")     ||
        t4c::find_factor(rdist, "rwq_x")     ||
        t4c::find_factor(rdist, "rwq_y")     ||
        t4c::find_factor(rdist, "rwq_z"))
    {
        lines.push_back({2, 0, 2, "const auto fm_ac_0 = fi_ab_0 * fexps_a[i] - fi_cd_0 * fexps_c[i];"});
        
        lines.push_back({2, 0, 2, "const auto fm_bd_0 = fi_ab_0 * fexps_b[i] - fi_cd_0 * fexps_d[i];"});
    }
    
    if (t4c::find_factor(rdist, "rwp_x"))
    {
        lines.push_back({2, 0, 2, "const auto rwp_x = -fe_cd_0 * fi_abcd_0 * (fm_ac_0 * ra_x[i] + fm_bd_0 * rb_x[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rwp_y"))
    {
        lines.push_back({2, 0, 2, "const auto rwp_y = -fe_cd_0 * fi_abcd_0 * (fm_ac_0 * ra_y[i] + fm_bd_0 * rb_y[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rwp_z"))
    {
        lines.push_back({2, 0, 2, "const auto rwp_z = -fe_cd_0 * fi_abcd_0 * (fm_ac_0 * ra_z[i] + fm_bd_0 * rb_z[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rwq_x"))
    {
        lines.push_back({2, 0, 2, "const auto rwq_x = fe_ab_0 * fi_abcd_0 * (fm_ac_0 * ra_x[i] + fm_bd_0 * rb_x[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rwq_y"))
    {
        lines.push_back({2, 0, 2, "const auto rwq_y = fe_ab_0 * fi_abcd_0 * (fm_ac_0 * ra_y[i] + fm_bd_0 * rb_y[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rwq_z"))
    {
        lines.push_back({2, 0, 2, "const auto rwq_z = fe_ab_0 * fi_abcd_0 * (fm_ac_0 * ra_z[i] + fm_bd_0 * rb_z[i]);"});
    }
}
