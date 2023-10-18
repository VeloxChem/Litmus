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
        _add_coords_compute(lines);
        
        for (const auto& label : _get_boys_vars_str(integral))
        {
            lines.push_back({1, 0, 2, label});
        }
        
        _add_boys_compute_lines(lines, integral); 
    }
    
    lines.push_back({1, 0, 2, "// set up pointer to integrals buffer"});
    
    lines.push_back({1, 0, 2, "auto fints = buffer.data();"});
    
    lines.push_back({1, 0, 2, "// compute electron repulsion integrals"});
    
    if (diagonal)
    {
        _add_func_pragma(lines, integral);
        
        _add_loop_start(lines, integral);
        
        _add_simd_code(lines, component, integral);
        
        _add_loop_end(lines);
    }
    else
    {
        _add_split_simd_code(lines, component, integral);
    }
    
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

void
T4CDiagPrimFuncBodyDriver::_add_coords_compute(VCodeLines&  lines) const
{
    lines.push_back({1, 0, 2, "// set up P and Q center coordinates"});
    
    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray coords_p_x;"});
    
    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray coords_p_y;"});
    
    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray coords_p_z;"});
    
    lines.push_back({1, 0, 2, "auto rp_x = coords_p_x.data();"});

    lines.push_back({1, 0, 2, "auto rp_y = coords_p_y.data();"});
    
    lines.push_back({1, 0, 2, "auto rp_z = coords_p_z.data();"});
    
    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray coords_q_x;"});
    
    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray coords_q_y;"});
    
    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray coords_q_z;"});
    
    lines.push_back({1, 0, 2, "auto rq_x = coords_q_x.data();"});

    lines.push_back({1, 0, 2, "auto rq_y = coords_q_y.data();"});
    
    lines.push_back({1, 0, 2, "auto rq_z = coords_q_z.data();"});
    
    lines.push_back({1, 0, 2, "// compute P and Q center coordinates"});
    
    lines.push_back({1, 0, 1, "#pragma omp simd aligned(rp_x, rp_y, rp_z, rq_x, rq_y, rq_z, ra_x, ra_y, ra_z, rb_x, rb_y, rb_z, fexps_a, fexps_b, fexps_c, fexps_d : 64)"});
    
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < ndim; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "const auto fi_ab_0 = 1.0 / (fexps_a[i] + fexps_b[i]);"});
    
    lines.push_back({2, 0, 2, "const auto fi_cd_0 = 1.0 / (fexps_c[i] + fexps_d[i]);"});
    
    lines.push_back({2, 0, 2, "rp_x[i] = fi_ab_0 * (fexps_a[i] * ra_x[i] + fexps_b[i] * rb_x[i]);"});

    lines.push_back({2, 0, 2, "rp_y[i] = fi_ab_0 * (fexps_a[i] * ra_y[i] + fexps_b[i] * rb_y[i]);"});

    lines.push_back({2, 0, 2, "rp_z[i] = fi_ab_0 * (fexps_a[i] * ra_z[i] + fexps_b[i] * rb_z[i]);"});
    
    lines.push_back({2, 0, 2, "rq_x[i] = fi_cd_0 * (fexps_c[i] * ra_x[i] + fexps_d[i] * rb_x[i]);"});
    
    lines.push_back({2, 0, 2, "rq_y[i] = fi_cd_0 * (fexps_c[i] * ra_y[i] + fexps_d[i] * rb_y[i]);"});
    
    lines.push_back({2, 0, 1, "rq_z[i] = fi_cd_0 * (fexps_c[i] * ra_z[i] + fexps_d[i] * rb_z[i]);"});
    
    lines.push_back({1, 0, 2, "}"});
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
    
    vstr.push_back("// set up overlap values");
    
    vstr.push_back("alignas(64) TDoubleArray fovl_ab_cd;");
    
    vstr.push_back("auto fss_abcd = fovl_ab_cd.data();");
    
    return vstr;
}

void
T4CDiagPrimFuncBodyDriver::_add_boys_compute_lines(      VCodeLines&  lines,
                                                   const I4CIntegral& integral) const
{
    lines.push_back({1, 0, 2, "// compute Boys function and overlap values"});
    
    lines.push_back({1, 0, 1, "#pragma omp simd aligned(targs, fss_abcd, ra_x, ra_y, ra_z, rb_x, rb_y, rb_z, rp_x, rp_y, rp_z, rq_x, rq_y, rq_z, fexps_a, fexps_b, fexps_c, fexps_d, bnorms, knorms : 64)"});
    
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < ndim; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "const auto fe_ab_0 = fexps_a[i] + fexps_b[i];"});
    
    lines.push_back({2, 0, 2, "const auto fe_cd_0 = fexps_c[i] + fexps_d[i];"});
    
    lines.push_back({2, 0, 2, "const auto fi_ab_0 = 1.0 / fe_ab_0;"});

    lines.push_back({2, 0, 2, "const auto fi_cd_0 = 1.0 / fe_cd_0;"});

    lines.push_back({2, 0, 2, "const auto rpq_x = rp_x[i] - rq_x[i];"});
    
    lines.push_back({2, 0, 2, "const auto rpq_y = rp_y[i] - rq_y[i];"});
    
    lines.push_back({2, 0, 2, "const auto rpq_z = rp_z[i] - rq_z[i];"});
    
    lines.push_back({2, 0, 2, "targs[i] = fe_ab_0 * fe_cd_0 * (rpq_x * rpq_x + rpq_y * rpq_y + rpq_z * rpq_z) / (fe_ab_0 + fe_cd_0);"});
    
    lines.push_back({2, 0, 2, "const auto ab_x = ra_x[i] - rb_x[i];"});
    
    lines.push_back({2, 0, 2, "const auto ab_y = ra_y[i] - rb_y[i];"});

    lines.push_back({2, 0, 2, "const auto ab_z = ra_z[i] - rb_z[i];"});

    lines.push_back({2, 0, 2, "fss_abcd[i] = 4.0 * bnorms[i] * knorms[i] * std::pow(fi_ab_0 * fi_cd_0 * fpi * fpi, 1.50)"});

    lines.push_back({2, 0, 2, "            * std::exp(-(fexps_a[i] * fexps_b[i] * fi_ab_0 + fexps_c[i] * fexps_d[i] * fi_cd_0) * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z))"});

    lines.push_back({2, 0, 1, "            * std::sqrt(invfpi * fe_ab_0 * fe_cd_0 / (fe_ab_0 + fe_cd_0));"});

    lines.push_back({1, 0, 2, "}"});
    
    const auto order = t4c::boys_order(integral);
    
    lines.push_back({1, 0, 2, "bf_table.compute<" + std::to_string(order + 1) + ">(bf_values, bf_args, ndim);"});
}

void
T4CDiagPrimFuncBodyDriver::_add_func_pragma(      VCodeLines&  lines,
                                            const I4CIntegral& integral) const
{
    std::vector<std::string> labels({"fints", "ra_x", "ra_y", "ra_z", "rb_x", "rb_y", "rb_z",
                                     "fexps_a", "fexps_b", "bnorms"});
    
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
                                           const I4CIntegral& integral) const
{
    
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < ndim; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "const auto ab_x = ra_x[i] - rb_x[i];"});
    
    lines.push_back({2, 0, 2, "const auto ab_y = ra_y[i] - rb_y[i];"});
    
    lines.push_back({2, 0, 2, "const auto ab_z = ra_z[i] - rb_z[i];"});
    
    lines.push_back({2, 0, 2, "const auto fe_ab_0 = fexps_a[i] + fexps_b[i];"});

    lines.push_back({2, 0, 2, "const auto fi_ab_0 = 1.0 / fe_ab_0;"});

    lines.push_back({2, 0, 2, "const auto fz_ab_0 = fexps_a[i] * fexps_b[i] * fi_ab_0;"});

    lines.push_back({2, 0, 2, "const auto fss_ab = bnorms[i] * std::pow(fi_ab_0 * fpi, 1.50) * std::exp(-fz_ab_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z));"});
}


void
T4CDiagPrimFuncBodyDriver::_add_simd_code(      VCodeLines&  lines,
                                          const T4CIntegral& component,
                                          const I4CIntegral& integral) const
{
    const auto rdist = (_generate_integral_group(component, integral, true))[0];
  
    _add_prefactors(lines, rdist);
        
    size_t index = 0;
        
    for (const auto& tint : rdist.unique_integrals())
    {
        const auto tdist = rdist.split(tint);
            
        _add_simd_lines_block(lines, tint, tdist, index, true);
            
        index++;
    }
}

void
T4CDiagPrimFuncBodyDriver::_add_split_simd_code(      VCodeLines&  lines,
                                                const T4CIntegral& component,
                                                const I4CIntegral& integral) const
{
    const auto rdist = (_generate_integral_group(component, integral, false))[0];
        
    for (const auto& tint : rdist.unique_integrals())
    {
        const auto tdist = rdist.split(tint);
            
        _add_split_simd_block(lines, tint, tdist);
    }
}

void
T4CDiagPrimFuncBodyDriver::_add_loop_end(VCodeLines& lines) const
{
    lines.push_back({2, 0, 1, "fints[i] += 2.0 * fss_ab * fss_ab * std::sqrt(0.5 * fe_ab_0 * invfpi) * fact;"});
    
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
            simd_str += t4c::get_factor_label(rdist[j], integral, j == sterm, diagonal);
        }
        
        if (simd_str.size() > 3)
        {
            if (simd_str[1] == '+') simd_str.erase(0, 3);
        }
        
        std::string var_str = ((index == 0) && (i == 0)) ? "auto fact = " : "fact += ";

        if (simd_str.empty())
        {
            lines.push_back({2, 0, 2, var_str + "1.0;"});
        }
        else
        {
            lines.push_back({2, 0, 2, var_str  + simd_str + ";"});
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

void
T4CDiagPrimFuncBodyDriver::_add_split_simd_block(      VCodeLines&  lines,
                                                 const T4CIntegral& integral,
                                                 const R4CDist&     rdist) const
{
    const auto order = integral.order();
    
    lines.push_back({1, 0, 2, "// add Boys order " + std::to_string(order) + " contributions"});
    
    _add_split_pragma(lines, integral, rdist);
    
    _add_split_loop_start(lines, integral, rdist);
    
    _add_simd_lines_block(lines, integral, rdist, 0, false);
    
    _add_split_loop_end(lines, integral);
}

void
T4CDiagPrimFuncBodyDriver::_add_split_pragma(      VCodeLines&  lines,
                                             const T4CIntegral& integral,
                                             const R4CDist&     rdist) const
{
    //const auto order = integral.order();
    
    std::string vars_str;
    
    if (t4c::find_factor(rdist, "fti_ab_0") ||
        t4c::find_factor(rdist, "fti_cd_0") ||
        t4c::find_factor(rdist, "fi_ab_0")  ||
        t4c::find_factor(rdist, "rwp_x")    ||
        t4c::find_factor(rdist, "rwp_y")    ||
        t4c::find_factor(rdist, "rwp_z")    ||
        t4c::find_factor(rdist, "rwq_x")    ||
        t4c::find_factor(rdist, "rwq_y")    ||
        t4c::find_factor(rdist, "rwq_z"))
    {
        vars_str += " fexps_a, fexps_b,";
    }
    
    if (t4c::find_factor(rdist, "fti_ab_0") ||
        t4c::find_factor(rdist, "fti_cd_0") ||
        t4c::find_factor(rdist, "fi_cd_0")  ||
        t4c::find_factor(rdist, "rwp_x")    ||
        t4c::find_factor(rdist, "rwp_y")    ||
        t4c::find_factor(rdist, "rwp_z")    ||
        t4c::find_factor(rdist, "rwq_x")    ||
        t4c::find_factor(rdist, "rwq_y")    ||
        t4c::find_factor(rdist, "rwq_z"))
    {
        vars_str += " fexps_c, fexps_d,";
    }
    
    if (t4c::find_factor(rdist, "rwp_x") || t4c::find_factor(rdist, "rwq_x"))
    {
        vars_str += " rp_x, rq_x,";
    }
    
    if (t4c::find_factor(rdist, "rwp_y") || t4c::find_factor(rdist, "rwq_y"))
    {
        vars_str += " rp_y, rq_y,";
    }
    
    if (t4c::find_factor(rdist, "rwp_z") || t4c::find_factor(rdist, "rwq_z"))
    {
        vars_str += " rp_z, rq_z,";
    }
    
    if (t4c::find_factor(rdist, "rpa_x") || t4c::find_factor(rdist, "rqc_x"))
    {
        vars_str += " ra_x,";
    }
    
    if (t4c::find_factor(rdist, "rpa_y") || t4c::find_factor(rdist, "rqc_y"))
    {
        vars_str += " ra_y,";
    }
    
    if (t4c::find_factor(rdist, "rpa_z") || t4c::find_factor(rdist, "rqc_z"))
    {
        vars_str += " ra_z,";
    }
    
    if (t4c::find_factor(rdist, "rpb_x") || t4c::find_factor(rdist, "rqd_x"))
    {
        vars_str += " rb_x,";
    }
    
    if (t4c::find_factor(rdist, "rpb_y") || t4c::find_factor(rdist, "rqd_y"))
    {
        vars_str += " rb_y,";
    }
    
    if (t4c::find_factor(rdist, "rpb_z") || t4c::find_factor(rdist, "rqd_z"))
    {
        vars_str += " rb_z,";
    }
    
    lines.push_back({1, 0, 1, "#pragma omp simd aligned(fints, " + vars_str + " fss_abcd : 64)"});
}

void
T4CDiagPrimFuncBodyDriver::_add_split_loop_start(      VCodeLines&  lines,
                                                 const T4CIntegral& integral,
                                                 const R4CDist&     rdist) const
{
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < ndim; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    if (t4c::find_factor(rdist, "fi_ab_0")  ||
        t4c::find_factor(rdist, "fi_abcd_0") ||
        t4c::find_factor(rdist, "fti_ab_0")  ||
        t4c::find_factor(rdist, "fti_cd_0")  ||
        t4c::find_factor(rdist, "rwp_x")     ||
        t4c::find_factor(rdist, "rwp_y")     ||
        t4c::find_factor(rdist, "rwp_z")     ||
        t4c::find_factor(rdist, "rwq_x")     ||
        t4c::find_factor(rdist, "rwq_y")     ||
        t4c::find_factor(rdist, "rwq_z"))
    {
        lines.push_back({2, 0, 2, "const auto fe_ab_0 = fexps_a[i] + fexps_b[i];"});
    }
    
    if (t4c::find_factor(rdist, "fi_cd_0")  ||
        t4c::find_factor(rdist, "fi_abcd_0") ||
        t4c::find_factor(rdist, "fti_ab_0")  ||
        t4c::find_factor(rdist, "fti_cd_0")  ||
        t4c::find_factor(rdist, "rwp_x")     ||
        t4c::find_factor(rdist, "rwp_y")     ||
        t4c::find_factor(rdist, "rwp_z")     ||
        t4c::find_factor(rdist, "rwq_x")     ||
        t4c::find_factor(rdist, "rwq_y")     ||
        t4c::find_factor(rdist, "rwq_z"))
    {
        lines.push_back({2, 0, 2, "const auto fe_cd_0 = fexps_c[i] + fexps_d[i];"});
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
    
    if (t4c::find_factor(rdist, "fi_ab_0") ||  t4c::find_factor(rdist, "fti_ab_0"))
    {
        lines.push_back({2, 0, 2, "const auto fi_ab_0 = 1.0 / fe_ab_0;"});
    }
    
    if (t4c::find_factor(rdist, "fi_cd_0") || t4c::find_factor(rdist, "fti_cd_0"))
    {
        lines.push_back({2, 0, 2, "const auto fi_cd_0 = 1.0 / fe_cd_0;"});
    }
    
    if (t4c::find_factor(rdist, "fti_ab_0"))
    {
        lines.push_back({2, 0, 2, "const auto fti_ab_0 = fe_cd_0 * fi_ab_0 * fi_abcd_0;"});
    }
    
    if (t4c::find_factor(rdist, "fti_cd_0"))
    {
        lines.push_back({2, 0, 2, "const auto fti_cd_0 = fe_ab_0 * fi_cd_0 * fi_abcd_0;"});
    }
    
    if (t4c::find_factor(rdist, "rwp_x"))
    {
        lines.push_back({2, 0, 2, "const auto rwp_x = fe_cd_0 * fi_abcd_0 * (rq_x[i] - rp_x[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rwp_y"))
    {
        lines.push_back({2, 0, 2, "const auto rwp_y = fe_cd_0 * fi_abcd_0 * (rq_y[i] - rp_y[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rwp_z"))
    {
        lines.push_back({2, 0, 2, "const auto rwp_z = fe_cd_0 * fi_abcd_0 * (rq_z[i] - rp_z[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rwq_x"))
    {
        lines.push_back({2, 0, 2, "const auto rwq_x = fe_ab_0 * fi_abcd_0 * (rp_x[i] - rq_x[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rwq_y"))
    {
        lines.push_back({2, 0, 2, "const auto rwq_y = fe_ab_0 * fi_abcd_0 * (rp_y[i] - rq_y[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rwq_z"))
    {
        lines.push_back({2, 0, 2, "const auto rwq_z = fe_ab_0 * fi_abcd_0 * (rp_z[i] - rq_z[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rpa_x"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_x = rp_x[i] - ra_x[i];"});
    }
    
    if (t4c::find_factor(rdist, "rpa_y"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_y = rp_y[i] - ra_y[i];"});
    }
    
    if (t4c::find_factor(rdist, "rpa_z"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_z = rp_z[i] - ra_z[i];"});
    }
    
    if (t4c::find_factor(rdist, "rpb_x"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_x = rp_x[i] - rb_x[i];"});
    }
    
    if (t4c::find_factor(rdist, "rpb_y"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_y = rp_y[i] - rb_y[i];"});
    }
    
    if (t4c::find_factor(rdist, "rpb_z"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_z = rp_z[i] - rb_z[i];"});
    }
    
    if (t4c::find_factor(rdist, "rqc_x"))
    {
        lines.push_back({2, 0, 2, "const auto rqc_x = rq_x[i] - ra_x[i];"});
    }
    
    if (t4c::find_factor(rdist, "rqc_y"))
    {
        lines.push_back({2, 0, 2, "const auto rqc_y = rq_y[i] - ra_y[i];"});
    }
    
    if (t4c::find_factor(rdist, "rqc_z"))
    {
        lines.push_back({2, 0, 2, "const auto rqc_z = rq_z[i] - ra_z[i];"});
    }
    
    if (t4c::find_factor(rdist, "rqd_x"))
    {
        lines.push_back({2, 0, 2, "const auto rqd_x = rq_x[i] - rb_x[i];"});
    }
    
    if (t4c::find_factor(rdist, "rqd_y"))
    {
        lines.push_back({2, 0, 2, "const auto rqd_y = rq_y[i] - rb_y[i];"});
    }
    
    if (t4c::find_factor(rdist, "rqd_z"))
    {
        lines.push_back({2, 0, 2, "const auto rqd_z = rq_z[i] - rb_z[i];"});
    }
}

void
T4CDiagPrimFuncBodyDriver::_add_split_loop_end(      VCodeLines&  lines,
                                               const T4CIntegral& integral) const
{
    const auto order = integral.order();
    
    lines.push_back({2, 0, 1, "fints[i] += fss_abcd[i] * fact * b" + std::to_string(order) + "_vals[i];"});
    
    lines.push_back({1, 0, 2, "}"});
}
