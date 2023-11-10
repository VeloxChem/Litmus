#include "t4c_full_prim_body.hpp"

#include "t4c_utils.hpp"
#include "t4c_full_eri_driver.hpp"

void
T4CFullPrimFuncBodyDriver::write_prim_func_body(      std::ofstream& fstream,
                                                const T4CIntegral&   component,
                                                const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_common_data_str())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_coords_compute(lines);
    
    for (const auto& label : _get_boys_vars_str(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_boys_compute_lines(lines, integral);
    
    lines.push_back({1, 0, 2, "// set up pointer to integrals buffer"});
    
    lines.push_back({1, 0, 2, "auto fints = buffer.data();"});
    
    lines.push_back({1, 0, 2, "// compute electron repulsion integrals"});
    
    _add_split_simd_code(lines, component, integral);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

void
T4CFullPrimFuncBodyDriver::write_vrr_func_body(      std::ofstream& fstream,
                                               const T4CIntegral&   component,
                                               const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_vrr_common_data_str())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_coords_compute(lines);
    
    for (const auto& label : _get_boys_vars_str(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_boys_compute_lines(lines, integral);
    
    lines.push_back({1, 0, 2, "// set up pointer to integrals buffer"});
    
    lines.push_back({1, 0, 2, "auto fints = buffer.data();"});
    
    lines.push_back({1, 0, 2, "// compute electron repulsion integrals"});
    
    _add_split_simd_code(lines, component, integral);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CFullPrimFuncBodyDriver::_get_common_data_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// set up math constants");
        
    vstr.push_back("const auto fpi = mathconst::getPiValue();");
    
    vstr.push_back("const auto invfpi = 1.0 / mathconst::getPiValue();");
    
    vstr.push_back("// set up coordinates for bra center A");

    vstr.push_back("const auto ra_x = coords_a[0];");

    vstr.push_back("const auto ra_y = coords_a[1];");

    vstr.push_back("const auto ra_z = coords_a[2];");

    vstr.push_back("// set up coordinates for bra center B");

    vstr.push_back("const auto rb_x = coords_b[0];");

    vstr.push_back("const auto rb_y = coords_b[1];");

    vstr.push_back("const auto rb_z = coords_b[2];");

    vstr.push_back("// set up coordinates for bra center C");

    vstr.push_back("const auto rc_x = coords_c_x.data();");

    vstr.push_back("const auto rc_y = coords_c_y.data();");

    vstr.push_back("const auto rc_z = coords_c_z.data();");

    vstr.push_back("// set up coordinates for bra center D");

    vstr.push_back("const auto rd_x = coords_d_x.data();");

    vstr.push_back("const auto rd_y = coords_d_y.data();");

    vstr.push_back("const auto rd_z = coords_d_z.data();");

    vstr.push_back("// set up ket side data");

    vstr.push_back("const auto fexps_c = ket_exps_c.data();");

    vstr.push_back("const auto fexps_d = ket_exps_d.data();");

    vstr.push_back("const auto knorms = ket_norms.data();");
    
    vstr.push_back("const auto kovls = ket_ovls.data();");

    vstr.push_back("// set up bra factors");

    vstr.push_back("const auto fe_ab_0 = bra_exp_a + bra_exp_b;");

    vstr.push_back("const auto fi_ab_0 = 1.0 / fe_ab_0;");

    vstr.push_back("// compute bra side overlap");

    vstr.push_back("const auto ab_x = ra_x - rb_x;");

    vstr.push_back("const auto ab_y = ra_y - rb_y;");

    vstr.push_back("const auto ab_z = ra_z - rb_z;");

    vstr.push_back("const auto fss_ab = bra_norm * bra_ovl;");
    
    return vstr;
}

std::vector<std::string>
T4CFullPrimFuncBodyDriver::_get_vrr_common_data_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// set up math constants");
        
    vstr.push_back("const auto fpi = mathconst::getPiValue();");
    
    vstr.push_back("const auto invfpi = 1.0 / mathconst::getPiValue();");
    
    vstr.push_back("// set up coordinates for bra center A");

    vstr.push_back("const auto ra_x = coords_a[0];");

    vstr.push_back("const auto ra_y = coords_a[1];");

    vstr.push_back("const auto ra_z = coords_a[2];");
    
    vstr.push_back("// set up coordinates for bra center B");

    vstr.push_back("const auto rb_x = coords_b[0];");

    vstr.push_back("const auto rb_y = coords_b[1];");

    vstr.push_back("const auto rb_z = coords_b[2];");
    
    vstr.push_back("// set up coordinates for bra center C");

    vstr.push_back("const auto rc_x = coords_c_x.data();");

    vstr.push_back("const auto rc_y = coords_c_y.data();");

    vstr.push_back("const auto rc_z = coords_c_z.data();");

    vstr.push_back("// set up coordinates for bra center D");

    vstr.push_back("const auto rd_x = coords_d_x.data();");

    vstr.push_back("const auto rd_y = coords_d_y.data();");

    vstr.push_back("const auto rd_z = coords_d_z.data();");

    vstr.push_back("// set up ket side data");

    vstr.push_back("const auto fexps_c = ket_exps_c.data();");

    vstr.push_back("const auto fexps_d = ket_exps_d.data();");

    vstr.push_back("const auto knorms = ket_norms.data();");
    
    vstr.push_back("const auto kovls = ket_ovls.data();");

    vstr.push_back("// set up bra factors");

    vstr.push_back("const auto fe_ab_0 = bra_exp_a + bra_exp_b;");

    vstr.push_back("const auto fi_ab_0 = 1.0 / fe_ab_0;");

    vstr.push_back("// compute bra side overlap");

    vstr.push_back("const auto fss_ab = bra_norm * bra_ovl;");
    
    return vstr;
}

void
T4CFullPrimFuncBodyDriver::_add_coords_compute(VCodeLines& lines) const
{
    lines.push_back({1, 0, 2, "// set up P center coordinates"});

    lines.push_back({1, 0, 2, "const auto rp_x = fi_ab_0 * (bra_exp_a * ra_x + bra_exp_b * rb_x);"});

    lines.push_back({1, 0, 2, "const auto rp_y = fi_ab_0 * (bra_exp_a * ra_y + bra_exp_b * rb_y);"});

    lines.push_back({1, 0, 2, "const auto rp_z = fi_ab_0 * (bra_exp_a * ra_z + bra_exp_b * rb_z);"});

    lines.push_back({1, 0, 2, "// compute Q center coordinates"});

    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray coords_q_x;"});
    
    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray coords_q_y;"});
    
    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray coords_q_z;"});

    lines.push_back({1, 0, 2, "auto rq_x = coords_q_x.data();"});

    lines.push_back({1, 0, 2, "auto rq_y = coords_q_y.data();"});

    lines.push_back({1, 0, 2, "auto rq_z = coords_q_z.data();"});

    lines.push_back({1, 0, 1, "#pragma omp simd aligned(rq_x, rq_y, rq_z, rc_x, rc_y, rc_z, rd_x, rd_y, rd_z, fexps_c, fexps_d : 64)"});
    
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < ket_dim; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "const auto fi_cd_0 = 1.0 / (fexps_c[i] + fexps_d[i]);"});

    lines.push_back({2, 0, 2, "rq_x[i] = fi_cd_0 * (fexps_c[i] * rc_x[i] + fexps_d[i] * rd_x[i]);"});

    lines.push_back({2, 0, 2, "rq_y[i] = fi_cd_0 * (fexps_c[i] * rc_y[i] + fexps_d[i] * rd_y[i]);"});

    lines.push_back({2, 0, 1, "rq_z[i] = fi_cd_0 * (fexps_c[i] * rc_z[i] + fexps_d[i] * rd_z[i]);"});
    
    lines.push_back({1, 0, 2, "}"});
}

std::vector<std::string>
T4CFullPrimFuncBodyDriver::_get_boys_vars_str(const I4CIntegral& integral) const
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
T4CFullPrimFuncBodyDriver::_add_boys_compute_lines(      VCodeLines&  lines,
                                                   const I4CIntegral& integral) const
{
    lines.push_back({1, 0, 2, "// compute Boys function and overlap values"});

    lines.push_back({1, 0, 1, "#pragma omp simd aligned(targs, fss_abcd, rc_x, rc_y, rc_z, rd_x, rd_y, rd_z, rq_x, rq_y, rq_z, fexps_c, fexps_d, knorms, kovls : 64)"});
    
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < ket_dim; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "const auto fe_cd_0 = fexps_c[i] + fexps_d[i];"});

    lines.push_back({2, 0, 2, "const auto fi_cd_0 = 1.0 / fe_cd_0;"});

    lines.push_back({2, 0, 2, "const auto rpq_x = rp_x - rq_x[i];"});

    lines.push_back({2, 0, 2, "const auto rpq_y = rp_y - rq_y[i];"});

    lines.push_back({2, 0, 2, "const auto rpq_z = rp_z - rq_z[i];"});

    lines.push_back({2, 0, 2, "targs[i] = fe_ab_0 * fe_cd_0 * (rpq_x * rpq_x + rpq_y * rpq_y + rpq_z * rpq_z) / (fe_ab_0 + fe_cd_0);"});

    //lines.push_back({2, 0, 2, "const auto cd_x = rc_x[i] - rd_x[i];"});

    //lines.push_back({2, 0, 2, "const auto cd_y = rc_y[i] - rd_y[i];"});

    //lines.push_back({2, 0, 2, "const auto cd_z = rc_z[i] - rd_z[i];"});

    lines.push_back({2, 0, 2, "fss_abcd[i] = 2.0 * fss_ab * knorms[i] * kovls[i] * std::sqrt(invfpi * fe_ab_0 * fe_cd_0 / (fe_ab_0 + fe_cd_0));"});
    
    lines.push_back({1, 0, 2, "}"});
    
    lines.push_back({1, 0, 2, "// rescale Boys function arguments and overlap for range sepatation"});
    
    lines.push_back({1, 0, 1, "if (use_rs)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 1, "#pragma omp simd aligned(targs, fss_abcd, fexps_c, fexps_d : 64)"});
    
    lines.push_back({2, 0, 1, "for (int64_t i = 0; i < ket_dim; i++)"});
    
    lines.push_back({2, 0, 1, "{"});
    
    lines.push_back({3, 0, 2, "const auto fe_cd_0 = fexps_c[i] + fexps_d[i];"});

    lines.push_back({3, 0, 2, "const auto frho = fe_ab_0 * fe_cd_0 / (fe_ab_0 + fe_cd_0);"});

    lines.push_back({3, 0, 2, "targs[i] *= omega * omega / (omega * omega + frho);"});

    lines.push_back({3, 0, 1, "fss_abcd[i] *= omega / std::sqrt(omega * omega + frho);"});
    
    lines.push_back({2, 0, 1, "}"});
    
    lines.push_back({1, 0, 2, "}"});
    
    const auto order = t4c::boys_order(integral);
    
    lines.push_back({1, 0, 2, "bf_table.compute<" + std::to_string(order + 1) + ">(bf_values, bf_args, ket_dim);"});
    
    lines.push_back({1, 0, 1, "if (use_rs)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 1, "#pragma omp simd aligned(fexps_c, fexps_d : 64)"});
    
    lines.push_back({2, 0, 1, "for (int64_t i = 0; i < ket_dim; i++)"});
    
    lines.push_back({2, 0, 1, "{"});
    
    lines.push_back({3, 0, 2, "const auto fe_cd_0 = fexps_c[i] + fexps_d[i];"});

    lines.push_back({3, 0, 2, "auto frho = fe_ab_0 * fe_cd_0 / (fe_ab_0 + fe_cd_0);"});
    
    lines.push_back({3, 0, 2, "const auto fact = omega * omega / (omega * omega + frho);"});
    
    for (size_t i = 1; i <= order; i++)
    {
        if (i == 1)
        {
            lines.push_back({3, 0, 2, "b" + std::to_string(i) + "_vals[i] *= fact;"});
        }
        else if (i == 2)
        {
            lines.push_back({3, 0, 2, "frho = fact * fact;"});
            
            lines.push_back({3, 0, 2, "b" + std::to_string(i) + "_vals[i] *= frho;"});
        }
        else
        {
            lines.push_back({3, 0, 2, "frho *= fact;"});
            
            lines.push_back({3, 0, 2, "b" + std::to_string(i) + "_vals[i] *= frho;"});
        }
    }
    
    lines.push_back({2, 0, 1, "}"});
    
    lines.push_back({1, 0, 2, "}"});
}

void
T4CFullPrimFuncBodyDriver::_add_split_simd_code(      VCodeLines&  lines,
                                                const T4CIntegral& component,
                                                const I4CIntegral& integral) const
{
    T4CFullElectronRepulsionDriver t4c_eri_drv;
    
    const auto rdist = (t4c_eri_drv.create_recursion({component, }))[0];
        
    for (const auto& tint : rdist.unique_integrals())
    {
        const auto tdist = rdist.split(tint);
            
        _add_split_simd_block(lines, tint, tdist);
    }
}

void
T4CFullPrimFuncBodyDriver::_add_split_simd_block(      VCodeLines&  lines,
                                                 const T4CIntegral& integral,
                                                 const R4CDist&     rdist) const
{
    const auto order = integral.order();
    
    lines.push_back({1, 0, 2, "// add Boys order " + std::to_string(order) + " contributions"});
    
    _add_split_pragma(lines, integral, rdist);
    
    _add_split_loop_start(lines, integral, rdist);
    
    _add_simd_lines_block(lines, integral, rdist, 0);
    
    _add_split_loop_end(lines, integral);
}

void
T4CFullPrimFuncBodyDriver::_add_split_pragma(      VCodeLines&  lines,
                                             const T4CIntegral& integral,
                                             const R4CDist&     rdist) const
{
    const auto order = integral.order();
    
    std::string vars_str;
    
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
        vars_str += " rq_x,";
    }
    
    if (t4c::find_factor(rdist, "rwp_y") || t4c::find_factor(rdist, "rwq_y"))
    {
        vars_str += " rq_y,";
    }
    
    if (t4c::find_factor(rdist, "rwp_z") || t4c::find_factor(rdist, "rwq_z"))
    {
        vars_str += " rq_z,";
    }
    
    if (t4c::find_factor(rdist, "rqc_x"))
    {
        vars_str += " rc_x,";
    }
    
    if (t4c::find_factor(rdist, "rqc_y"))
    {
        vars_str += " rc_y,";
    }
    
    if (t4c::find_factor(rdist, "rqc_z"))
    {
        vars_str += " rc_z,";
    }
    
    if (t4c::find_factor(rdist, "rqd_x"))
    {
        vars_str += " rd_x,";
    }
    
    if (t4c::find_factor(rdist, "rqd_y"))
    {
        vars_str += " rd_y,";
    }
    
    if (t4c::find_factor(rdist, "rqd_z"))
    {
        vars_str += " rd_z,";
    }
    
    //vars_str += " b" + std::to_string(order) + "_vals";
    
    lines.push_back({1, 0, 1, "#pragma omp simd aligned(fints, " + vars_str + " fss_abcd : 64)"});
}

void
T4CFullPrimFuncBodyDriver::_add_split_loop_start(      VCodeLines&  lines,
                                                 const T4CIntegral& integral,
                                                 const R4CDist&     rdist) const
{
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < ket_dim; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    
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
        lines.push_back({2, 0, 2, "const auto rwp_x = fe_cd_0 * fi_abcd_0 * (rq_x[i] - rp_x);"});
    }
    
    if (t4c::find_factor(rdist, "rwp_y"))
    {
        lines.push_back({2, 0, 2, "const auto rwp_y = fe_cd_0 * fi_abcd_0 * (rq_y[i] - rp_y);"});
    }
    
    if (t4c::find_factor(rdist, "rwp_z"))
    {
        lines.push_back({2, 0, 2, "const auto rwp_z = fe_cd_0 * fi_abcd_0 * (rq_z[i] - rp_z);"});
    }
    
    if (t4c::find_factor(rdist, "rwq_x"))
    {
        lines.push_back({2, 0, 2, "const auto rwq_x = fe_ab_0 * fi_abcd_0 * (rp_x - rq_x[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rwq_y"))
    {
        lines.push_back({2, 0, 2, "const auto rwq_y = fe_ab_0 * fi_abcd_0 * (rp_y - rq_y[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rwq_z"))
    {
        lines.push_back({2, 0, 2, "const auto rwq_z = fe_ab_0 * fi_abcd_0 * (rp_z - rq_z[i]);"});
    }
    
    if (t4c::find_factor(rdist, "rpa_x"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_x = rp_x - ra_x;"});
    }
    
    if (t4c::find_factor(rdist, "rpa_y"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_y = rp_y - ra_y;"});
    }
    
    if (t4c::find_factor(rdist, "rpa_z"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_z = rp_z - ra_z;"});
    }
    
    if (t4c::find_factor(rdist, "rpb_x"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_x = rp_x - rb_x;"});
    }
    
    if (t4c::find_factor(rdist, "rpb_y"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_y = rp_y - rb_y;"});
    }
    
    if (t4c::find_factor(rdist, "rpb_z"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_z = rp_z - rb_z;"});
    }
    
    if (t4c::find_factor(rdist, "rqc_x"))
    {
        lines.push_back({2, 0, 2, "const auto rqc_x = rq_x[i] - rc_x[i];"});
    }
    
    if (t4c::find_factor(rdist, "rqc_y"))
    {
        lines.push_back({2, 0, 2, "const auto rqc_y = rq_y[i] - rc_y[i];"});
    }
    
    if (t4c::find_factor(rdist, "rqc_z"))
    {
        lines.push_back({2, 0, 2, "const auto rqc_z = rq_z[i] - rc_z[i];"});
    }
    
    if (t4c::find_factor(rdist, "rqd_x"))
    {
        lines.push_back({2, 0, 2, "const auto rqd_x = rq_x[i] - rd_x[i];"});
    }
    
    if (t4c::find_factor(rdist, "rqd_y"))
    {
        lines.push_back({2, 0, 2, "const auto rqd_y = rq_y[i] - rd_y[i];"});
    }
    
    if (t4c::find_factor(rdist, "rqd_z"))
    {
        lines.push_back({2, 0, 2, "const auto rqd_z = rq_z[i] - rd_z[i];"});
    }
}

void
T4CFullPrimFuncBodyDriver::_add_split_loop_end(      VCodeLines&  lines,
                                               const T4CIntegral& integral) const
{
    const auto order = integral.order();
    
    lines.push_back({2, 0, 1, "fints[i] += fss_abcd[i] * fact * b" + std::to_string(order) + "_vals[i];"});
    
    lines.push_back({1, 0, 2, "}"});
}

void
T4CFullPrimFuncBodyDriver::_add_simd_lines_block(      VCodeLines&  lines,
                                                 const T4CIntegral& integral,
                                                 const R4CDist&     rdist,
                                                 const size_t       index) const
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
            simd_str += t4c::get_factor_label(rdist[j], integral, j == sterm, false);
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
