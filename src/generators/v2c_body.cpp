#include "v2c_body.hpp"

#include "string_formater.hpp"
#include "angular_components.hpp"
#include "v2i_ovl_driver.hpp"

void
V2CFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                   const SI2CIntegrals& integrals,
                                   const I2CIntegral&   integral,
                                   const bool           sum_form,
                                   const bool           diagonal) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_gtos_def(diagonal))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_ket_variables_def(diagonal))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_coordinates_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_buffers_def(integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_loop_start(lines, integral);
    
    _add_ket_loop_start(lines, integral, diagonal);
    
    _add_auxilary_integrals(lines, integrals);
    
    _add_call_tree(lines, integrals);
    
    _add_ket_loop_end(lines, integral, diagonal);
    
    _add_loop_end(lines, integral, diagonal);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
V2CFuncBodyDriver::_get_gtos_def(const bool diagonal) const
{
    std::vector<std::string> vstr;
    
    if (diagonal)
    {
        vstr.push_back("// intialize GTOs data");

        vstr.push_back("const auto gto_coords_x = gto_block.coordinates_x();");

        vstr.push_back("const auto gto_coords_y = gto_block.coordinates_y();");

        vstr.push_back("const auto gto_coords_z = gto_block.coordinates_z();");

        vstr.push_back("const auto gto_exps = gto_block.exponents();");

        vstr.push_back("const auto gto_norms = gto_block.normalization_factors();");

        vstr.push_back("const auto gto_indices = gto_block.orbital_indices();");

        vstr.push_back("const auto ncgtos = gto_block.number_of_basis_functions();");

        vstr.push_back("const auto npgtos = gto_block.number_of_primitives();");
    }
    else
    {
        vstr.push_back("// intialize GTOs data on bra side");
        
        vstr.push_back("const auto bra_gto_coords_x = bra_gto_block.coordinates_x();");
        
        vstr.push_back("const auto bra_gto_coords_y = bra_gto_block.coordinates_y();");
        
        vstr.push_back("const auto bra_gto_coords_z = bra_gto_block.coordinates_z();");
        
        vstr.push_back("const auto bra_gto_exps = bra_gto_block.exponents();");
        
        vstr.push_back("const auto bra_gto_norms = bra_gto_block.normalization_factors();");
       
        vstr.push_back("const auto bra_gto_indices = bra_gto_block.orbital_indices();");
        
        vstr.push_back("const auto bra_ncgtos = bra_gto_block.number_of_basis_functions();");

        vstr.push_back("const auto bra_npgtos = bra_gto_block.number_of_primitives();");
        
        vstr.push_back("// intialize GTOs data on ket side");
        
        vstr.push_back("const auto ket_gto_coords_x = ket_gto_block.coordinates_x();");
        
        vstr.push_back("const auto ket_gto_coords_y = ket_gto_block.coordinates_y();");
        
        vstr.push_back("const auto ket_gto_coords_z = ket_gto_block.coordinates_z();");
        
        vstr.push_back("const auto ket_gto_exps = ket_gto_block.exponents();");
        
        vstr.push_back("const auto ket_gto_norms = ket_gto_block.normalization_factors();");
       
        vstr.push_back("const auto ket_gto_indices = ket_gto_block.orbital_indices();");
        
        vstr.push_back("const auto ket_ncgtos = ket_gto_block.number_of_basis_functions();");

        vstr.push_back("const auto ket_npgtos = ket_gto_block.number_of_primitives();");
    }
    
    return vstr;
}

std::vector<std::string>
V2CFuncBodyDriver::_get_ket_variables_def(const bool diagonal) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned 2D arrays for ket side");
    
    vstr.push_back("const auto ket_dim = ket_indices[1] - ket_indices[0];");
    
    if (diagonal)
    {
        vstr.push_back("const auto ket_pdim = ket_dim * npgtos;");
    }
    else
    {
        vstr.push_back("const auto ket_pdim = ket_dim * ket_npgtos;");
    }
    
    vstr.push_back("CSimdArray<double> b_x(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> b_y(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> b_z(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> b_exps(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> b_norms(1, ket_pdim);");

    vstr.push_back(" // load GTOs data for ket side");

    if (diagonal)
    {
        vstr.push_back("b_x.replicate(gto_coords_x, ket_indices, npgtos);");

        vstr.push_back("b_y.replicate(gto_coords_y, ket_indices, npgtos);");

        vstr.push_back("b_z.replicate(gto_coords_z, ket_indices, npgtos);");

        vstr.push_back("b_exps.load(gto_exps, ket_indices, npgtos);");

        vstr.push_back("b_norms.load(gto_norms, ket_indices, npgtos);");
    }
    else
    {
        vstr.push_back("b_x.replicate(ket_gto_coords_x, ket_indices, ket_npgtos);");

        vstr.push_back("b_y.replicate(ket_gto_coords_y, ket_indices, ket_npgtos);");

        vstr.push_back("b_z.replicate(ket_gto_coords_z, ket_indices, ket_npgtos);");

        vstr.push_back("b_exps.load(ket_gto_exps, ket_indices, ket_npgtos);");

        vstr.push_back("b_norms.load(ket_gto_norms, ket_indices, ket_npgtos);");
    }
    
    return vstr;
}

std::vector<std::string>
V2CFuncBodyDriver::_get_coordinates_def(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned distances R(AB) = A - B");

    vstr.push_back("CSimdArray<double> ab_x(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> ab_y(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> ab_z(1, ket_pdim);");
    
    if (integral[1] > 0)
    {
        vstr.push_back("// allocate aligned distances R(PB) = P - B");
        
        vstr.push_back("CSimdArray<double> pb_x(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> pb_y(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> pb_z(1, ket_pdim);");
    }
    
    if (integral[0] > 0)
    {
        vstr.push_back("// allocate aligned distances R(PA) = P - A");
        
        vstr.push_back("CSimdArray<double> pa_x(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> pa_y(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> pa_z(1, ket_pdim);");
    }
    
    return vstr;
}

std::vector<std::string>
V2CFuncBodyDriver::_get_buffers_def(const SI2CIntegrals& integrals,
                                    const I2CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned primitive integrals");
    
    for (const auto& tint : integrals)
    {
        std::string label = "CSimdArray<double> ";
        
        label += t2c::get_buffer_label(tint, {"prim"});
        
        const auto angpair = std::array<int, 2>({tint[0], tint[1]});
        
        const auto tcomps = ten::number_of_cartesian_components(angpair);
        
        label += "(" + std::to_string(tcomps) + ", ket_pdim);";
        
        vstr.push_back(label);
    }
    
    vstr.push_back("// allocate aligned contracted integrals");
    
    const auto angpair = std::array<int, 2>({integral[0], integral[1]});
    
    auto icomps = ten::number_of_cartesian_components(angpair);
    
    std::string label = "CSimdArray<double> ";
    
    label += t2c::get_buffer_label(integral, {"cart"});
    
    label += "(" + std::to_string(icomps) + ", ket_dim);";
    
    vstr.push_back(label);
    
    if ((integral[0] > 1) || (integral[1] > 1))
    {
        icomps = ten::number_of_spherical_components(angpair);
        
        std::string label = "CSimdArray<double> ";
        
        label += t2c::get_buffer_label(integral, {"spher"});
        
        label += "(" + std::to_string(icomps) + ", ket_dim);";
        
        vstr.push_back(label);
    }
    
    return vstr;
}

void
V2CFuncBodyDriver::_add_loop_start(      VCodeLines&  lines,
                                   const I2CIntegral& integral) const
{
    lines.push_back({1, 0, 2, "// loop over contracted GTOs on bra side"});
        
    lines.push_back({1, 0, 1, "for (auto i = bra_indices[0]; i < bra_indices[1]; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, t2c::get_buffer_label(integral, {"cart"}) +".zero();"});

    if ((integral[0] > 1) || (integral[1] > 1))
    {
        lines.push_back({2, 0, 2, t2c::get_buffer_label(integral, {"spher"}) +".zero();"});
    }

    lines.push_back({2, 0, 2, "const auto a_x = gto_coords_x[i];"});

    lines.push_back({2, 0, 2, "const auto a_y = gto_coords_y[i];"});

    lines.push_back({2, 0, 2, "const auto a_z = gto_coords_z[i];"});

    lines.push_back({2, 0, 2, "t2cfunc::comp_distances_ab(ab_x[0], ab_y[0], ab_z[0], a_x, a_y, a_z, b_x[0], b_y[0], b_z[0], ket_pdim);"});
}

void
V2CFuncBodyDriver::_add_ket_loop_start(      VCodeLines&  lines,
                                       const I2CIntegral& integral,
                                       const bool         diagonal) const
{
    if (diagonal)
    {
        lines.push_back({2, 0, 1, "for (int j = 0; j < npgtos; j++))"});
    }
    else
    {
        lines.push_back({2, 0, 1, "for (int j = 0; j < bra_npgtos; j++)"});
    }
    
    lines.push_back({2, 0, 1, "{"});
    
    if (diagonal)
    {
        lines.push_back({3, 0, 2, "const auto a_exp = gto_exps[j * ncgtos + i];"});
            
        lines.push_back({3, 0, 2, "const auto a_norm = gto_norms[j * ncgtos + i];"});
    }
    else
    {
        lines.push_back({3, 0, 2, "const auto a_exp = bra_gto_exps[j * bra_ncgtos + i];"});
            
        lines.push_back({3, 0, 2, "const auto a_norm = bra_gto_norms[j * bra_ncgtos + i];"});
    }
    
    if (integral[0] > 1)
    {
        lines.push_back({3, 0, 2, "t2cfunc::comp_distances_pa(pa_x[0], pa_y[0], pa_z[0], ab_x[0], ab_y[0], ab_z[0], a_exp, b_exps[0], ket_pdim);"});
    }

    if (integral[1] > 1)
    {
        lines.push_back({3, 0, 2, "t2cfunc::comp_distances_pb(pb_x[0], pb_y[0], pb_z[0], ab_x[0], ab_y[0], ab_z[0], a_exp, b_exps[0], ket_pdim);"});
    }
}

void
V2CFuncBodyDriver::_add_auxilary_integrals(      VCodeLines&  lines,
                                           const SI2CIntegrals& integrals) const
{
    for (const auto& tint : integrals)
    {
        if ((tint[0] == 0) && (tint[1] == 0))
        {
            if (tint.integrand().name() == "1")
            {
                lines.push_back({3, 0, 2, "ovlrec::comp_prim_overlap_s_s(prim_buffer_ovl_ss, ab_x[0], ab_y[0], ab_z[0], a_exp, b_exps[0], a_norm, b_norms[0]);"});
            }
            
            // TODO: other integrals...
        }
    }
}

void
V2CFuncBodyDriver::_add_call_tree(      VCodeLines&  lines,
                                  const SI2CIntegrals& integrals) const
{
    for (const auto& tint : integrals)
    {
        if ((tint[0] != 0) || (tint[1] != 0))
        {
            const auto [nsize, name] = t2c::prim_compute_func_name(tint, false);
            
            auto label = t2c::namespace_label(tint) + "::" + name + "(";
            
            label += _get_arguments(tint);
            
            if (tint[0] > 0)
            {
                label += "pa_x[0], pa_y[0], pa_z[0], ";
            }
            
            if ((tint[1] > 0) && (tint[0] == 0))
            {
                label += "pb_x[0], pb_y[0], pb_z[0], ";
            }
            
            if ((tint[0] + tint[1]) > 1)
            {
                label += "a_exp, b_exps[0]";
            }
            
            label += ");";
            
            lines.push_back({3, 0, 2, label});
        }
    }
}

std::string
V2CFuncBodyDriver::_get_arguments(const I2CIntegral& integral) const
{
    std::string label;
    
    
    for (const auto& tint : t2c::get_integrals(integral))
    {
        if (tint != integral)
        {
            label += t2c::get_buffer_label(tint, {"prim"}) + ", ";
        }
    }
    
    return label;
}

void
V2CFuncBodyDriver::_add_ket_loop_end(      VCodeLines&  lines,
                                     const I2CIntegral& integral,
                                     const bool         diagonal) const
{
    std::string label = "t2cfunc::reduce(";
    
    label += t2c::get_buffer_label(integral, "cart") + ", ";
    
    label += t2c::get_buffer_label(integral, "prim") + ", ";
    
    if (diagonal)
    {
        label += "ket_dim, npgtos);";
    }
    else
    {
        label += "ket_dim, ket_npgtos);";
    }
    
    lines.push_back({3, 0, 1, label});
    
    lines.push_back({2, 0, 2, "}"});
}

void
V2CFuncBodyDriver::_add_loop_end(      VCodeLines&  lines,
                                 const I2CIntegral& integral,
                                 const bool         diagonal) const
{
    std::string label;
    
    if ((integral[0] > 1) || (integral[1] > 1))
    {
        label = "t2cfunc::transform<"  + std::to_string(integral[0]);
        
        label += ", " + std::to_string(integral[1]) + ">(";
        
        label += t2c::get_buffer_label(integral, "spher") + ", ";
        
        label += t2c::get_buffer_label(integral, "cart") + ");";
        
        lines.push_back({2, 0, 2, label});
        
        label = "t2cfunc::distribute(matrix," + t2c::get_buffer_label(integral, "spher") + ", ";
       
        if (diagonal)
        {
            label += "gto_indices, ";
        }
        else
        {
            label += "bra_gto_indices, ket_gto_indices, ";
        }
        
        label += std::to_string(integral[0]) + ", ";
        
        label += std::to_string(integral[1]) + ", ";
        
        if (diagonal)
        {
            label += "i, ket_indices);";
        }
        else
        {
            if (integral[0] == integral[1])
            {
                label += "i, ket_indices, mat_type);";
            }
            else
            {
                label += "i, ket_indices, ang_order);";
            }
        }
        
        lines.push_back({2, 0, 1, label});
    }
    
    label = t2c::get_buffer_label(integral, "cart");
    
    if ((integral[0] == 0) && (integral[1] == 0))
    {
        if (diagonal)
        {
            lines.push_back({2, 0, 1, "t2cfunc::distribute(matrix, " + label + "[0], gto_indices, 0, 0, i, ket_indices);"});
        }
        else
        {
            lines.push_back({2, 0, 1, "t2cfunc::distribute(matrix, " + label  + "[0], bra_gto_indices, ket_gto_indices, 0, 0, i, ket_indices, mat_type);"});
        }
    }
    
    if ((integral[0] == 0) && (integral[1] == 1))
    {
        lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label + "[0], bra_gto_indices, ket_gto_indices, 0, 2, i, ket_indices, ang_order);"});
        
        lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label + "[1], bra_gto_indices, ket_gto_indices, 0, 0, i, ket_indices, ang_order);"});

        lines.push_back({2, 0, 1, "t2cfunc::distribute(matrix, " + label + "[2], bra_gto_indices, ket_gto_indices, 0, 1, i, ket_indices, ang_order);"});
    }
    
    if ((integral[0] == 1) && (integral[1] == 1))
    {
        lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label + "[0], bra_gto_indices, ket_gto_indices, 2, 0, i, ket_indices, ang_order);"});
        
        lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label + "[1], bra_gto_indices, ket_gto_indices, 0, 0, i, ket_indices, ang_order);"});

        lines.push_back({2, 0, 1, "t2cfunc::distribute(matrix, " + label + "[2], bra_gto_indices, ket_gto_indices, 1, 0, i, ket_indices, ang_order);"});
    }
    
    if ((integral[0] == 1) && (integral[1] == 1))
    {
        if (diagonal)
        {
            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[0], gto_indices, 2, 2, i, ket_indices);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[1], gto_indices, 2, 0, i, ket_indices);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[2], gto_indices, 2, 1, i, ket_indices);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[3], gto_indices, 0, 2, i, ket_indices);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[4], gto_indices, 0, 0, i, ket_indices);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[5], gto_indices, 0, 1, i, ket_indices);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[6], gto_indices, 1, 2, i, ket_indices);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[7], gto_indices, 1, 0, i, ket_indices);"});

            lines.push_back({2, 0, 1, "t2cfunc::distribute(matrix, " + label  + "[8], gto_indices, 1, 1, i, ket_indices);"});
        }
        else
        {
            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[0], bra_gto_indices, ket_gto_indices, 2, 2, i, ket_indices, mat_type);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[1], bra_gto_indices, ket_gto_indices, 2, 0, i, ket_indices, mat_type);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[2], bra_gto_indices, ket_gto_indices, 2, 1, i, ket_indices, mat_type);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[3], bra_gto_indices, ket_gto_indices, 0, 2, i, ket_indices, mat_type);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[4], bra_gto_indices, ket_gto_indices, 0, 0, i, ket_indices, mat_type);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[5], bra_gto_indices, ket_gto_indices, 0, 1, i, ket_indices, mat_type);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[6], bra_gto_indices, ket_gto_indices, 1, 2, i, ket_indices, mat_type);"});

            lines.push_back({2, 0, 2, "t2cfunc::distribute(matrix, " + label  + "[7], bra_gto_indices, ket_gto_indices, 1, 0, i, ket_indices, mat_type);"});

            lines.push_back({2, 0, 1, "t2cfunc::distribute(matrix, " + label  + "[8], bra_gto_indices, ket_gto_indices, 1, 1, i, ket_indices, mat_type);"});
        }
    }
    
    lines.push_back({1, 0, 1, "}"});
}
