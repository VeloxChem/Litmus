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

#include "t4c_body.hpp"

#include "t2c_utils.hpp"
#include "t4c_utils.hpp"

void
T4CFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                   const SI4CIntegrals& bra_integrals,
                                   const SI4CIntegrals& ket_integrals,
                                   const SI4CIntegrals& vrr_integrals,
                                   const I4CIntegral&   integral,
                                   const bool           diagonal) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_gto_pairs_def(diagonal))
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
    
    for (const auto& label : _get_prim_buffers_def(vrr_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_cart_buffers_def(bra_integrals, ket_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_contr_buffers_def(bra_integrals, ket_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_half_spher_buffers_def(bra_integrals, ket_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_spher_buffers_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_boys_function_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_loop_start(lines, bra_integrals, ket_integrals, integral, diagonal);
    
    _add_ket_loop_start(lines, integral, diagonal);

    _add_auxilary_integrals(lines, vrr_integrals);

    _add_vrr_call_tree(lines, vrr_integrals);

    _add_ket_loop_end(lines, bra_integrals, ket_integrals, integral, diagonal);
    
    _add_ket_hrr_call_tree(lines, ket_integrals);
    
    _add_ket_trafo_call_tree(lines, bra_integrals, ket_integrals, integral);
    
    _add_bra_hrr_call_tree(lines, bra_integrals);
    
    _add_bra_trafo_call_tree(lines, integral);
    
    _add_loop_end(lines, integral, diagonal);
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

void
T4CFuncBodyDriver::write_diag_func_body(      std::ofstream& fstream,
                                        const SI4CIntegrals& bra_integrals,
                                        const SI4CIntegrals& ket_integrals,
                                        const SI4CIntegrals& vrr_integrals,
                                        const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_gto_pairs_def(true))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_diag_ket_variables_def())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_diag_coordinates_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_diag_prim_buffers_def(vrr_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_diag_cart_buffers_def(bra_integrals, ket_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_diag_contr_buffers_def(bra_integrals, ket_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_diag_half_spher_buffers_def(bra_integrals, ket_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_diag_spher_buffers_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_diag_boys_function_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
//    
//    _add_loop_start(lines, bra_integrals, ket_integrals, integral, diagonal);
//    
//    _add_ket_loop_start(lines, integral, diagonal);
//
//    _add_auxilary_integrals(lines, vrr_integrals);
//
//    _add_vrr_call_tree(lines, vrr_integrals);
//
//    _add_ket_loop_end(lines, bra_integrals, ket_integrals, integral, diagonal);
//    
//    _add_ket_hrr_call_tree(lines, ket_integrals);
//    
//    _add_ket_trafo_call_tree(lines, bra_integrals, ket_integrals, integral);
//    
//    _add_bra_hrr_call_tree(lines, bra_integrals);
//    
//    _add_bra_trafo_call_tree(lines, integral);
//    
//    _add_loop_end(lines, integral, diagonal);
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

void
T4CFuncBodyDriver::write_geom_func_body(      std::ofstream& fstream,
                                        const SI4CIntegrals& geom_integrals,
                                        const SI4CIntegrals& vrr_integrals,
                                        const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_gto_pairs_def(false))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_ket_variables_def(false))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_full_coordinates_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_full_prim_buffers_def(vrr_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_full_cart_buffers_def(geom_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_full_spher_buffers_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_boys_function_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_full_loop_start(lines, geom_integrals, integral);
    
    _add_full_ket_loop_start(lines, integral);

    _add_auxilary_integrals(lines, vrr_integrals);

    _add_full_vrr_call_tree(lines, vrr_integrals);
    
    _add_geom_call_tree(lines, geom_integrals, integral);

    _add_full_ket_loop_end(lines, integral);

    _add_full_trafo(lines, integral);

    _add_full_loop_end(lines, integral);
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CFuncBodyDriver::_get_gto_pairs_def(const bool diagonal) const
{
    std::vector<std::string> vstr;
    
    if (diagonal)
    {
        vstr.push_back("// intialize GTOs pair data");

        vstr.push_back("const auto a_coords_x = gto_pair_block.bra_coordinates_x();");
        
        vstr.push_back("const auto a_coords_y = gto_pair_block.bra_coordinates_y();");
        
        vstr.push_back("const auto a_coords_z = gto_pair_block.bra_coordinates_z();");
        
        vstr.push_back("const auto b_coords_x = gto_pair_block.ket_coordinates_x();");
        
        vstr.push_back("const auto b_coords_y = gto_pair_block.ket_coordinates_y();");
        
        vstr.push_back("const auto b_coords_z = gto_pair_block.ket_coordinates_z();");
        
        vstr.push_back("const auto a_vec_exps = gto_pair_block.bra_exponents();");
        
        vstr.push_back("const auto b_vec_exps = gto_pair_block.ket_exponents();");
        
        vstr.push_back("const auto ab_vec_norms = gto_pair_block.normalization_factors();");
        
        vstr.push_back("const auto ab_vec_ovls = gto_pair_block.overlap_factors();");
        
        vstr.push_back("const auto a_indices = gto_pair_block.bra_orbital_indices();");
        
        vstr.push_back("const auto b_indices = gto_pair_block.ket_orbital_indices();");
        
        vstr.push_back("const auto ncgtos = gto_pair_block.number_of_contracted_pairs();");
        
        vstr.push_back("const auto npgtos = gto_pair_block.number_of_primitive_pairs();");
    }
    else
    {
        vstr.push_back("// intialize GTOs pair data on bra side");

        vstr.push_back("const auto a_coords_x = bra_gto_pair_block.bra_coordinates_x();");
        
        vstr.push_back("const auto a_coords_y = bra_gto_pair_block.bra_coordinates_y();");
        
        vstr.push_back("const auto a_coords_z = bra_gto_pair_block.bra_coordinates_z();");
        
        vstr.push_back("const auto b_coords_x = bra_gto_pair_block.ket_coordinates_x();");
        
        vstr.push_back("const auto b_coords_y = bra_gto_pair_block.ket_coordinates_y();");
        
        vstr.push_back("const auto b_coords_z = bra_gto_pair_block.ket_coordinates_z();");
        
        vstr.push_back("const auto a_vec_exps = bra_gto_pair_block.bra_exponents();");
        
        vstr.push_back("const auto b_vec_exps = bra_gto_pair_block.ket_exponents();");
        
        vstr.push_back("const auto ab_vec_norms = bra_gto_pair_block.normalization_factors();");
        
        vstr.push_back("const auto ab_vec_ovls = bra_gto_pair_block.overlap_factors();");
        
        vstr.push_back("const auto a_indices = bra_gto_pair_block.bra_orbital_indices();");
        
        vstr.push_back("const auto b_indices = bra_gto_pair_block.ket_orbital_indices();");
        
        //vstr.push_back("const auto bra_ang_mom = bra_gto_pair_block.angular_momentums();");
        
        vstr.push_back("const auto bra_ncgtos = bra_gto_pair_block.number_of_contracted_pairs();");
        
        vstr.push_back("const auto bra_npgtos = bra_gto_pair_block.number_of_primitive_pairs();");
        
        vstr.push_back("// intialize GTOs data on ket side");
        
        vstr.push_back("const auto c_coords_x = ket_gto_pair_block.bra_coordinates_x();");
        
        vstr.push_back("const auto c_coords_y = ket_gto_pair_block.bra_coordinates_y();");
        
        vstr.push_back("const auto c_coords_z = ket_gto_pair_block.bra_coordinates_z();");
        
        vstr.push_back("const auto d_coords_x = ket_gto_pair_block.ket_coordinates_x();");
        
        vstr.push_back("const auto d_coords_y = ket_gto_pair_block.ket_coordinates_y();");
        
        vstr.push_back("const auto d_coords_z = ket_gto_pair_block.ket_coordinates_z();");
        
        vstr.push_back("const auto c_vec_exps = ket_gto_pair_block.bra_exponents();");
        
        vstr.push_back("const auto d_vec_exps = ket_gto_pair_block.ket_exponents();");
        
        vstr.push_back("const auto cd_vec_norms = ket_gto_pair_block.normalization_factors();");
        
        vstr.push_back("const auto cd_vec_ovls = ket_gto_pair_block.overlap_factors();");
        
        vstr.push_back("const auto c_indices = ket_gto_pair_block.bra_orbital_indices();");
        
        vstr.push_back("const auto d_indices = ket_gto_pair_block.ket_orbital_indices();");
        
        vstr.push_back("const auto ket_npgtos = ket_gto_pair_block.number_of_primitive_pairs();");
    }
    
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_ket_variables_def(const bool diagonal) const
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
    
    vstr.push_back("CSimdArray<double> c_x(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> c_y(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> c_z(1, ket_pdim);");
    
    vstr.push_back("CSimdArray<double> d_x(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> d_y(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> d_z(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> c_exps(1, ket_pdim);");
    
    vstr.push_back("CSimdArray<double> d_exps(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> cd_norms(1, ket_pdim);");
    
    vstr.push_back("CSimdArray<double> cd_ovls(1, ket_pdim);");

    vstr.push_back(" // load GTOs data for ket side");

    if (diagonal)
    {
        vstr.push_back("c_x.replicate(a_coords_x, ket_indices, npgtos);");

        vstr.push_back("c_y.replicate(a_coords_y, ket_indices, npgtos);");

        vstr.push_back("c_z.replicate(a_coords_z, ket_indices, npgtos);");
        
        vstr.push_back("d_x.replicate(b_coords_x, ket_indices, npgtos);");

        vstr.push_back("d_y.replicate(b_coords_y, ket_indices, npgtos);");

        vstr.push_back("d_z.replicate(b_coords_z, ket_indices, npgtos);");

        vstr.push_back("c_exps.load(a_vec_exps, ket_indices, npgtos);");
        
        vstr.push_back("d_exps.load(b_vec_exps, ket_indices, npgtos);");

        vstr.push_back("cd_norms.load(ab_vec_norms, ket_indices, npgtos);");
        
        vstr.push_back("cd_ovls.load(ab_vec_ovls, ket_indices, npgtos);");
    }
    else
    {
        vstr.push_back("c_x.replicate(c_coords_x, ket_indices, ket_npgtos);");

        vstr.push_back("c_y.replicate(c_coords_y, ket_indices, ket_npgtos);");

        vstr.push_back("c_z.replicate(c_coords_z, ket_indices, ket_npgtos);");
        
        vstr.push_back("d_x.replicate(d_coords_x, ket_indices, ket_npgtos);");

        vstr.push_back("d_y.replicate(d_coords_y, ket_indices, ket_npgtos);");

        vstr.push_back("d_z.replicate(d_coords_z, ket_indices, ket_npgtos);");

        vstr.push_back("c_exps.load(c_vec_exps, ket_indices, ket_npgtos);");
        
        vstr.push_back("d_exps.load(d_vec_exps, ket_indices, ket_npgtos);");

        vstr.push_back("cd_norms.load(cd_vec_norms, ket_indices, ket_npgtos);");
        
        vstr.push_back("cd_ovls.load(cd_vec_ovls, ket_indices, ket_npgtos);");
    }
    
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_diag_ket_variables_def() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned 2D arrays for ket side");
    
    vstr.push_back("CSimdArray<double> c_x(1, npgtos);");

    vstr.push_back("CSimdArray<double> c_y(1, npgtos);");

    vstr.push_back("CSimdArray<double> c_z(1, npgtos);");
    
    vstr.push_back("CSimdArray<double> d_x(1, npgtos);");

    vstr.push_back("CSimdArray<double> d_y(1, npgtos);");

    vstr.push_back("CSimdArray<double> d_z(1, npgtos);");

    vstr.push_back("CSimdArray<double> c_exps(1, npgtos);");
    
    vstr.push_back("CSimdArray<double> d_exps(1, npgtos);");

    vstr.push_back("CSimdArray<double> cd_norms(1, npgtos);");
    
    vstr.push_back("CSimdArray<double> cd_ovls(1, npgtos);");

    return vstr;
}


std::vector<std::string>
T4CFuncBodyDriver::_get_coordinates_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    const auto integrand = integral.integrand();
    
    vstr.push_back("// allocate aligned coordinates of Q center");

    vstr.push_back("CSimdArray<double> q_x(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> q_y(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> q_z(1, ket_pdim);");
        
    if ((integral[0] + integral[1] + integral[2] + integral[3]) > 0)
    {
        vstr.push_back("// allocate aligned coordinates of W center");

        vstr.push_back("CSimdArray<double> w_x(1, ket_pdim);");

        vstr.push_back("CSimdArray<double> w_y(1, ket_pdim);");

        vstr.push_back("CSimdArray<double> w_z(1, ket_pdim);");
    }
    
    vstr.push_back("// allocate aligned distances R(PQ) = P - Q");

    vstr.push_back("CSimdArray<double> pq_x(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> pq_y(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> pq_z(1, ket_pdim);");
    
    if ((integral[2] + integral[3]) > 0)
    {
        vstr.push_back("// allocate aligned distances R(QD) = Q - D");
        
        vstr.push_back("CSimdArray<double> qd_x(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> qd_y(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> qd_z(1, ket_pdim);");
        
        vstr.push_back("// allocate aligned distances R(WQ) = W - Q");
        
        vstr.push_back("CSimdArray<double> wq_x(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> wq_y(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> wq_z(1, ket_pdim);");
    }
    
    if ((integral[0] + integral[1]) > 0)
    {
        vstr.push_back("// allocate aligned distances R(WP) = W - P");
        
        vstr.push_back("CSimdArray<double> wp_x(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> wp_y(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> wp_z(1, ket_pdim);");
    }
    
    vstr.push_back("// allocate combined overlap factor");
    
    vstr.push_back("CSimdArray<double> fss_abcd(1, ket_pdim);");
    
    if (integral[2] > 0)
    {
        vstr.push_back("// allocate and initialize aligned distances R(CD) = C - D");
        
        vstr.push_back("CSimdArray<double> cd_x(1, ket_dim);");
        
        vstr.push_back("CSimdArray<double> cd_y(1, ket_dim);");
        
        vstr.push_back("CSimdArray<double> cd_z(1, ket_dim);");
        
        vstr.push_back("t4cfunc::comp_distances_cd(cd_x[0], cd_y[0], cd_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], ket_dim);");
    }
    
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_diag_coordinates_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    const auto integrand = integral.integrand();
    
    vstr.push_back("// allocate aligned coordinates of Q center");

    vstr.push_back("CSimdArray<double> q_x(1, npgtos);");

    vstr.push_back("CSimdArray<double> q_y(1, npgtos);");

    vstr.push_back("CSimdArray<double> q_z(1, npgtos);");
        
    if ((integral[0] + integral[1] + integral[2] + integral[3]) > 0)
    {
        vstr.push_back("// allocate aligned coordinates of W center");

        vstr.push_back("CSimdArray<double> w_x(1, npgtos);");

        vstr.push_back("CSimdArray<double> w_y(1, npgtos);");

        vstr.push_back("CSimdArray<double> w_z(1, npgtos);");
    }
    
    vstr.push_back("// allocate aligned distances R(PQ) = P - Q");

    vstr.push_back("CSimdArray<double> pq_x(1, npgtos);");

    vstr.push_back("CSimdArray<double> pq_y(1, npgtos);");

    vstr.push_back("CSimdArray<double> pq_z(1, npgtos);");
    
    if ((integral[2] + integral[3]) > 0)
    {
        vstr.push_back("// allocate aligned distances R(QD) = Q - D");
        
        vstr.push_back("CSimdArray<double> qd_x(1, npgtos);");
        
        vstr.push_back("CSimdArray<double> qd_y(1, npgtos);");
        
        vstr.push_back("CSimdArray<double> qd_z(1, npgtos);");
        
        vstr.push_back("// allocate aligned distances R(WQ) = W - Q");
        
        vstr.push_back("CSimdArray<double> wq_x(1, npgtos);");
        
        vstr.push_back("CSimdArray<double> wq_y(1, npgtos);");
        
        vstr.push_back("CSimdArray<double> wq_z(1, npgtos);");
    }
    
    if ((integral[0] + integral[1]) > 0)
    {
        vstr.push_back("// allocate aligned distances R(WP) = W - P");
        
        vstr.push_back("CSimdArray<double> wp_x(1, npgtos);");
        
        vstr.push_back("CSimdArray<double> wp_y(1, npgtos);");
        
        vstr.push_back("CSimdArray<double> wp_z(1, npgtos);");
    }
    
    vstr.push_back("// allocate combined overlap factor");
    
    vstr.push_back("CSimdArray<double> fss_abcd(1, npgtos);");
    
    if (integral[2] > 0)
    {
        vstr.push_back("// allocate and initialize aligned distances R(CD) = C - D");
        
        vstr.push_back("CSimdArray<double> cd_x(1, 1);");
        
        vstr.push_back("CSimdArray<double> cd_y(1, 1);");
        
        vstr.push_back("CSimdArray<double> cd_z(1, 1);");
    }
    
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_full_coordinates_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    const auto integrand = integral.integrand();
    
    auto a_angmom = integral[0];
    
    auto b_angmom = integral[1];
    
    auto c_angmom = integral[2];
    
    auto d_angmom = integral[3];
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        a_angmom += prefixes[0].shape().order();
        
        b_angmom += prefixes[1].shape().order();
        
        c_angmom += prefixes[2].shape().order();
        
        d_angmom += prefixes[3].shape().order();
    }
    
    vstr.push_back("// allocate aligned coordinates of Q center");

    vstr.push_back("CSimdArray<double> q_x(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> q_y(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> q_z(1, ket_pdim);");
        
    if ((a_angmom + b_angmom + c_angmom + d_angmom) > 0)
    {
        vstr.push_back("// allocate aligned coordinates of W center");

        vstr.push_back("CSimdArray<double> w_x(1, ket_pdim);");

        vstr.push_back("CSimdArray<double> w_y(1, ket_pdim);");

        vstr.push_back("CSimdArray<double> w_z(1, ket_pdim);");
    }
    
    vstr.push_back("// allocate aligned distances R(PQ) = P - Q");

    vstr.push_back("CSimdArray<double> pq_x(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> pq_y(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> pq_z(1, ket_pdim);");
    
    if (c_angmom > 0)
    {
        vstr.push_back("// allocate aligned distances R(QC) = Q - C");
        
        vstr.push_back("CSimdArray<double> qc_x(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> qc_y(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> qc_z(1, ket_pdim);");
    }
    
    if (d_angmom > 0)
    {
        vstr.push_back("// allocate aligned distances R(QD) = Q - D");
        
        vstr.push_back("CSimdArray<double> qd_x(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> qd_y(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> qd_z(1, ket_pdim);");
    }
    
    if ((c_angmom + d_angmom) > 0)
    {
        vstr.push_back("// allocate aligned distances R(WQ) = W - Q");
        
        vstr.push_back("CSimdArray<double> wq_x(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> wq_y(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> wq_z(1, ket_pdim);");
    }
    
    if ((a_angmom + b_angmom) > 0)
    {
        vstr.push_back("// allocate aligned distances R(WP) = W - P");
        
        vstr.push_back("CSimdArray<double> wp_x(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> wp_y(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> wp_z(1, ket_pdim);");
    }
    
    vstr.push_back("// allocate combined overlap factor");
    
    vstr.push_back("CSimdArray<double> fss_abcd(1, ket_pdim);");
    
    return vstr;
}

SI4CIntegrals
T4CFuncBodyDriver::_get_cart_buffer_integrals(const SI4CIntegrals& bra_integrals,
                                              const SI4CIntegrals& ket_integrals) const
{
    SI4CIntegrals tints;
    
    for (const auto& tint : ket_integrals)
    {
        if ((tint[0] + tint[2]) == 0)
        {
            tints.insert(tint);
        }
    }
    
    for (const auto& tint : bra_integrals)
    {
        if ((tint[0] + tint[2]) == 0)
        {
            tints.insert(tint);
        }
    }
    
    return tints;
}

SI4CIntegrals
T4CFuncBodyDriver::_get_half_spher_buffers_integrals(const SI4CIntegrals& bra_integrals,
                                                     const SI4CIntegrals& ket_integrals,
                                                     const I4CIntegral&   integral) const
{
    SI4CIntegrals tints;
    
    std::vector<std::string> vstr;
    
    std::set<std::string> buffers;
    
    for (const auto& tint : ket_integrals)
    {
        if ((tint[2] == integral[2]) && (tint[3] == integral[3]))
        {
            tints.insert(tint);
        }
    }
    
    if (integral[0] > 0)
    {
        for (const auto& tint : bra_integrals)
        {
            if ((tint[0] >= 0) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
            {
                tints.insert(tint);
            }
        }
    }
    
    tints.insert(integral);
    
    return tints;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_prim_buffers_def(const SI4CIntegrals& integrals,
                                         const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned primitive integrals");
    
    for (const auto& tint : integrals)
    {
        if ((tint[0] + tint[2]) == 0)
        {
            std::string label = "CSimdArray<double> ";
            
            label += t4c::get_buffer_label(tint, "prim");
            
            const auto angpair = std::array<int, 2>({tint[1], tint[3]});
            
            const auto tcomps = t2c::number_of_cartesian_components(angpair);
            
            label += "(" + std::to_string(tcomps) + ", ket_pdim);";
            
            vstr.push_back(label);
        }
    }
    
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_diag_prim_buffers_def(const SI4CIntegrals& integrals,
                                              const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned primitive integrals");
    
    for (const auto& tint : integrals)
    {
        if ((tint[0] + tint[2]) == 0)
        {
            std::string label = "CSimdArray<double> ";
            
            label += t4c::get_buffer_label(tint, "prim");
            
            const auto angpair = std::array<int, 2>({tint[1], tint[3]});
            
            const auto tcomps = t2c::number_of_cartesian_components(angpair);
            
            label += "(" + std::to_string(tcomps) + ", npgtos);";
            
            vstr.push_back(label);
        }
    }
    
    return vstr;
}


std::vector<std::string>
T4CFuncBodyDriver::_get_full_prim_buffers_def(const SI4CIntegrals& integrals,
                                              const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned primitive integrals");
    
    for (const auto& tint : integrals)
    {
        std::string label = "CSimdArray<double> ";
            
        label += t4c::get_buffer_label(tint, "prim");
            
        auto angpair = std::array<int, 2>({tint[0], tint[1]});
            
        auto tcomps = t2c::number_of_cartesian_components(angpair);
        
        angpair = std::array<int, 2>({tint[2], tint[3]});
            
        tcomps *= t2c::number_of_cartesian_components(angpair);
            
        label += "(" + std::to_string(tcomps) + ", ket_pdim);";
            
        vstr.push_back(label);
    }
    
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_cart_buffers_def(const SI4CIntegrals& bra_integrals,
                                         const SI4CIntegrals& ket_integrals,
                                         const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned Cartesian integrals");
    
    for (const auto& tint : _get_cart_buffer_integrals(bra_integrals, ket_integrals))
    {
        if ((tint[0] + tint[2]) == 0)
        {
            std::string label = "CSimdArray<double> ";
            
            label += t4c::get_buffer_label(tint, "cart");
            
            const auto angpair = std::array<int, 2>({tint[1], tint[3]});
            
            const auto tcomps = t2c::number_of_cartesian_components(angpair);
            
            label += "(" + std::to_string(tcomps) + ", ket_dim);";
            
            vstr.push_back(label);
        }
    }
    
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_diag_cart_buffers_def(const SI4CIntegrals& bra_integrals,
                                              const SI4CIntegrals& ket_integrals,
                                              const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned Cartesian integrals");
    
    for (const auto& tint : _get_cart_buffer_integrals(bra_integrals, ket_integrals))
    {
        if ((tint[0] + tint[2]) == 0)
        {
            std::string label = "CSimdArray<double> ";
            
            label += t4c::get_buffer_label(tint, "cart");
            
            const auto angpair = std::array<int, 2>({tint[1], tint[3]});
            
            const auto tcomps = t2c::number_of_cartesian_components(angpair);
            
            label += "(" + std::to_string(tcomps) + ", 1);";
            
            vstr.push_back(label);
        }
    }
    
    return vstr;
}


std::vector<std::string>
T4CFuncBodyDriver::_get_full_cart_buffers_def(const SI4CIntegrals& integrals,
                                              const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned Cartesian integrals");
    
    std::string label = "CSimdArray<double> ";
            
    label += t4c::get_buffer_label(integral, "cart");
            
    const auto tcomps = integral.components<T2CPair, T2CPair>().size();
        
    label += "(" + std::to_string(tcomps) + ", ket_dim);";
            
    vstr.push_back(label);
    
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_contr_buffers_def(const SI4CIntegrals& bra_integrals,
                                          const SI4CIntegrals& ket_integrals,
                                          const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    for (const auto& tint : ket_integrals)
    {
        if ((tint[0] == 0) && (tint[2] > 0))
        {
            std::string label = "CSimdArray<double> ";
            
            label += t4c::get_buffer_label(tint, "contr");
            
            const auto angpair = std::array<int, 2>({tint[2], tint[3]});
            
            auto tcomps = t2c::number_of_cartesian_components(angpair);
            
            tcomps *= t2c::number_of_cartesian_components(tint[1]);
            
            label += "(" + std::to_string(tcomps) + ", ket_dim);";
            
            vstr.push_back(label);
        }
    }
    
    if (!vstr.empty())
    {
        vstr.insert(vstr.begin(), "// allocate aligned contracted integrals");
    }
    
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_diag_contr_buffers_def(const SI4CIntegrals& bra_integrals,
                                               const SI4CIntegrals& ket_integrals,
                                               const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    for (const auto& tint : ket_integrals)
    {
        if ((tint[0] == 0) && (tint[2] > 0))
        {
            std::string label = "CSimdArray<double> ";
            
            label += t4c::get_buffer_label(tint, "contr");
            
            const auto angpair = std::array<int, 2>({tint[2], tint[3]});
            
            auto tcomps = t2c::number_of_cartesian_components(angpair);
            
            tcomps *= t2c::number_of_cartesian_components(tint[1]);
            
            label += "(" + std::to_string(tcomps) + ", 1);";
            
            vstr.push_back(label);
        }
    }
    
    if (!vstr.empty())
    {
        vstr.insert(vstr.begin(), "// allocate aligned contracted integrals");
    }
    
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_half_spher_buffers_def(const SI4CIntegrals& bra_integrals,
                                               const SI4CIntegrals& ket_integrals,
                                               const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned half transformed integrals");
    
    for (const auto& tint : _get_half_spher_buffers_integrals(bra_integrals, ket_integrals, integral))
    {
        std::string label = "CSimdArray<double> ";
                
        label += t4c::get_buffer_label(tint, "ket_spher");
                
        auto angpair = std::array<int, 2>({tint[2], tint[3]});
                
        auto tcomps = t2c::number_of_spherical_components(angpair);
            
        angpair = std::array<int, 2>({tint[0], tint[1]});
                
        tcomps *= t2c::number_of_cartesian_components(angpair);
                
        label += "(" + std::to_string(tcomps) + ", ket_dim);";
                
        vstr.push_back(label);
    }
        
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_diag_half_spher_buffers_def(const SI4CIntegrals& bra_integrals,
                                                    const SI4CIntegrals& ket_integrals,
                                                    const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned half transformed integrals");
    
    for (const auto& tint : _get_half_spher_buffers_integrals(bra_integrals, ket_integrals, integral))
    {
        std::string label = "CSimdArray<double> ";
                
        label += t4c::get_buffer_label(tint, "ket_spher");
                
        auto angpair = std::array<int, 2>({tint[2], tint[3]});
                
        auto tcomps = t2c::number_of_spherical_components(angpair);
            
        angpair = std::array<int, 2>({tint[0], tint[1]});
                
        tcomps *= t2c::number_of_cartesian_components(angpair);
                
        label += "(" + std::to_string(tcomps) + ", 1);";
                
        vstr.push_back(label);
    }
        
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_spher_buffers_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned spherical integrals");
     
    std::string label = "CSimdArray<double> ";
                    
    label += t4c::get_buffer_label(integral, "spher");
                    
    auto angpair = std::array<int, 2>({integral[2], integral[3]});
                    
    auto tcomps = t2c::number_of_spherical_components(angpair);
                    
    angpair = std::array<int, 2>({integral[0], integral[1]});
                    
    tcomps *= t2c::number_of_spherical_components(angpair);
                    
    label += "(" + std::to_string(tcomps) + ", ket_dim);";
                    
    vstr.push_back(label);
   
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_diag_spher_buffers_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned spherical integrals");
     
    std::string label = "CSimdArray<double> ";
                    
    label += t4c::get_buffer_label(integral, "spher");
                    
    auto angpair = std::array<int, 2>({integral[2], integral[3]});
                    
    auto tcomps = t2c::number_of_spherical_components(angpair);
                    
    angpair = std::array<int, 2>({integral[0], integral[1]});
                    
    tcomps *= t2c::number_of_spherical_components(angpair);
                    
    label += "(" + std::to_string(tcomps) + ", 1);";
                    
    vstr.push_back(label);
   
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_full_spher_buffers_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned spherical integrals");
     
    std::string label = "CSimdArray<double> ";
                    
    label += t4c::get_buffer_label(integral, "spher");
                    
    auto angpair = std::array<int, 2>({integral[2], integral[3]});
                    
    auto tcomps = t2c::number_of_spherical_components(angpair);
                    
    angpair = std::array<int, 2>({integral[0], integral[1]});
                    
    tcomps *= t2c::number_of_spherical_components(angpair);
    
    for (const auto& prefix : integral.prefixes())
    {
        tcomps *= prefix.components().size();
    }
                    
    label += "(" + std::to_string(tcomps) + ", ket_dim);";
                    
    vstr.push_back(label);
   
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_boys_function_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto order = integral[0] + integral[1] + integral[2] + integral[3];
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        for (const auto& prefix : prefixes)
        {
            order += prefix.shape().order(); 
        }
    }
        
    vstr.push_back("// setup Boys fuction data");
        
    vstr.push_back("const CBoysFunc<" + std::to_string(order) + "> bf_table;");
        
    vstr.push_back("CSimdArray<double> bf_args(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> bf_values(" + std::to_string(order + 1) + ", ket_pdim);");
   
    return vstr;
}

std::vector<std::string>
T4CFuncBodyDriver::_get_diag_boys_function_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto order = integral[0] + integral[1] + integral[2] + integral[3];
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        for (const auto& prefix : prefixes)
        {
            order += prefix.shape().order();
        }
    }
        
    vstr.push_back("// setup Boys fuction data");
        
    vstr.push_back("const CBoysFunc<" + std::to_string(order) + "> bf_table;");
        
    vstr.push_back("CSimdArray<double> bf_args(1, npgtos);");

    vstr.push_back("CSimdArray<double> bf_values(" + std::to_string(order + 1) + ", npgtos);");
   
    return vstr;
}

void
T4CFuncBodyDriver::_add_loop_start(      VCodeLines&    lines,
                                   const SI4CIntegrals& bra_integrals,
                                   const SI4CIntegrals& ket_integrals,
                                   const I4CIntegral&   integral,
                                   const bool           diagonal) const
{
    lines.push_back({1, 0, 2, "// loop over contracted GTOs on bra side"});
        
    lines.push_back({1, 0, 1, "for (auto i = bra_indices[0]; i < bra_indices[1]; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    for (const auto& tint : _get_cart_buffer_integrals(bra_integrals, ket_integrals))
    {
        std::string label = t4c::get_buffer_label(tint, "cart")  + ".zero();";
            
        lines.push_back({2, 0, 2, label});
    }
    
    for (const auto& tint : _get_half_spher_buffers_integrals(bra_integrals, ket_integrals, integral))
    {
        std::string label = t4c::get_buffer_label(tint, "ket_spher")  + ".zero();";
            
        lines.push_back({2, 0, 2, label});
    }
    
    std::string label = t4c::get_buffer_label(integral, "spher")  + ".zero();";
        
    lines.push_back({2, 0, 2, label});

    lines.push_back({2, 0, 2, "const auto a_x = a_coords_x[i];"});

    lines.push_back({2, 0, 2, "const auto a_y = a_coords_y[i];"});

    lines.push_back({2, 0, 2, "const auto a_z = a_coords_z[i];"});
    
    lines.push_back({2, 0, 2, "const auto b_x = b_coords_x[i];"});

    lines.push_back({2, 0, 2, "const auto b_y = b_coords_y[i];"});

    lines.push_back({2, 0, 2, "const auto b_z = b_coords_z[i];"});
    
    if (integral[0] > 0)
    {
        lines.push_back({2, 0, 2, "const auto ab_x = a_x - b_x;"});

        lines.push_back({2, 0, 2, "const auto ab_y = a_y - b_y;"});

        lines.push_back({2, 0, 2, "const auto ab_z = a_z - b_z;"});
    }
}

void
T4CFuncBodyDriver::_add_full_loop_start(      VCodeLines&    lines,
                                        const SI4CIntegrals& integrals,
                                        const I4CIntegral&   integral) const
{
    lines.push_back({1, 0, 2, "// loop over contracted GTOs on bra side"});
        
    lines.push_back({1, 0, 1, "for (auto i = bra_indices[0]; i < bra_indices[1]; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    std::string label = t4c::get_buffer_label(integral, "cart")  + ".zero();";
        
    lines.push_back({2, 0, 2, label});
    
    label = t4c::get_buffer_label(integral, "spher")  + ".zero();";
        
    lines.push_back({2, 0, 2, label});

    lines.push_back({2, 0, 2, "const auto a_x = a_coords_x[i];"});

    lines.push_back({2, 0, 2, "const auto a_y = a_coords_y[i];"});

    lines.push_back({2, 0, 2, "const auto a_z = a_coords_z[i];"});
    
    lines.push_back({2, 0, 2, "const auto b_x = b_coords_x[i];"});

    lines.push_back({2, 0, 2, "const auto b_y = b_coords_y[i];"});

    lines.push_back({2, 0, 2, "const auto b_z = b_coords_z[i];"});
}

void
T4CFuncBodyDriver::_add_loop_end(      VCodeLines&  lines,
                                 const I4CIntegral& integral,
                                 const bool         diagonal) const
{
    std::string label = "distributor->distribute(";

    label +=  t4c::get_buffer_label(integral, "spher") + ", ";
        
    if (diagonal)
    {
        label += "a_indices, b_indices, ";
    }
    else
    {
        label += "a_indices, b_indices, c_indices, d_indices, ";
    }
            
    label += std::to_string(integral[0]) + ", ";
            
    label += std::to_string(integral[1]) + ", ";
    
    if (!diagonal)
    {
        label += std::to_string(integral[2]) + ", ";
                
        label += std::to_string(integral[3]) + ", ";
    }
            
    label += "i, ket_indices);";
            
    lines.push_back({2, 0, 1, label});
   
    lines.push_back({1, 0, 1, "}"});
}

void
T4CFuncBodyDriver::_add_full_loop_end(      VCodeLines&  lines,
                                      const I4CIntegral& integral) const
{
    std::string label = "distributor->distribute(";

    label +=  t4c::get_buffer_label(integral, "spher") + ", ";
        
    label += "a_indices, b_indices, c_indices, d_indices, ";
            
    label += std::to_string(integral[0]) + ", ";
            
    label += std::to_string(integral[1]) + ", ";
    
    label += std::to_string(integral[2]) + ", ";
                
    label += std::to_string(integral[3]) + ", ";
    
    const auto prefixes = integral.prefixes();
    
    if (!prefixes.empty())
    {
        label += std::to_string(prefixes[0].shape().order()) + ", ";
                
        label += std::to_string(prefixes[1].shape().order()) + ", ";
        
        label += std::to_string(prefixes[2].shape().order()) + ", ";
                    
        label += std::to_string(prefixes[3].shape().order()) + ", ";
    }
            
    label += "i, ket_indices);";
            
    lines.push_back({2, 0, 1, label});
   
    lines.push_back({1, 0, 1, "}"});
}

void
T4CFuncBodyDriver::_add_ket_loop_start(      VCodeLines&  lines,
                                       const I4CIntegral& integral,
                                       const bool         diagonal) const
{
    if (diagonal)
    {
        lines.push_back({2, 0, 1, "for (int j = 0; j < npgtos; j++)"});
    }
    else
    {
        lines.push_back({2, 0, 1, "for (int j = 0; j < bra_npgtos; j++)"});
    }
    
    lines.push_back({2, 0, 1, "{"});
    
    if (diagonal && (integral[0] == integral[2]) && (integral[1] == integral[3]))
    {
        lines.push_back({3, 0, 2, "const auto a_exp = a_vec_exps[j * ncgtos + i];"});
        
        lines.push_back({3, 0, 2, "const auto b_exp = b_vec_exps[j * ncgtos + i];"});
            
        lines.push_back({3, 0, 2, "const auto ab_norm = ab_vec_norms[j * ncgtos + i];"});
        
        lines.push_back({3, 0, 2, "const auto ab_ovl = ab_vec_ovls[j * ncgtos + i];"});
    }
    else
    {
        lines.push_back({3, 0, 2, "const auto a_exp = a_vec_exps[j * bra_ncgtos + i];"});
        
        lines.push_back({3, 0, 2, "const auto b_exp = b_vec_exps[j * bra_ncgtos + i];"});
            
        lines.push_back({3, 0, 2, "const auto ab_norm = ab_vec_norms[j * bra_ncgtos + i];"});
        
        lines.push_back({3, 0, 2, "const auto ab_ovl = ab_vec_ovls[j * bra_ncgtos + i];"});
    }
    
    lines.push_back({3, 0, 2, "const auto p_x = (a_x * a_exp + b_x * b_exp) / (a_exp + b_exp);"});
    
    lines.push_back({3, 0, 2, "const auto p_y = (a_y * a_exp + b_y * b_exp) / (a_exp + b_exp);"});
    
    lines.push_back({3, 0, 2, "const auto p_z = (a_z * a_exp + b_z * b_exp) / (a_exp + b_exp);"});
    
    if ((integral[0] + integral[1]) > 0)
    {
        lines.push_back({3, 0, 2, "const auto pb_x = p_x - b_x;"});
        
        lines.push_back({3, 0, 2, "const auto pb_y = p_y - b_y;"});
        
        lines.push_back({3, 0, 2, "const auto pb_z = p_z - b_z;"});
    }
    
    lines.push_back({3, 0, 2, "t4cfunc::comp_coordinates_q(q_x[0], q_y[0], q_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], c_exps[0], d_exps[0], ket_pdim);"});
    
    if ((integral[0] + integral[1] + integral[2] + integral[3]) > 0)
    {
        lines.push_back({3, 0, 2, "t4cfunc::comp_coordinates_w(w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], a_exp, b_exp, c_exps[0], d_exps[0], ket_pdim);"});
    }
    
    lines.push_back({3, 0, 2, "t4cfunc::comp_distances_pq(pq_x[0], pq_y[0], pq_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], ket_pdim);"});
    
    if ((integral[2] + integral[3]) > 0)
    {
        lines.push_back({3, 0, 2, "t4cfunc::comp_distances_wq(wq_x[0], wq_y[0], wq_z[0], w_x[0], w_y[0], w_z[0], q_x[0], q_y[0], q_z[0], ket_pdim);"});
        
        lines.push_back({3, 0, 2, "t4cfunc::comp_distances_qd(qd_x[0], qd_y[0], qd_z[0], q_x[0], q_y[0], q_z[0], d_x[0], d_y[0], d_z[0], ket_pdim);"});
    }
    
    if ((integral[0] + integral[1]) > 0)
    {
        lines.push_back({3, 0, 2, "t4cfunc::comp_distances_wp(wp_x[0], wp_y[0], wp_z[0], w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, ket_pdim);"});
    }
    
    lines.push_back({3, 0, 2, "t4cfunc::comp_boys_args(bf_args, pq_x[0], pq_y[0], pq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);"});
    
    lines.push_back({3, 0, 2, "bf_table.compute(bf_values, bf_args);"});
    
    lines.push_back({3, 0, 2, "t4cfunc::comp_ovl_factors(fss_abcd, ab_ovl, cd_ovls[0], ab_norm, cd_norms[0], a_exp, b_exp, c_exps[0], d_exps[0]);"});
}

void
T4CFuncBodyDriver::_add_full_ket_loop_start(      VCodeLines&  lines,
                                            const I4CIntegral& integral) const
{
    lines.push_back({2, 0, 1, "for (int j = 0; j < bra_npgtos; j++)"});
    
    lines.push_back({2, 0, 1, "{"});
    
    lines.push_back({3, 0, 2, "const auto a_exp = a_vec_exps[j * bra_ncgtos + i];"});
    
    lines.push_back({3, 0, 2, "const auto b_exp = b_vec_exps[j * bra_ncgtos + i];"});
        
    lines.push_back({3, 0, 2, "const auto ab_norm = ab_vec_norms[j * bra_ncgtos + i];"});
    
    lines.push_back({3, 0, 2, "const auto ab_ovl = ab_vec_ovls[j * bra_ncgtos + i];"});
    
    lines.push_back({3, 0, 2, "const auto p_x = (a_x * a_exp + b_x * b_exp) / (a_exp + b_exp);"});
    
    lines.push_back({3, 0, 2, "const auto p_y = (a_y * a_exp + b_y * b_exp) / (a_exp + b_exp);"});
    
    lines.push_back({3, 0, 2, "const auto p_z = (a_z * a_exp + b_z * b_exp) / (a_exp + b_exp);"});
    
    auto a_angmom = integral[0];
    
    auto b_angmom = integral[1];
    
    auto c_angmom = integral[2];
    
    auto d_angmom = integral[3];
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        a_angmom += prefixes[0].shape().order();
        
        b_angmom += prefixes[1].shape().order();
        
        c_angmom += prefixes[2].shape().order();
        
        d_angmom += prefixes[3].shape().order();
    }
    
    if (a_angmom > 0)
    {
        lines.push_back({3, 0, 2, "const auto pa_x = p_x - a_x;"});
        
        lines.push_back({3, 0, 2, "const auto pa_y = p_y - a_y;"});
        
        lines.push_back({3, 0, 2, "const auto pa_z = p_z - a_z;"});
    }
    
    if (b_angmom > 0)
    {
        lines.push_back({3, 0, 2, "const auto pb_x = p_x - b_x;"});
        
        lines.push_back({3, 0, 2, "const auto pb_y = p_y - b_y;"});
        
        lines.push_back({3, 0, 2, "const auto pb_z = p_z - b_z;"});
    }
    
    lines.push_back({3, 0, 2, "t4cfunc::comp_coordinates_q(q_x[0], q_y[0], q_z[0], c_x[0], c_y[0], c_z[0], d_x[0], d_y[0], d_z[0], c_exps[0], d_exps[0], ket_pdim);"});
    
    if ((a_angmom + b_angmom + c_angmom + d_angmom) > 0)
    {
        lines.push_back({3, 0, 2, "t4cfunc::comp_coordinates_w(w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], a_exp, b_exp, c_exps[0], d_exps[0], ket_pdim);"});
    }
    
    lines.push_back({3, 0, 2, "t4cfunc::comp_distances_pq(pq_x[0], pq_y[0], pq_z[0], p_x, p_y, p_z, q_x[0], q_y[0], q_z[0], ket_pdim);"});
    
    if ((c_angmom + d_angmom) > 0)
    {
        lines.push_back({3, 0, 2, "t4cfunc::comp_distances_wq(wq_x[0], wq_y[0], wq_z[0], w_x[0], w_y[0], w_z[0], q_x[0], q_y[0], q_z[0], ket_pdim);"});
    }
    
    if (c_angmom > 0)
    {
        lines.push_back({3, 0, 2, "t4cfunc::comp_distances_qc(qc_x[0], qc_y[0], qc_z[0], q_x[0], q_y[0], q_z[0], c_x[0], c_y[0], c_z[0], ket_pdim);"});
    }
    
    if (d_angmom > 0)
    {
        lines.push_back({3, 0, 2, "t4cfunc::comp_distances_qd(qd_x[0], qd_y[0], qd_z[0], q_x[0], q_y[0], q_z[0], d_x[0], d_y[0], d_z[0], ket_pdim);"});
    }
    
    if ((a_angmom + b_angmom) > 0)
    {
        lines.push_back({3, 0, 2, "t4cfunc::comp_distances_wp(wp_x[0], wp_y[0], wp_z[0], w_x[0], w_y[0], w_z[0], p_x, p_y, p_z, ket_pdim);"});
    }
    
    lines.push_back({3, 0, 2, "t4cfunc::comp_boys_args(bf_args, pq_x[0], pq_y[0], pq_z[0], a_exp, b_exp, c_exps[0], d_exps[0]);"});
    
    lines.push_back({3, 0, 2, "bf_table.compute(bf_values, bf_args);"});
    
    lines.push_back({3, 0, 2, "t4cfunc::comp_ovl_factors(fss_abcd, ab_ovl, cd_ovls[0], ab_norm, cd_norms[0], a_exp, b_exp, c_exps[0], d_exps[0]);"});
}

void
T4CFuncBodyDriver::_add_ket_loop_end(      VCodeLines&  lines,
                                     const SI4CIntegrals& bra_integrals,
                                     const SI4CIntegrals& ket_integrals,
                                     const I4CIntegral& integral,
                                     const bool         diagonal) const
{
    for (const auto& tint : ket_integrals)
    {
        if ((tint[0] + tint[2]) == 0)
        {
            std::string label = "t2cfunc::reduce(";
            
            label += t4c::get_buffer_label(tint, "cart") + ", ";
            
            label += t4c::get_buffer_label(tint, "prim") + ", ";
            
            if (diagonal)
            {
                label += "ket_dim, npgtos);";
            }
            else
            {
                label += "ket_dim, ket_npgtos);";
            }
            
            lines.push_back({3, 0, 2, label});
        }
    }
    
    for (const auto& tint : bra_integrals)
    {
        if ((tint[0] + tint[2]) == 0)
        {
            std::string label = "t2cfunc::reduce(";
            
            label += t4c::get_buffer_label(tint, "cart") + ", ";
            
            label += t4c::get_buffer_label(tint, "prim") + ", ";
            
            if (diagonal)
            {
                label += "ket_dim, npgtos);";
            }
            else
            {
                label += "ket_dim, ket_npgtos);";
            }
            
            lines.push_back({3, 0, 2, label});
        }
    }
    
    lines.push_back({2, 0, 2, "}"});
}

void
T4CFuncBodyDriver::_add_full_ket_loop_end(      VCodeLines&  lines,
                                          const I4CIntegral& integral) const
{
    std::string label = "t2cfunc::reduce(";
    
    label += t4c::get_buffer_label(integral, "cart") + ", ";
    
    label += t4c::get_buffer_label(integral, "prim") + ", ";
    
    label += "ket_dim, ket_npgtos);";
    
    lines.push_back({3, 0, 1, label});
    
    lines.push_back({2, 0, 2, "}"});
}

void
T4CFuncBodyDriver::_add_auxilary_integrals(      VCodeLines&    lines,
                                           const SI4CIntegrals& integrals) const
{
    for (const auto& tint : integrals)
    {
        if ((tint[0] + tint[1] + tint[2] + tint[3]) == 0)
        {
            const auto label = std::to_string(tint.order());
                    
            lines.push_back({3, 0, 2, "erirec::comp_prim_electron_repulsion_ssss(prim_buffer_" + label + "_ssss, fss_abcd[0], bf_values[" + label + "]);"});
            
            // TODO: other integrals...
        }
    }
}

void
T4CFuncBodyDriver::_add_vrr_call_tree(      VCodeLines&  lines,
                                      const SI4CIntegrals& integrals) const
{
    for (const auto& tint : integrals)
    {
        if (((tint[0] + tint[2]) == 0) && ((tint[1] + tint[3]) > 0))
        {
            const auto name = t4c::prim_compute_func_name(tint);
            
            auto label = t4c::namespace_label(tint) + "::" + name + "(";
            
            label += _get_vrr_arguments(tint);
            
            if (tint[1] > 0)
            {
                if ((tint[1] == 1) && (tint[3] == 0))
                {
                    label += "pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]";
                }
                else
                {
                    label += "pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], ";
                }
            }
            
            if ((tint[3] > 0) && (tint[1] == 0))
            {
                if ((tint[1] == 0) && (tint[3] == 1))
                {
                    label += "qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0]";
                }
                else
                {
                    label += "qd_x[0], qd_y[0], qd_z[0], wq_x[0], wq_y[0], wq_z[0], ";
                }
            }
            
            if ((tint[1] + tint[3]) > 1)
            {
                label += "a_exp, b_exp, c_exps[0], d_exps[0]";
            }
            
            label += ");";
            
            lines.push_back({3, 0, 2, label});
        }
    }
}

void
T4CFuncBodyDriver::_add_full_vrr_call_tree(      VCodeLines&  lines,
                                           const SI4CIntegrals& integrals) const
{
    for (const auto& tint : integrals)
    {
        if ((tint[0] + tint[1] + tint[2] + tint[3]) > 0)
        {
            const auto name = t4c::prim_compute_func_name(tint);
            
            auto label = t4c::namespace_label(tint) + "::" + name + "(";
            
            label += _get_full_vrr_arguments(tint);
            
            if (tint[0] > 0)
            {
                if ((tint[0] == 1) && ((tint[1] + tint[2] + tint[3]) == 0))
                {
                    label += "pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0]";
                }
                else
                {
                    label += "pa_x, pa_y, pa_z, wp_x[0], wp_y[0], wp_z[0], ";
                }
            }
            
            if ((tint[1] > 0) && (tint[0] == 0))
            {
                if ((tint[1] == 1) && ((tint[2] + tint[3]) == 0))
                {
                    label += "pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0]";
                }
                else
                {
                    label += "pb_x, pb_y, pb_z, wp_x[0], wp_y[0], wp_z[0], ";
                }
            }
            
            if ((tint[2] > 0) && ((tint[0] + tint[1]) == 0))
            {
                if ((tint[2] == 1) && (tint[3] == 0))
                {
                    label += "qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0]";
                }
                else
                {
                    label += "qc_x, qc_y, qc_z, wq_x[0], wq_y[0], wq_z[0], ";
                }
            }
            
            if ((tint[3] > 0) && ((tint[0] + tint[1] + tint[2]) == 0))
            {
                if (tint[3] == 1)
                {
                    label += "qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0]";
                }
                else
                {
                    label += "qd_x, qd_y, qd_z, wq_x[0], wq_y[0], wq_z[0], ";
                }
            }
            
            if ((tint[0] + tint[1] + tint[2] + tint[3]) > 1)
            {
                label += "a_exp, b_exp, c_exps[0], d_exps[0]";
            }
            
            label += ");";
            
            lines.push_back({3, 0, 2, label});
        }
    }
}


std::string
T4CFuncBodyDriver::_get_vrr_arguments(const I4CIntegral& integral) const
{
    std::string label = t4c::get_buffer_label(integral, {"prim"}) + ", ";;
    
    for (const auto& tint : t4c::get_vrr_integrals(integral))
    {
        label += t4c::get_buffer_label(tint, {"prim"}) + ", ";
    }
    
    return label;
    
}

std::string
T4CFuncBodyDriver::_get_full_vrr_arguments(const I4CIntegral& integral) const
{
    std::string label = t4c::get_buffer_label(integral, {"prim"}) + ", ";;
    
    for (const auto& tint : t4c::get_full_vrr_integrals(integral))
    {
        label += t4c::get_buffer_label(tint, {"prim"}) + ", ";
    }
    
    return label;
}

void
T4CFuncBodyDriver::_add_geom_call_tree(      VCodeLines&    lines,
                                       const SI4CIntegrals& integrals,
                                       const I4CIntegral&   integral) const
{
    const auto name = t4c::geom_compute_func_name(integral);
    
    auto label = t4c::geom_namespace_label() + "::" + name + "(";
    
    label += _get_geom_arguments(integrals, integral);
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        if (prefixes[0].shape().order() > 0) label += "a_exp, ";
        
        if (prefixes[1].shape().order() > 0) label += "b_exp, ";
        
        if (prefixes[2].shape().order() > 0) label += "c_exps[0], ";
        
        if (prefixes[3].shape().order() > 0) label += "d_exps[0]";
    }
    
    if (label[label.size() - 2] == ',') label.erase(label.end() - 2);
    
    label += ");";
    
    lines.push_back({3, 0, 2, label});
}

std::string
T4CFuncBodyDriver::_get_geom_arguments(const SI4CIntegrals& integrals,
                                       const I4CIntegral&   integral) const
{
    std::string label = t4c::get_buffer_label(integral, {"prim"}) + ", ";
    
    SI4CIntegrals ref_tints;
    
    for (const auto& tint : integrals)
    {
        ref_tints.insert(tint.base());
    }
    
    for (const auto& tint : ref_tints)
    {
        label += t4c::get_buffer_label(tint, {"prim"}) + ", ";
    }
    
    return label;
}

void
T4CFuncBodyDriver::_add_ket_hrr_call_tree(      VCodeLines&  lines,
                                          const SI4CIntegrals& integrals) const
{
    for (const auto& tint : integrals)
    {
        if ((tint[0] == 0) && (tint[2] > 0))
        {
            const auto name = t4c::ket_hrr_compute_func_name(tint);
            
            auto label = t4c::namespace_label(tint) + "::" + name + "(";
            
            label += _get_ket_hrr_arguments(tint);
            
            label += "cd_x[0], cd_y[0], cd_z[0], ";
            
            label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]);
            
            label += ");";
            
            lines.push_back({2, 0, 2, label});
        }
    }
}

std::string
T4CFuncBodyDriver::_get_ket_hrr_arguments(const I4CIntegral& integral) const
{
    std::string label = t4c::get_buffer_label(integral, {"contr"}) + ", ";;
    
    for (const auto& tint : t4c::get_ket_hrr_integrals(integral))
    {
        if (tint[2] == 0)
        {
            label += t4c::get_buffer_label(tint, {"cart"}) + ", ";
        }
        else
        {
            label += t4c::get_buffer_label(tint, {"contr"}) + ", ";
        }
       
    }
    
    return label;
}

void
T4CFuncBodyDriver::_add_ket_trafo_call_tree(      VCodeLines&  lines,
                                            const SI4CIntegrals& bra_integrals,
                                            const SI4CIntegrals& ket_integrals,
                                            const I4CIntegral&   integral) const
{
    for (const auto& tint : ket_integrals)
    {
        if ((tint[0] == 0) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
        {
            std::string label = "t4cfunc::ket_transform<" + std::to_string(tint[2]) + ", " + std::to_string(tint[3]) + ">";
                
            label += "("  +  t4c::get_buffer_label(tint, "ket_spher") + ", ";
            
            label += t4c::get_buffer_label(tint, "contr") + ", ";
            
            label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]) + ");";
                
            lines.push_back({2, 0, 2, label});
        }
    }
    
    if ((integral[0] > 0) && (integral[2] == 0))
    {
        for (const auto& tint : bra_integrals)
        {
            if ((tint[0] == 0) && (tint[2] == 0))
            {
                std::string label = "t4cfunc::ket_transform<" + std::to_string(tint[2]) + ", " + std::to_string(tint[3]) + ">";
                    
                label += "("  +  t4c::get_buffer_label(tint, "ket_spher") + ", ";
                
                label += t4c::get_buffer_label(tint, "cart") + ", ";
                
                label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]) + ");";
                    
                lines.push_back({2, 0, 2, label});
            }
        }
    }
    
    if ((integral[0] == 0) && (integral[2] == 0))
    {
        std::string label = "t4cfunc::ket_transform<" + std::to_string(integral[2]) + ", " + std::to_string(integral[3]) + ">";
            
        label += "("  +  t4c::get_buffer_label(integral, "ket_spher") + ", ";
        
        label += t4c::get_buffer_label(integral, "cart") + ", ";
        
        label += std::to_string(integral[0]) + ", " + std::to_string(integral[1]) + ");";
            
        lines.push_back({2, 0, 2, label});
    }
}

void
T4CFuncBodyDriver::_add_bra_hrr_call_tree(      VCodeLines&  lines,
                                          const SI4CIntegrals& integrals) const
{
    for (const auto& tint : integrals)
    {
        if (tint[0] > 0)
        {
            const auto name = t4c::bra_hrr_compute_func_name(tint);
            
            auto label = t4c::namespace_label(tint) + "::" + name + "(";
            
            label += _get_bra_hrr_arguments(tint);
            
            label += "ab_x, ab_y, ab_z, ";
            
            label += std::to_string(tint[2]) + ", " + std::to_string(tint[3]);
            
            label += ");";
            
            lines.push_back({2, 0, 2, label});
        }
    }
}

std::string
T4CFuncBodyDriver::_get_bra_hrr_arguments(const I4CIntegral& integral) const
{
    std::string label = t4c::get_buffer_label(integral, {"ket_spher"}) + ", ";;
    
    for (const auto& tint : t4c::get_bra_hrr_integrals(integral))
    {
        label += t4c::get_buffer_label(tint, {"ket_spher"}) + ", ";
    }
    
    return label;
}

void
T4CFuncBodyDriver::_add_bra_trafo_call_tree(      VCodeLines&  lines,
                                            const I4CIntegral& integral) const
{
    std::string label = "t4cfunc::bra_transform<" + std::to_string(integral[0]) + ", " + std::to_string(integral[1]) + ">";
        
    label += "("  +  t4c::get_buffer_label(integral, "spher") + ", ";
    
    label += t4c::get_buffer_label(integral, "ket_spher") + ", ";
    
    label += std::to_string(integral[2]) + ", " + std::to_string(integral[3]) + ");";
        
    lines.push_back({2, 0, 2, label});
}

void
T4CFuncBodyDriver::_add_full_trafo(      VCodeLines&  lines,
                                   const I4CIntegral& integral) const
{
    std::string label = "t4cfunc::full_transform<" + std::to_string(integral[0]) + ", ";
    
    label += std::to_string(integral[1]) + ", ";
    
    label += std::to_string(integral[2]) + ", ";
    
    label += std::to_string(integral[3]) + ">";
        
    label += "("  +  t4c::get_buffer_label(integral, "spher") + ", ";
    
    label += t4c::get_buffer_label(integral, "cart_spher") + ");";
        
    lines.push_back({2, 0, 2, label});
}
