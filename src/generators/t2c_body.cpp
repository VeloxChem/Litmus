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

#include <iostream>
#include "t2c_body.hpp"

#include "t2c_utils.hpp"

void
T2CFuncBodyDriver::write_func_body(      std::ofstream&         fstream,
                                   const SI2CIntegrals&         geom_integrals,
                                   const SI2CIntegrals&         vrr_integrals,
                                   const I2CIntegral&           integral,
                                   const std::array<int, 3>& geom_drvs, 
                                   const std::pair<bool, bool>& rec_form,
                                   const bool                   use_rs) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_external_data_def(integral, rec_form))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_gtos_def())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_ket_variables_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_buffers_def(vrr_integrals, integral, geom_drvs))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_boys_function_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_loop_start(lines, integral);
    
    _add_ket_loop_start(lines, integral, rec_form);
    
    _add_auxilary_integrals(lines, vrr_integrals, integral, rec_form, false);
    
    _add_sum_loop_start(lines, integral, rec_form, use_rs); 
    
    _add_auxilary_integrals(lines, vrr_integrals, integral, rec_form, true);
    
    _add_call_tree(lines, vrr_integrals, integral, rec_form);
    
    _add_geom_call_tree(lines, geom_integrals, vrr_integrals, integral, geom_drvs, rec_form);
    
    _add_sum_loop_end(lines, vrr_integrals, integral, rec_form);
    
    _add_ket_loop_end(lines, vrr_integrals, integral, rec_form);
    
    _add_loop_end(lines, integral, rec_form);
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CFuncBodyDriver::_get_external_data_def(const I2CIntegral&           integral,
                                          const std::pair<bool, bool>& rec_form) const
{
    std::vector<std::string> vstr;
    
    if (_need_external_coords(integral))
    {
        vstr.push_back("// intialize external coordinate(s)");
        
        if (rec_form.first)
        {
            vstr.push_back("const auto coords = distributor.coordinates();");
        }
        else
        {
            vstr.push_back("const auto r_c = distributor.coordinates()[0];");
        }
    }
    
    if ((integral.integrand().name() == "G(r)"))
    {
        vstr.push_back("// intialize external Gaussian(s)");
        
        vstr.push_back("const auto exgtos = distributor.data();");
    }
    
    if ((integral.integrand().name() == "A"))
    {
        vstr.push_back("// intialize external charge(s)");
        
        if (rec_form.first)
        {
            vstr.push_back("const auto charges = distributor.data();");
        }
        else
        {
            vstr.push_back("const auto charge = distributor.data()[0];");
        }
    }
    
    if ((integral.integrand().name() == "AG"))
    {
        const auto iorder = integral.integrand().shape().order();

        if (rec_form.first)
        {
            if (iorder == 1)
            {
                vstr.push_back("// intialize external dipoles data");

                vstr.push_back("const auto dipoles = distributor.data();");
            }

            if (iorder == 2)
            {
                vstr.push_back("// intialize external quadrupoles data");

                vstr.push_back("const auto quadrupoles = distributor.data();");
            }
//
//            if (iorder == 3)
//            {
//                vstr.push_back("// intialize external octupoles data");
//
//                vstr.push_back("const auto quadrupoles = distributor->data();");
//            }
//
//            vstr.push_back("const auto coords_x = distributor->coordinates_x();");
//
//            vstr.push_back("const auto coords_y = distributor->coordinates_y();");
//
//            vstr.push_back("const auto coords_z = distributor->coordinates_z();");
       }
        else
        {
//            if (iorder == 1)
//            {
//                vstr.push_back("// intialize external dipole data");
//
//                vstr.push_back("const auto dipole = distributor->data();");
//            }
//
//            if (iorder == 2)
//            {
//                vstr.push_back("// intialize external quadrupole data");
//
//                vstr.push_back("const auto quadrupole = distributor->data();");
//            }
//
//            if (iorder == 3)
//            {
//                vstr.push_back("// intialize external octupole data");
//
//                vstr.push_back("const auto quadrupole = distributor->data();");
//            }
//
//            vstr.push_back("const double coord_x = distributor->coordinates_x()[0];");
//
//            vstr.push_back("const double coord_y = distributor->coordinates_y()[0];");
//
//            vstr.push_back("const double coord_z = distributor->coordinates_z()[0];");
        }
   }
    
    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_gtos_def() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// intialize GTOs data on bra side");
        
    vstr.push_back("const auto bra_gto_coords = bra_gto_block.coordinates();");
        
    vstr.push_back("const auto bra_gto_exps = bra_gto_block.exponents();");
        
    vstr.push_back("const auto bra_gto_norms = bra_gto_block.normalization_factors();");
       
    vstr.push_back("const auto bra_gto_indices = bra_gto_block.orbital_indices();");
        
    vstr.push_back("const auto bra_ncgtos = bra_gto_block.number_of_basis_functions();");

    vstr.push_back("const auto bra_npgtos = bra_gto_block.number_of_primitives();");
        
    vstr.push_back("// intialize GTOs data on ket side");
        
    vstr.push_back("const auto ket_gto_coords = ket_gto_block.coordinates();");
        
    vstr.push_back("const auto ket_gto_exps = ket_gto_block.exponents();");
        
    vstr.push_back("const auto ket_gto_norms = ket_gto_block.normalization_factors();");
       
    vstr.push_back("const auto ket_gto_indices = ket_gto_block.orbital_indices();");

    vstr.push_back("const auto ket_npgtos = ket_gto_block.number_of_primitives();");
    
    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_ket_variables_def(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned 2D arrays for ket side");
    
    size_t nelems = 8;
    
    if (_need_center_p(integral)) nelems += 3;
    
    if (t2c::get_effective_order(integral, 0) > 0) nelems += 3;
    
    if (t2c::get_effective_order(integral, 1) > 0) nelems += 3;
    
    if (_need_distances_pc(integral)) nelems += 3;
    
    vstr.push_back("CSimdArray<double> factors(" + std::to_string(nelems) +  ", ket_npgtos);");

    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_buffers_def(const SI2CIntegrals& integrals,
                                    const I2CIntegral&   integral,
                                    const std::array<int, 3>& geom_drvs) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned primitive integrals");
    
    size_t icomps = 0;
    
    for (const auto& tint : integrals)
    {
        icomps += tint.components<T1CPair, T1CPair>().size();
    }
    
    if (_need_geom_drvs(geom_drvs))
    {
        icomps += integral.components<T1CPair, T1CPair>().size();
    }
    
    auto label = "CSimdArray<double> pbuffer(" + std::to_string(icomps) + ", ket_npgtos);";
    
    vstr.push_back(label);
    
    vstr.push_back("// allocate aligned contracted integrals");
    
    icomps = integral.components<T1CPair, T1CPair>().size();
    
    label = "CSimdArray<double> cbuffer(" + std::to_string(icomps) + ", 1);";
    
    vstr.push_back(label);
    
    if ((integral[0] + integral[1]) > 0)
    {
        const auto angpair = std::array<int, 2>({integral[0], integral[1]});
        
        icomps = t2c::number_of_spherical_components(angpair);
        
        icomps *= integral.integrand().components().size();
        
        if (const auto prefixes = integral.prefixes(); !prefixes.empty())
        {
            icomps *= make_components<OperatorComponent>(prefixes).size();
        }
        
        label = "CSimdArray<double> sbuffer(" + std::to_string(icomps) + ", 1);";
        
        vstr.push_back(label);
    }
    
    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_boys_function_def(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (_need_boys_func(integral))
    {
        auto order = integral[0] + integral[1] + integral.integrand().shape().order();
        
        if (auto prefixes = integral.prefixes(); !prefixes.empty())
        {
            for (const auto& prefix : prefixes)
            {
                order += prefix.shape().order();
            }
        }
        
        vstr.push_back("// setup Boys function data");
        
        vstr.push_back("const CBoysFunc<" + std::to_string(order) + "> bf_table;");

        vstr.push_back("CSimdArray<double> bf_data(" + std::to_string(order + 2) + ", ket_npgtos);");
    }
    
    return vstr;
}

void
T2CFuncBodyDriver::_add_loop_start(      VCodeLines&  lines,
                                   const I2CIntegral& integral) const
{
    lines.push_back({1, 0, 2, "// set up ket partitioning"});

    lines.push_back({1, 0, 2, "const auto ket_dim = ket_indices.second - ket_indices.first;"});

    lines.push_back({1, 0, 2, "const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());"});

    lines.push_back({1, 0, 1, "for (size_t i = 0; i < ket_blocks; i++)"});
                    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), ket_indices.first);"});

    lines.push_back({2, 0, 2, "factors.load(ket_gto_exps, ket_range, 0, ket_npgtos);"});

    lines.push_back({2, 0, 2, "factors.load(ket_gto_norms, ket_range, 1, ket_npgtos);"});

    lines.push_back({2, 0, 2, "factors.replicate_points(ket_gto_coords, ket_range, 2, ket_npgtos);"});
    
    lines.push_back({2, 0, 2, "// set up active SIMD width"});
    
    lines.push_back({2, 0, 2, "const auto ket_width = ket_range.second - ket_range.first;"});
    
    if ((integral[0] + integral[1]) > 0)
    {
        lines.push_back({2, 0, 2, "sbuffer.set_active_width(ket_width);"});
    }
    
    lines.push_back({2, 0, 2, "cbuffer.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "pbuffer.set_active_width(ket_width);"});
    
    if (_need_boys_func(integral))
    {
        lines.push_back({2, 0, 2, "bf_data.set_active_width(ket_width);"});
    }

    lines.push_back({2, 0, 2, "// loop over contracted basis functions on bra side"});
    
    lines.push_back({2, 0, 1, "for (auto j = bra_indices.first; j < bra_indices.second; j++)"});
    
    lines.push_back({2, 0, 1, "{"});
        
    lines.push_back({3, 0, 2, "cbuffer.zero();"});
    
    if ((integral[0] + integral[1]) > 0)
    {
        lines.push_back({3, 0, 2, "sbuffer.zero();"});
    }

    lines.push_back({3, 0, 2, "const auto r_a = bra_gto_coords[j];"});

    lines.push_back({3, 0, 2, "t2cfunc::comp_distances_ab(factors, 5, 2, r_a);"});
}

void
T2CFuncBodyDriver::_add_loop_end(      VCodeLines&  lines,
                                 const I2CIntegral& integral,
                                 const std::pair<bool, bool>& rec_form) const
{
    std::string label;
    
    if ((integral[0] + integral[1]) > 0)
    {
        label = "t2cfunc::transform<"  + std::to_string(integral[0]);
            
        label += ", " + std::to_string(integral[1]) + ">(sbuffer, cbuffer);";
            
        lines.push_back({3, 0, 2, label});
    }
    
    if ((integral[0] + integral[1]) > 0)
    {
        label = "distributor.distribute(sbuffer, ";
    }
    else
    {
        label = "distributor.distribute(cbuffer, ";
    }
    
    label += "bra_gto_indices, ket_gto_indices, ";
            
    label += std::to_string(integral[0]) + ", ";
            
    label += std::to_string(integral[1]) + ", ";
    
    label += "j, ket_range, bra_eq_ket);";
    
    lines.push_back({3, 0, 1, label});
   
    lines.push_back({2, 0, 1, "}"});
    
    lines.push_back({1, 0, 1, "}"});
}

void
T2CFuncBodyDriver::_add_ket_loop_start(      VCodeLines&            lines,
                                       const I2CIntegral&           integral,
                                       const std::pair<bool, bool>& rec_form) const
{
    lines.push_back({3, 0, 1, "for (size_t k = 0; k < bra_npgtos; k++)"});
    
    lines.push_back({3, 0, 1, "{"});
    
    lines.push_back({4, 0, 2, "const auto a_exp = bra_gto_exps[k * bra_ncgtos + j];"});

    lines.push_back({4, 0, 2, "const auto a_norm = bra_gto_norms[k * bra_ncgtos + j];"});

    if (_need_center_p(integral))
    {
        lines.push_back({4, 0, 2, "t2cfunc::comp_coordinates_p(factors, 8, 2, r_a, a_exp);"});
    }

    if ((t2c::get_effective_order(integral, 0) > 0) && (integral.integrand().name() != "G(r)"))
    {
        
        const auto label = std::to_string(_get_index_pa(integral));
        
        if (_need_center_p(integral))
        {
            lines.push_back({4, 0, 2, "t2cfunc::comp_distances_pa_from_p(factors, " + label + " , 8, r_a);"});
        }
        else
        {
            lines.push_back({4, 0, 2, "t2cfunc::comp_distances_pa(factors, " + label +  ", 5, a_exp);"});
        }
    }

    if ((t2c::get_effective_order(integral, 1) > 0) && (integral.integrand().name() != "G(r)"))
    {
        const auto label = std::to_string(_get_index_pb(integral));
        
        if (_need_center_p(integral))
        {
            lines.push_back({4, 0, 2, "t2cfunc::comp_distances_pb_from_p(factors, " +  label + " , 8, 2);"});
        }
        else
        {
            lines.push_back({4, 0, 2, "t2cfunc::comp_distances_pb(factors, " +  label + ", 5, a_exp);"});
        }
    }

    if (_need_distances_pc(integral) && !rec_form.first)
    {
        const auto label = std::to_string(_get_index_pc(integral));
        
        lines.push_back({4, 0, 2, "t2cfunc::comp_distances_pc(factors, " +  label + ", 8, r_c);"});
    }
}

void
T2CFuncBodyDriver::_add_ket_loop_end(      VCodeLines&            lines,
                                     const SI2CIntegrals&         integrals,
                                     const I2CIntegral&           integral,
                                     const std::pair<bool, bool>& rec_form) const
{
    if (!rec_form.first)
    {
        std::string label = "t2cfunc::reduce(cbuffer, pbuffer, ";
        
        if (integral.is_simple())
        {
            label += std::to_string(_get_position(integral, integrals)) + ", ";
        }
        else
        {
            size_t icomps = 0;
            
            for (const auto& tint : integrals)
            {
                icomps += tint.components<T1CPair, T1CPair>().size();
            }
            
            label += std::to_string(icomps)  + ", ";
        }
        
        label += "ket_width, ket_npgtos);";
        
        lines.push_back({4, 0, 1, label});
    }
    
    lines.push_back({3, 0, 2, "}"});
}

void
T2CFuncBodyDriver::_add_sum_loop_start(      VCodeLines&            lines,
                                       const I2CIntegral&           integral,
                                       const std::pair<bool, bool>& rec_form,
                                       const bool                   use_rs) const
{
    if (rec_form.first)
    {
        const auto integrand = integral.integrand();
        
        if (integrand.name() == "G(r)")
        {
            lines.push_back({4, 0, 2, "const size_t npoints = coords.size();"});
                
            lines.push_back({4, 0, 1, "for (size_t l = 0; l < npoints; l++)"});
        }
        else
        {
            lines.push_back({4, 0, 1, "for (size_t l = 0; l < coords.size(); l++)"});
        }
        
        lines.push_back({4, 0, 1, "{"});
       
        if (_need_distances_pc(integral))
        {
            const auto label = std::to_string(_get_index_pc(integral));
        
            lines.push_back({5, 0, 2, "t2cfunc::comp_distances_pc(factors, " + label + ", 8, coords[l]);"});
        }
        
        if (_need_distances_ga(integral))
        {
            const auto label = std::to_string(_get_index_ga(integral));
        
            lines.push_back({5, 0, 2, "t2cfunc::comp_distances_ga(factors, " + label + ", 8, r_a, coords[l], a_exp, exgtos[l]);"});
        }
        
        if (_need_distances_gb(integral))
        {
            const auto label = std::to_string(_get_index_gb(integral));
        
            lines.push_back({5, 0, 2, "t2cfunc::comp_distances_gb(factors, " + label + ", 8, 2, coords[l], a_exp, exgtos[l]);"});
        }
        
        if (_need_boys_func(integral))
        {
            auto order = integral[0] + integral[1] + integral.integrand().shape().order();
            
            const auto label = std::to_string(_get_index_pc(integral));
            
            if (use_rs)
            {
                lines.push_back({5, 0, 2, "t2cfunc::comp_boys_args(bf_data, " + std::to_string(order + 1) + ", factors, " + label + ", a_exp, omegas[l]);"});
                
                lines.push_back({5, 0, 2, "bf_table.compute(bf_data, 0, " + std::to_string(order + 1) + ", factors, a_exp, omegas[l]);"});
            }
            else
            {
                lines.push_back({5, 0, 2, "t2cfunc::comp_boys_args(bf_data, " + std::to_string(order + 1) + ", factors, " + label + ", a_exp);"});
                
                lines.push_back({5, 0, 2, "bf_table.compute(bf_data, 0, " + std::to_string(order + 1) + ");"});
            }
        }
    }
}

void
T2CFuncBodyDriver::_add_sum_loop_end(      VCodeLines&            lines,
                                     const SI2CIntegrals&         integrals,
                                     const I2CIntegral&           integral,
                                     const std::pair<bool, bool>& rec_form) const
{
    if (rec_form.first)
    {
        std::string label = "t2cfunc::reduce(cbuffer, pbuffer, ";
        
        if (integral.is_simple())
        {
            label += std::to_string(_get_position(integral, integrals)) + ", ";
        }
        else
        {
            size_t icomps = 0;
            
            for (const auto& tint : integrals)
            {
                icomps += tint.components<T1CPair, T1CPair>().size();
            }
            
            label += std::to_string(icomps)  + ", ";
        }
        
        if (integral.integrand().name() == "A")
        {
            label += "charges[l], ";
        }
        
        if (integral.integrand().name() == "AG")
        {
            const auto iorder = integral.integrand().shape().order();
            
            if (iorder == 1) label += "dipoles, 3, l, ";
            
            if (iorder == 2) label += "quadrupoles, 6, l, ";
            
            if (iorder == 3) label += "octupoles, 10, l, ";
            
            if (iorder == 4) label += "hexadecapoles, 15, l, ";
        }
        
        label += "ket_width, ket_npgtos);";
    
        lines.push_back({5, 0, 1, label});
        
        lines.push_back({4, 0, 1, "}"});
    }
}

void
T2CFuncBodyDriver::_add_auxilary_integrals(      VCodeLines&            lines,
                                           const SI2CIntegrals&         integrals,
                                           const I2CIntegral&           integral,
                                           const std::pair<bool, bool>& rec_form,
                                           const bool                   in_sum_loop) const
{
    const auto spacer = (rec_form.first) ? 5 : 4;
    
    for (const auto& tint : integrals)
    {
        if (!tint.is_simple()) continue;
        
        if ((tint[0] == 0) && (tint[1] == 0))
        {
            if ((tint.integrand().name() == "1") && (!in_sum_loop))
            {
                lines.push_back({4, 0, 2, "ovlrec::comp_prim_overlap_ss(pbuffer, " + std::to_string(_get_position(tint, integrals)) + ", factors, a_exp, a_norm);"});
            }
            
            if ((tint.integrand().name() == "T") && (!in_sum_loop))
            {
                const auto sint = tint.replace(Operator("1"));
                
                const auto label = std::to_string(_get_position(tint, integrals)) + ", " + std::to_string(_get_position(sint, integrals));
                
                lines.push_back({4, 0, 2, "kinrec::comp_prim_kinetic_energy_ss(pbuffer, " + label + ", factors, a_exp);"});
            }

            if ((tint.integrand().name() == "r") && (!in_sum_loop))
            {
                const auto sint = tint.replace(Operator("1"));
                
                const auto label = std::to_string(_get_position(tint, integrals)) + ", " + std::to_string(_get_position(sint, integrals));
                
                lines.push_back({4, 0, 2, "diprec::comp_prim_electric_dipole_momentum_ss(pbuffer, " + label + ", factors, " +  std::to_string(_get_index_pc(integral)) + ");"});
            }

            if ((tint.integrand().name() == "p") && (!in_sum_loop))
            {
                lines.push_back({3, 0, 2, "linmomrec::comp_prim_dipole_ss(prim_buffer_dip_ss, prim_buffer_ovl_ss);"});
            }
            
            if ((tint.integrand().name() == "A"))
            {
                const auto sint = tint.replace(Operator("1"));
                
                if (rec_form.first && in_sum_loop)
                {
                    const auto label = std::to_string(_get_position(tint, integrals)) + ", " + std::to_string(_get_position(sint, integrals));
                    
                    lines.push_back({spacer, 0, 2, "npotrec::comp_prim_nuclear_potential_ss(pbuffer, " + label + ", bf_data, " + std::to_string(tint.order()) + ", factors, a_exp);"});
                }
            }
            
            if ((tint.integrand().name() == "G(r)"))
            {
                const auto sint = tint.replace(Operator("1"));
                
                if (rec_form.first && in_sum_loop)
                {
                    const auto label = std::to_string(_get_position(tint, integrals)) + ", " + std::to_string(_get_position(sint, integrals));
                    
                    lines.push_back({spacer, 0, 2, "t3ovlrec::comp_prim_overlap_ss(pbuffer, " + label + ", factors, " + std::to_string(_get_index_pc(integral)) + ", a_exp, exgtos[l], exgtos[npoints + l]);"});
                }
            }
            
            
            if ((tint.integrand().name() == "AG"))
            {
                const auto iorder = tint.integrand().shape().order();
                
                if (rec_form.first && in_sum_loop)
                {
                    if (iorder == 1)
                    {
                        auto xint = tint.replace(Operator("A"));
                        
                        auto label = std::to_string(_get_position(tint, integrals)) + ", ";
                        
                        xint.set_order(tint.order() + 1);
                        
                        label += std::to_string(_get_position(xint, integrals));
                        
                        lines.push_back({spacer, 0, 2, "npotrec::comp_prim_nuclear_potential_geom_010_ss(pbuffer, " + label + ", factors, " +  std::to_string(_get_index_pc(integral)) + ", a_exp);"});
                    }
                    
                    if (iorder == 2)
                    {
                        auto xint = tint.replace(Operator("A"));
                        
                        auto label = std::to_string(_get_position(tint, integrals)) + ", ";
                        
                        xint.set_order(tint.order() + 1);
                        
                        label += std::to_string(_get_position(xint, integrals)) + ", ";
                        
                        xint.set_order(tint.order() + 2);
                        
                        label += std::to_string(_get_position(xint, integrals));
                        
                        lines.push_back({spacer, 0, 2, "npotrec::comp_prim_nuclear_potential_geom_020_ss(pbuffer, " + label + ", factors, " +  std::to_string(_get_index_pc(integral)) + ", a_exp);"});
                    }
                }
            }
            
            // TODO: other integrals...
        }
    }
}

// May need to change this for new integral situations
void
T2CFuncBodyDriver::_add_call_tree(      VCodeLines&            lines,
                                  const SI2CIntegrals&         integrals,
                                  const I2CIntegral&           integral,
                                  const std::pair<bool, bool>& rec_form) const
{
    const int spacer = (rec_form.first) ? 5 : 4;
    
    for (const auto& tint : integrals)
    {
        if (!tint.is_simple()) continue;
        
        if ((tint[0] != 0) || (tint[1] != 0))
        {
            const auto name = t2c::prim_compute_func_name(tint);
            
            auto label = t2c::namespace_label(tint) + "::" + name + "(pbuffer, ";
            
            label += _get_arguments(tint, integrals);
            
            label += "factors, "; 
            
            if (tint[0] > 0)
            {
                label += std::to_string(_get_index_pa(integral)) + ", ";
            }
            
            if ((tint[1] > 0) && (tint[0] == 0))
            {
                label += std::to_string(_get_index_pb(integral)) + ", ";
            }
            
            if (_need_distances_pc_in_call_tree(tint))
            {
                label += std::to_string(_get_index_pc(integral)) + ", ";
            }
           
            if (_need_exponents(tint))
            {
                label += "a_exp";
            }
            else
            {
                label.pop_back();
                
                label.pop_back();
            }

            label += ");";
            
            lines.push_back({spacer, 0, 2, label});
        }
    }
}

void
T2CFuncBodyDriver::_add_geom_call_tree(      VCodeLines&            lines,
                                       const SI2CIntegrals&         geom_integrals,
                                       const SI2CIntegrals&         vrr_integrals,
                                       const I2CIntegral&           integral,
                                       const std::array<int, 3>&    geom_drvs,
                                       const std::pair<bool, bool>& rec_form) const
{
    if (_need_geom_drvs(geom_drvs))
    {
        const int spacer = (rec_form.first) ? 5 : 4;
        
        const auto tint = integral.replace(Operator("R")); 
        
        const auto name = t2c::prim_compute_func_name(tint);
        
        auto label = "t2cgeom::" + name + "(pbuffer, ";
        
        size_t icomps = 0;
        
        for (const auto& tint : vrr_integrals)
        {
            icomps += tint.components<T1CPair, T1CPair>().size();
        }
        
        label += std::to_string(icomps)  + ", ";
        
        for (auto cint : geom_integrals)
        {
            label += std::to_string(_get_position(cint, vrr_integrals)) + ", ";
        }
        
        label += std::to_string(integral.integrand().shape().components().size()) + ", ";
        
        if (geom_drvs[2] == 0)
        {
            label += std::to_string(Tensor(integral[1]).components().size()) + ", ";
        }
        
        if (geom_drvs[2] > 0) label += "factors, ";
        
        label += "a_exp);";
        
        lines.push_back({spacer, 0, 2, label});
    }
}

std::string
T2CFuncBodyDriver::_get_arguments(const I2CIntegral& integral) const
{
    std::string label = t2c::get_buffer_label(integral, {"prim"}) + ", ";;
   
    for (const auto& tint : t2c::get_integrals(integral))
    {
        label += t2c::get_buffer_label(tint, {"prim"}) + ", ";
    }
    
    return label;
}

std::string
T2CFuncBodyDriver::_get_arguments(const I2CIntegral&   integral,
                                  const SI2CIntegrals& integrals) const
{
    auto label = std::to_string(_get_position(integral, integrals)) + ", ";
    
    for (const auto& tint : t2c::get_integrals(integral))
    {
        label += std::to_string(_get_position(tint, integrals)) + ", ";
    }
    
    return label;
}

size_t
T2CFuncBodyDriver::_get_position(const I2CIntegral&   integral,
                                 const SI2CIntegrals& integrals) const
{
    size_t pos = 0;
    
    for (const auto& tint : integrals)
    {
        if (tint == integral) return pos;
        
        pos += tint.components<T1CPair, T1CPair>().size();
    }
    
    return 0;
}

bool
T2CFuncBodyDriver::_need_center_p(const I2CIntegral& integral) const
{
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "r") return true;
    
    if (integrand.name() == "G(r)") return true;
    
    if (integrand.name() == "A") return true;
    
    if (integrand.name() == "AG") return true;
    
    return false;
}

bool
T2CFuncBodyDriver::_need_distances_pc(const I2CIntegral& integral) const
{
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "r") return true;
    
    if (integrand.name() == "G(r)") return true;
    
    if (integrand.name() == "A") return true;
    
    if (integrand.name() == "AG") return true;
    
    return false;
}

bool
T2CFuncBodyDriver::_need_distances_pc_in_call_tree(const I2CIntegral& integral) const
{
    const auto integrand = integral.integrand();
    
    if (((integrand.name() == "r") && ((integral[0] + integral[1]) == 0))) return true;
    
    if (integrand.name() == "A") return true;
    
    if (integrand.name() == "AG") return true;
    
    return false;
}

bool
T2CFuncBodyDriver::_need_distances_pa(const I2CIntegral& integral) const
{
    if (integral.integrand().name() == "G(r)") return false;
    
    if (integral.is_simple())
    {
        return integral[0] > 0;
    }
    else
    {
        return (integral[0] + (integral.prefixes())[0].shape().order()) > 0;
    }
}

bool
T2CFuncBodyDriver::_need_distances_pb(const I2CIntegral& integral) const
{
    if (integral.integrand().name() == "G(r)") return false;
    
    if (integral.is_simple())
    {
        return integral[1] > 0;
    }
    else
    {
        if (auto prefixes = integral.prefixes(); prefixes.size() == 2)
        {
            return (integral[1] + prefixes[1].shape().order()) > 0;
        }
        else
        {
            return integral[1] > 0;
        }
    }
}

bool
T2CFuncBodyDriver::_need_distances_ga(const I2CIntegral& integral) const
{
    if (integral.integrand().name() != "G(r)") return false;
    
    if (integral.is_simple())
    {
        return integral[0] > 0;
    }
    else
    {
        return (integral[0] + (integral.prefixes())[0].shape().order()) > 0;
    }
}

bool
T2CFuncBodyDriver::_need_distances_gb(const I2CIntegral& integral) const
{
    if (integral.integrand().name() != "G(r)") return false;
    
    if (integral.is_simple())
    {
        return integral[1] > 0;
    }
    else
    {
        if (auto prefixes = integral.prefixes(); prefixes.size() == 2)
        {
            return (integral[1] + prefixes[1].shape().order()) > 0;
        }
        else
        {
            return integral[1] > 0;
        }
    }
}

bool
T2CFuncBodyDriver::_need_exponents(const I2CIntegral& integral) const
{
    if (integral.integrand().name() == "T") return true;
    
    if ((integral.integrand().name() == "r") && ((integral[0] + integral[1]) > 0)) return true; 
    
    int order = 0;
    
    for (const auto& prefix : integral.prefixes())
    {
        order += prefix.shape().order();
    }
    
    return (order + integral[0] + integral[1]) > 1;
}

bool
T2CFuncBodyDriver::_need_boys_func(const I2CIntegral& integral) const
{
    const auto integrand = integral.integrand();
        
    if (integrand.name() == "A") return true;
    
    if (integrand.name() == "AG") return true;
    
    return false;
}

bool
T2CFuncBodyDriver::_need_external_coords(const I2CIntegral& integral) const
{
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "A") return true;
    
    if (integrand.name() == "AG") return true;
    
    if (integrand.name() == "r") return true;
    
    if (integrand.name() == "G(r)") return true;
    
    return false;
}

bool
T2CFuncBodyDriver::_need_geom_drvs(const std::array<int, 3>& geom_drvs) const
{
    return (geom_drvs[0] + geom_drvs[2]) > 0;
}

int
T2CFuncBodyDriver::_get_index_pa(const I2CIntegral& integral) const
{
    if (_need_center_p(integral))
    {
        return 11;
    }
    else
    {
        return 8;
    }
}

int
T2CFuncBodyDriver::_get_index_pb(const I2CIntegral& integral) const
{
    if (_need_distances_pa(integral))
    {
       return _get_index_pa(integral) + 3;
    }
    else
    {
        return _get_index_pa(integral);
    }
}

int
T2CFuncBodyDriver::_get_index_pc(const I2CIntegral& integral) const
{
    if (_need_distances_pb(integral))
    {
       return _get_index_pb(integral) + 3;
    }
    else
    {
        if (_need_distances_pa(integral))
        {
            return _get_index_pa(integral) + 3;
        }
        else
        {
            return 11;
        }
    }
}

int
T2CFuncBodyDriver::_get_index_ga(const I2CIntegral& integral) const
{
    return _get_index_pc(integral) + 3;
}

int
T2CFuncBodyDriver::_get_index_gb(const I2CIntegral& integral) const
{
    if (_need_distances_ga(integral))
    {
       return _get_index_ga(integral) + 3;
    }
    else
    {
        return _get_index_pc(integral) + 3;
    }
}
