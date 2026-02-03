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

#include "t2c_ecp_body.hpp"

#include "t2c_utils.hpp"

void
T2CECPFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                      const SI2CIntegrals& hrr_integrals,
                                      const SI2CIntegrals& vrr_integrals,
                                      const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_gtos_def())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_ket_variables_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    auto ctints = _filter_contracted(vrr_integrals, integral);
    
    for (const auto& tint : hrr_integrals)
    {
        ctints.insert(tint);
    }
    
    for (const auto& label : _get_buffers_def(ctints, vrr_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_loop_start(lines, integral);
    
    _add_ket_loop_start(lines, integral);
//    
//    _add_auxilary_integrals(lines, vrr_integrals, integral, rec_form, false);
//    
//    _add_sum_loop_start(lines, integral, rec_form, use_rs);
//    
//    _add_auxilary_integrals(lines, vrr_integrals, integral, rec_form, true);
//    
//    _add_call_tree(lines, vrr_integrals, integral, rec_form);
//    
//    _add_geom_call_tree(lines, geom_integrals, vrr_integrals, integral, geom_drvs, rec_form);
//    
//    _add_sum_loop_end(lines, vrr_integrals, integral, rec_form);
//    
    _add_ket_loop_end(lines, vrr_integrals, integral);
    
    _add_loop_end(lines, integral);
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CECPFuncBodyDriver::_get_gtos_def() const
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
    
    vstr.push_back("// intialize basic ECP data");
    
    vstr.push_back("const auto ecp_nppt = ecp_potential.number_of_primitive_potentials();");
    
    vstr.push_back("const auto ecp_exps = ecp_potential.get_exponents();");
        
    vstr.push_back("const auto ecp_facts = ecp_potential.get_factors();");
    
    return vstr;
}

std::vector<std::string>
T2CECPFuncBodyDriver::_get_ket_variables_def(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned 2D arrays for ket side");
    
    size_t nelems = 8;
    
    if (_need_distances_ra(integral)) nelems += 3;
    
    if (_need_distances_rb(integral)) nelems += 3;
    
    vstr.push_back("CSimdArray<double> pfactors(" + std::to_string(nelems) +  ", ket_npgtos);");
    
    vstr.push_back("CSimdArray<double> cfactors(6, 1);");

    return vstr;
}

std::vector<std::string>
T2CECPFuncBodyDriver::_get_buffers_def(const SI2CIntegrals& hrr_integrals,
                                       const SI2CIntegrals& vrr_integrals,
                                       const I2CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned primitive integrals");
    
    size_t icomps = 0;
    
    for (const auto& tint : vrr_integrals)
    {
        icomps += tint.components<T1CPair, T1CPair>().size();
    }
    
    auto label = "CSimdArray<double> pbuffer(" + std::to_string(icomps) + ", ket_npgtos);";
    
    vstr.push_back(label);
    
    vstr.push_back("// allocate aligned contracted integrals");
    
    icomps = 0;
    
    for (const auto& tint : hrr_integrals)
    {
        icomps += tint.components<T1CPair, T1CPair>().size();
    }
    
    label = "CSimdArray<double> cbuffer(" + std::to_string(icomps) + ", 1);";
    
    vstr.push_back(label);
    
    if ((integral[0] + integral[1]) > 0)
    {
        const auto angpair = std::array<int, 2>({integral[0], integral[1]});
        
        icomps = t2c::number_of_spherical_components(angpair);
        
        label = "CSimdArray<double> sbuffer(" + std::to_string(icomps) + ", 1);";
        
        vstr.push_back(label);
    }
    
    return vstr;
}

SI2CIntegrals
T2CECPFuncBodyDriver::_filter_contracted(const SI2CIntegrals& integrals,
                                         const I2CIntegral&   integral) const
{
    SI2CIntegrals tints;
    
    if (integral[0] > integral[1])
    {
        for (const auto &tint : integrals)
        {
            if (tint[0] >= integral[0]) tints.insert(tint);
        }
    } else {
        for (const auto &tint : integrals)
        {
            if (tint[1] >= integral[1]) tints.insert(tint);
        }
    }
    
    return tints;
}

void
T2CECPFuncBodyDriver::_add_loop_start(      VCodeLines&  lines,
                                      const I2CIntegral& integral) const
{
    lines.push_back({1, 0, 2, "// set up ket partitioning"});

    lines.push_back({1, 0, 2, "const auto ket_dim = ket_indices.second - ket_indices.first;"});

    lines.push_back({1, 0, 2, "const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());"});

    lines.push_back({1, 0, 1, "for (size_t i = 0; i < ket_blocks; i++)"});
                    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), ket_indices.first);"});

    lines.push_back({2, 0, 2, "pfactors.load(ket_gto_exps, ket_range, 0, ket_npgtos);"});

    lines.push_back({2, 0, 2, "pfactors.load(ket_gto_norms, ket_range, 1, ket_npgtos);"});

    lines.push_back({2, 0, 2, "pfactors.replicate_points(ket_gto_coords, ket_range, 2, ket_npgtos);"});
    
    lines.push_back({2, 0, 2, "cfactors.replicate_points(ket_gto_coords, ket_range, 0, 1);"});
    
    lines.push_back({2, 0, 2, "// set up active SIMD width"});
    
    lines.push_back({2, 0, 2, "const auto ket_width = ket_range.second - ket_range.first;"});
    
    lines.push_back({2, 0, 2, "sbuffer.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "cbuffer.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "pbuffer.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "// loop over contracted basis functions on bra side"});
    
    lines.push_back({2, 0, 1, "for (auto j = bra_indices.first; j < bra_indices.second; j++)"});
    
    lines.push_back({2, 0, 1, "{"});
        
    lines.push_back({3, 0, 2, "cbuffer.zero();"});
    
    lines.push_back({3, 0, 2, "sbuffer.zero();"});

    lines.push_back({3, 0, 2, "const auto r_a = bra_gto_coords[j];"});

    lines.push_back({3, 0, 2, "t2cfunc::comp_distances_ab(cfactors, 3, 0, r_a);"});
}

void
T2CECPFuncBodyDriver::_add_loop_end(      VCodeLines&  lines,
                                    const I2CIntegral& integral) const
{
    std::string label;

    label = "t2cfunc::transform<"  + std::to_string(integral[0]);
            
    label += ", " + std::to_string(integral[1]) + ">(sbuffer, cbuffer);";
            
    lines.push_back({3, 0, 2, label});
    
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
T2CECPFuncBodyDriver::_add_ket_loop_start(      VCodeLines&  lines,
                                          const I2CIntegral& integral) const
{
    lines.push_back({3, 0, 1, "for (size_t k = 0; k < bra_npgtos; k++)"});
    
    lines.push_back({3, 0, 1, "{"});
    
    lines.push_back({4, 0, 2, "const auto a_exp = bra_gto_exps[k * bra_ncgtos + j];"});

    lines.push_back({4, 0, 2, "const auto a_norm = bra_gto_norms[k * bra_ncgtos + j];"});
    
    lines.push_back({4, 0, 1, "for (size_t l = 0; l < ecp_nppt; l++)"});
    
    lines.push_back({4, 0, 1, "{"});
    
    lines.push_back({5, 0, 2, "const auto c_exp = ecp_exps[l];"});

    lines.push_back({5, 0, 2, "const auto c_norm = ecp_facts[l];"});
    
    lines.push_back({5, 0, 2, "t2cfunc::comp_coordinates_r(factors, 5, 2, r_a, a_exp, c_exp);"});
    
    
    

//    if (_need_center_p(integral))
//    {
//        lines.push_back({4, 0, 2, "t2cfunc::comp_coordinates_p(factors, 8, 2, r_a, a_exp);"});
//    }
//
//    if (
//        (t2c::get_effective_order(integral, 0) > 0) &&
//        (integral.integrand().name() != "G(r)")     &&
//        (integral.integrand().name() != "GX(r)")    &&
//        (integral.integrand().name() != "GR2(r)")   &&
//        (integral.integrand().name() != "GR.R2(r)")
//        )
//    {
//        const auto label = std::to_string(_get_index_pa(integral));
//        
//        if (_need_center_p(integral))
//        {
//            lines.push_back({4, 0, 2, "t2cfunc::comp_distances_pa_from_p(factors, " + label + " , 8, r_a);"});
//        }
//        else
//        {
//            lines.push_back({4, 0, 2, "t2cfunc::comp_distances_pa(factors, " + label +  ", 5, a_exp);"});
//        }
//    }
//
//    if (
//        (t2c::get_effective_order(integral, 1) > 0) &&
//        (integral.integrand().name() != "G(r)")     &&
//        (integral.integrand().name() != "GX(r)")    &&
//        (integral.integrand().name() != "GR2(r)")   &&
//        (integral.integrand().name() != "GR.R2(r)")
//        )
//    {
//        const auto label = std::to_string(_get_index_pb(integral));
//        
//        if (_need_center_p(integral))
//        {
//            lines.push_back({4, 0, 2, "t2cfunc::comp_distances_pb_from_p(factors, " +  label + " , 8, 2);"});
//        }
//        else
//        {
//            lines.push_back({4, 0, 2, "t2cfunc::comp_distances_pb(factors, " +  label + ", 5, a_exp);"});
//        }
//    }
//
//    if (_need_distances_pc(integral) && !rec_form.first)
//    {
//        const auto label = std::to_string(_get_index_pc(integral));
//        
//        lines.push_back({4, 0, 2, "t2cfunc::comp_distances_pc(factors, " +  label + ", 8, r_c);"});
//    }
//    
//    if (_need_boys_func(integral))
//    {
//        if (!rec_form.first)
//        {
//            auto order = integral[0] + integral[1] + integral.integrand().shape().order();
//            
//            lines.push_back({4, 0, 2, "t2cfunc::comp_boys_args_with_rho(bf_data, " + std::to_string(order + 1) + ", factors, 5, a_exp);"});
//            
//            lines.push_back({4, 0, 2, "bf_table.compute(bf_data, 0, " + std::to_string(order + 1) + ");"});
//        }
//    }
}

void
T2CECPFuncBodyDriver::_add_ket_loop_end(      VCodeLines&            lines,
                                        const SI2CIntegrals&         integrals,
                                        const I2CIntegral&           integral) const
{
//    if (!rec_form.first)
//    {
//        std::string label = "t2cfunc::reduce(cbuffer, pbuffer, ";
//        
//        if (integral.is_simple())
//        {
//            label += std::to_string(_get_position(integral, integrals)) + ", ";
//        }
//        else
//        {
//            size_t icomps = 0;
//            
//            for (const auto& tint : integrals)
//            {
//                icomps += tint.components<T1CPair, T1CPair>().size();
//            }
//            
//            label += std::to_string(icomps)  + ", ";
//        }
//        
//        label += "ket_width, ket_npgtos);";
//        
//        lines.push_back({4, 0, 1, label});
//    }
    
    lines.push_back({3, 0, 2, "}"});
}

bool
T2CECPFuncBodyDriver::_need_distances_ra(const I2CIntegral& integral) const
{
    return (integral[0] > integral[1]) && ((integral[0] + integral[1]) > 0);
}

bool
T2CECPFuncBodyDriver::_need_distances_rb(const I2CIntegral& integral) const
{
    return (integral[0] <= integral[1]) && ((integral[0] + integral[1]) > 0);
}
