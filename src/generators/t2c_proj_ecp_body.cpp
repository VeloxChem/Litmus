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

#include "t2c_proj_ecp_body.hpp"

#include "t2c_utils.hpp"

void
T2CProjECPFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                          const SM2Integrals&  integrals,
                                          const M2Integral&    integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_gtos_def())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_ket_variables_def(integrals))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_buffers_def(integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_loop_start(lines, integral);
    
    _add_ket_loop_start(lines, integral);
//    
//    _add_vrr_call_tree(lines, vrr_integrals, integral);
//    
//    _add_reduce_call_tree(lines, vrr_integrals, ctints, integral);
//    
    _add_ket_loop_end(lines, integrals, integral);
  
    _add_loop_end(lines, integral);
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CProjECPFuncBodyDriver::_get_gtos_def() const
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
T2CProjECPFuncBodyDriver::_get_ket_variables_def(const SM2Integrals& integrals) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned 2D arrays for ket side");
    
    vstr.push_back("CSimdArray<double> pfactors(9, ket_npgtos);");
    
    vstr.push_back("// allpcate I_n and L_n values");
    
    const auto imax = _get_max_bessel(integrals);
    
    const auto lmax = _get_max_momentum(integrals);

    vstr.push_back("CSimdArray<double> i_values(" +  std::to_string(imax + 1) + ", ket_npgtos);");
    
    vstr.push_back("CSimdArray<double> l_values(" + std::to_string(lmax + 1) + ", ket_npgtos);"); 

    return vstr;
}

int
T2CProjECPFuncBodyDriver::_get_max_momentum(const SM2Integrals& integrals) const
{
    int order = 0;
    
    for (const auto& [pref, tint] : integrals)
    {
        if (const auto torder = tint.order(); torder > order)
        {
            order = torder;
        }
    }
    
    return order;
}

int
T2CProjECPFuncBodyDriver::_get_max_bessel(const SM2Integrals& integrals) const
{
    int order = 0;
    
    for (const auto& [pref, tint] : integrals)
    {
        if (const auto torder = tint.order() + pref[0] + pref[1] + pref[2]; torder > order)
        {
            order = torder;
        }
    }
    
    return order;
}

std::vector<std::string>
T2CProjECPFuncBodyDriver::_get_buffers_def(const SM2Integrals& integrals,
                                           const M2Integral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned primitive integrals");
    
    size_t icomps = 0;
    
    for (const auto& [pref, tint] : integrals)
    {
        icomps += tint.components<T1CPair, T1CPair>().size();
    }
    
    auto label = "CSimdArray<double> pbuffer(" + std::to_string(icomps) + ", ket_npgtos);";
    
    vstr.push_back(label);
    
    vstr.push_back("// allocate aligned contracted integrals");
    
    icomps = integral.second.components<T1CPair, T1CPair>().size();
    
    label = "CSimdArray<double> cbuffer(" + std::to_string(icomps) + ", 1);";
    
    vstr.push_back(label);
    
    label = "CSimdArray<double> sbuffer(" + std::to_string(icomps) + ", 1);";
    
    vstr.push_back(label);
        
    return vstr;
}

void
T2CProjECPFuncBodyDriver::_add_loop_start(      VCodeLines& lines,
                                          const M2Integral& integral) const
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
    
    lines.push_back({2, 0, 2, "i_values.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "l_values.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "sbuffer.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "cbuffer.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "pbuffer.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "// loop over contracted basis functions on bra side"});
    
    lines.push_back({2, 0, 1, "for (auto j = bra_indices.first; j < bra_indices.second; j++)"});
    
    lines.push_back({2, 0, 1, "{"});
        
    lines.push_back({3, 0, 2, "cbuffer.zero();"});
    
    lines.push_back({3, 0, 2, "sbuffer.zero();"});

    lines.push_back({3, 0, 2, "const auto r_a = bra_gto_coords[j];"});
}

void
T2CProjECPFuncBodyDriver::_add_loop_end(      VCodeLines& lines,
                                        const M2Integral& integral) const
{
    std::string label;

    label = "t2cfunc::transform<"  + std::to_string(integral.second[0]);
            
    label += ", " + std::to_string(integral.second[1]) + ">(sbuffer, cbuffer);";
            
    lines.push_back({3, 0, 2, label});
    
    label = "distributor.distribute(sbuffer, ";
    
    label += "bra_gto_indices, ket_gto_indices, ";
            
    label += std::to_string(integral.second[0]) + ", ";
            
    label += std::to_string(integral.second[1]) + ", ";
    
    label += "j, ket_range, bra_eq_ket);";
    
    lines.push_back({3, 0, 1, label});
   
    lines.push_back({2, 0, 1, "}"});
    
    lines.push_back({1, 0, 1, "}"});
}

void
T2CProjECPFuncBodyDriver::_add_ket_loop_start(      VCodeLines& lines,
                                              const M2Integral& integral) const
{
    lines.push_back({3, 0, 1, "for (size_t k = 0; k < bra_npgtos; k++)"});
    
    lines.push_back({3, 0, 1, "{"});
    
    lines.push_back({4, 0, 2, "const auto a_exp = bra_gto_exps[k * bra_ncgtos + j];"});

    lines.push_back({4, 0, 2, "const auto a_norm = bra_gto_norms[k * bra_ncgtos + j];"});
    
    lines.push_back({4, 0, 1, "for (size_t l = 0; l < ecp_nppt; l++)"});
    
    lines.push_back({4, 0, 1, "{"});
    
    lines.push_back({5, 0, 2, "const auto c_exp = ecp_exps[l];"});

    lines.push_back({5, 0, 2, "const auto c_norm = ecp_facts[l];"});
    
//    lines.push_back({5, 0, 2, "t2cfunc::comp_coordinates_r(pfactors, 5, 2, r_a, a_exp, c_exp);"});
//    
//    if (_need_distances_ra(integral))
//    {
//        lines.push_back({5, 0, 2, "t2cfunc::comp_distances_ra(pfactors, 8, 5, r_a);"});
//    }
//    else
//    {
//        lines.push_back({5, 0, 2, "t2cfunc::comp_distances_rb(pfactors, 8, 5, 2);"});
//    }
//    
//    if ((integral[0] + integral[1]) > 1)
//    {
//        lines.push_back({5, 0, 2, "t2cfunc::comp_inverted_zeta(pfactors, 11, a_exp, c_exp);"});
//    }
}

void
T2CProjECPFuncBodyDriver::_add_ket_loop_end(      VCodeLines&   lines,
                                            const SM2Integrals& integrals,
                                            const M2Integral&   integral) const
{
    lines.push_back({4, 0, 1, "}"});
    
    lines.push_back({3, 0, 2, "}"});
}
