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
    
    _add_ket_loop_start(lines, integrals, integral);
    
    _add_aux_call_tree(lines, integrals, integral);

    _add_vrr_call_tree(lines, integrals, integral);
    
    _add_reduce_call_tree(lines, integrals, integral);
    
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
    
    const auto angpair = std::array<int, 2>({integral.second[0], integral.second[1]});
    
    icomps = t2c::number_of_spherical_components(angpair);
    
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
T2CProjECPFuncBodyDriver::_add_ket_loop_start(      VCodeLines&   lines,
                                              const SM2Integrals& integrals,
                                              const M2Integral&   integral) const
{
    lines.push_back({3, 0, 1, "for (size_t k = 0; k < bra_npgtos; k++)"});
    
    lines.push_back({3, 0, 1, "{"});
    
    lines.push_back({4, 0, 2, "const auto a_exp = bra_gto_exps[k * bra_ncgtos + j];"});

    lines.push_back({4, 0, 2, "const auto a_norm = bra_gto_norms[k * bra_ncgtos + j];"});
    
    lines.push_back({4, 0, 1, "for (size_t l = 0; l < ecp_nppt; l++)"});
    
    lines.push_back({4, 0, 1, "{"});
    
    lines.push_back({5, 0, 2, "const auto c_exp = ecp_exps[l];"});

    lines.push_back({5, 0, 2, "const auto c_norm = ecp_facts[l];"});
    
    lines.push_back({5, 0, 2, "t2cfunc::comp_coordinates_norm(pfactors, 5, 2);"});
                        
    lines.push_back({5, 0, 2, "t2cfunc::comp_legendre_args(pfactors, 6, 2, 5, r_a);"});
                        
    lines.push_back({5, 0, 2, "t2cfunc::comp_gamma_factors(pfactors, 7, 5, r_a, a_exp, c_exp);"});
                        
    lines.push_back({5, 0, 2, "t2cfunc::comp_bessel_args(pfactors, 8, 5, r_a, a_exp, c_exp);"});
    
    const auto imax = _get_max_bessel(integrals);
    
    const auto lmax = _get_max_momentum(integrals);
    
    lines.push_back({5, 0, 2, "t2cfunc::comp_i_vals(i_values, " + std::to_string(imax) + ", pfactors, 8);"});
                       
    lines.push_back({5, 0, 2, "t2cfunc::comp_l_vals(l_values, " + std::to_string(lmax) + ", pfactors, 8, 6);"});
}

void
T2CProjECPFuncBodyDriver::_add_ket_loop_end(      VCodeLines&   lines,
                                            const SM2Integrals& integrals,
                                            const M2Integral&   integral) const
{
    lines.push_back({4, 0, 1, "}"});
    
    lines.push_back({3, 0, 2, "}"});
}

void
T2CProjECPFuncBodyDriver::_add_reduce_call_tree(      VCodeLines&   lines,
                                                const SM2Integrals& integrals,
                                                const M2Integral&   integral) const
{
    const int spacer = 5;
    
    std::string label = "t2cfunc::reduce(cbuffer, 0, pbuffer, ";
        
    label += std::to_string(_get_position(integral, integrals)) + ", ";
        
    label += std::to_string(integral.second.components<T1CPair, T1CPair>().size()) + ", ";
        
    label += "ket_width, ket_npgtos);";
        
    lines.push_back({spacer, 0, 1, label});
}

size_t
T2CProjECPFuncBodyDriver::_get_position(const M2Integral&   integral,
                                        const SM2Integrals& integrals) const
{
    size_t pos = 0;
    
    for (const auto& tint : integrals)
    {
        if (tint == integral) return pos;
        
        pos += tint.second.components<T1CPair, T1CPair>().size();
    }
    
    return 0;
}

void
T2CProjECPFuncBodyDriver::_add_aux_call_tree(      VCodeLines&   lines,
                                             const SM2Integrals& integrals,
                                             const M2Integral&   integral) const
{
    const int spacer = 5;
    
    for (const auto& [pref, tint] : integrals)
    {
        if (!tint.is_simple()) continue;
        
        if ((tint[0] + tint[1]) == 0)
        {
            std::string label = "t2pecp::comp_prim_projected_core_potential_ss(";
            
            label += std::to_string(tint.order()) + ", " + std::to_string(pref[0]) + ", ";
            
            label += std::to_string(pref[1]) + ", " + std::to_string(pref[2]) + ", ";
            
            label += "pbuffer, " + std::to_string(_get_position(M2Integral(pref, tint), integrals)) + ", ";
            
            label += "i_values, l_values,  pfactors, 7, 5, r_a, a_norm, c_norm);";
            
            lines.push_back({spacer, 0, 2, label});
        }
    }
}

void
T2CProjECPFuncBodyDriver::_add_vrr_call_tree(      VCodeLines&   lines,
                                             const SM2Integrals& integrals,
                                             const M2Integral&   integral) const
{
    const int spacer = 5;
    
    SI2CIntegrals rints;
    
    // select non-auxilary integrals
    
    for (const auto& [pref, tint] : integrals)
    {
        if ((tint[0] + tint[1]) > 0)
        {
            rints.insert(tint);
        }
    }
    
    // write ordered VRR call tree
    
    for (const auto& rint : rints)
    {
        for (const auto& tint : integrals)
        {
            if (rint == tint.second)
            {
                auto label = "t2pecp::" + t2c::prim_compute_func_name(tint) + "(pbuffer, ";
                
                label += _get_vrr_arguments(tint, integrals); 
                
                label += "a_exp, c_exp);";
                
                lines.push_back({spacer, 0, 2, label});
            }
        }
    }
}

std::string
T2CProjECPFuncBodyDriver::_get_vrr_arguments(const M2Integral&   integral,
                                             const SM2Integrals& integrals) const
{
    auto label = std::to_string(_get_position(integral, integrals)) + ", ";
    
    for (const auto& tint : t2c::get_common_integrals(integral))
    {
        label += std::to_string(_get_position(tint, integrals)) + ", ";
    }
    
    if (integral.second[0] > 0)
    {
        label += std::to_string(integral.first[1]) + ", ";
        
    }
    else
    {
        label += std::to_string(integral.first[0]) + ", ";
    }
    
    const auto mrefint = (integral.second[0] > 0) ? M2Integral({0,1,0}, integral.second) : M2Integral({1,0,0}, integral.second);
    
    const auto mints = t2c::get_special_integrals(mrefint);
    
    const auto rints = t2c::get_special_integrals(integral);
    
    if (mints.size() == rints.size())
    {
        for (const auto& tint : rints)
        {
            label += std::to_string(_get_position(tint, integrals)) + ", ";
        }
    }
    else
    {
        for (const auto& tint : mints)
        {
            label += "-1, ";
        }
    }
    
    label += "pfactors, ";
    
    if (integral.second[0] > 0)
    {
        if (integral.second.order() > 0)
        {
            label += "2, ";
        }
        
        label += "r_a, ";
    }
    else
    {
        label += "2, ";
        
        if (integral.second.order() > 0)
        {
            label += "r_a, ";
        }
    }
    
    return label;
}
