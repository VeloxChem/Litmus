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

#include "t3c_geom_body.hpp"

#include "t2c_utils.hpp"
#include "t3c_utils.hpp"

void
T3CGeomFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                       const SG3Terms&      cterms,
                                       const SG3Terms&      skterms,
                                       const SI3CIntegrals& vrr_integrals,
                                       const I3CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_gto_pairs_def())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_ket_variables_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_prim_buffers_def(vrr_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_cart_buffers_def(cterms, integral))
    {
        lines.push_back({1, 0, 2, label});
    }

    for (const auto& label : _get_half_spher_buffers_def(skterms, integral))
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
    
    _add_loop_start(lines, integral);
    
    _add_ket_loop_start(lines, integral);
    
    _add_auxilary_integrals(lines, vrr_integrals, integral, 4);
    
    _add_vrr_call_tree(lines, vrr_integrals, integral, 4);
    
    _add_ket_loop_end(lines, cterms, vrr_integrals, integral);
    
    _add_bra_geom_call_tree(lines, cterms, integral);
    
    _add_bra_trafo_call_tree(lines, cterms, skterms, integral);
    
    _add_hrr_call_tree(lines, skterms, integral);
    
    _add_ket_trafo_call_tree(lines, skterms, integral);
    
    _add_loop_end(lines, integral);
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T3CGeomFuncBodyDriver::_get_gto_pairs_def() const
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
        
    vstr.push_back("const auto c_coords = ket_gto_pair_block.bra_coordinates();");
      
    vstr.push_back("const auto d_coords = ket_gto_pair_block.ket_coordinates();");
        
    vstr.push_back("const auto c_vec_exps = ket_gto_pair_block.bra_exponents();");
        
    vstr.push_back("const auto d_vec_exps = ket_gto_pair_block.ket_exponents();");
        
    vstr.push_back("const auto cd_vec_norms = ket_gto_pair_block.normalization_factors();");
        
    vstr.push_back("const auto cd_vec_ovls = ket_gto_pair_block.overlap_factors();");
        
    vstr.push_back("const auto c_indices = ket_gto_pair_block.bra_orbital_indices();");
        
    vstr.push_back("const auto d_indices = ket_gto_pair_block.ket_orbital_indices();");
        
    vstr.push_back("const auto ket_npgtos = ket_gto_pair_block.number_of_primitive_pairs();");
    
    return vstr;
}

std::vector<std::string>
T3CGeomFuncBodyDriver::_get_ket_variables_def(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned 2D arrays for ket side");
  
    // c_exps, d_exps, cd_ovls, cd_norms, c_coords, d_coords, q_coords, pq_coords, f_ss
    
    size_t nelems = 17;
    
    if (_need_center_w(integral)) nelems += 3;
    
    if (_need_distances_qd(integral)) nelems += 3;
    
    if (_need_distances_wq(integral)) nelems += 3;
    
    if (_need_distances_wa(integral)) nelems += 3;
        
    vstr.push_back("CSimdArray<double> pfactors(" + std::to_string(nelems) +  ", ket_npgtos);");
  
    if (_need_hrr(integral))
    {
        vstr.push_back("CSimdArray<double> cfactors(9, 1);");
    }
    
    return vstr;
}

bool
T3CGeomFuncBodyDriver::_need_center_w(const I3CIntegral& integral) const
{
    const auto orders = integral.prefixes_order();
    
    if (orders.empty())
    {
        return (integral[0] + integral[1] + integral[2]) > 0;
    }
    else
    {
        return (integral[0] + integral[1] + integral[2] + orders[0] + orders[1] + orders[2]) > 0;
    }
    
}

bool
T3CGeomFuncBodyDriver::_need_distances_qd(const I3CIntegral& integral) const
{
    const auto orders = integral.prefixes_order();
    
    if (orders.empty())
    {
        return (integral[1] + integral[2]) > 0;
    }
    else
    {
        return (integral[1] + integral[2] + orders[1] + orders[2]) > 0;
    }
}

bool
T3CGeomFuncBodyDriver::_need_distances_wq(const I3CIntegral& integral) const
{
    const auto orders = integral.prefixes_order();
    
    if (orders.empty())
    {
        return (integral[1] + integral[2]) > 0;
    }
    else
    {
        return (integral[1] + integral[2] + orders[1] + orders[2]) > 0;
    }
}

bool
T3CGeomFuncBodyDriver::_need_distances_wa(const I3CIntegral& integral) const
{
    const auto orders = integral.prefixes_order();
    
    if (orders.empty())
    {
        return integral[0]  > 0;
    }
    else
    {
        return (integral[0] + orders[0])  > 0;
    }
}

bool
T3CGeomFuncBodyDriver::_need_hrr(const I3CIntegral& integral) const
{
    const auto orders = integral.prefixes_order();
    
    if (orders.empty())
    {
        return integral[1] > 0;
    }
    else
    {
        return (integral[1] + orders[1]) > 0;
    }
}

size_t
T3CGeomFuncBodyDriver::_get_all_components(const SI3CIntegrals& integrals) const
{
    size_t tcomps= 0;
    
    for (const auto& tint : integrals)
    {
        tcomps += tint.components<T1CPair, T2CPair>().size();
    }
    
    return tcomps;
}

std::vector<std::string>
T3CGeomFuncBodyDriver::_get_prim_buffers_def(const SI3CIntegrals& integrals,
                                             const I3CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned primitive integrals");
    
    auto tcomps = _get_all_components(integrals);
    
    std::string label = "CSimdArray<double> pbuffer";
    
    label += "(" + std::to_string(tcomps) + ", ket_npgtos);";
    
    vstr.push_back(label);
    
    return vstr;
}

std::vector<std::string>
T3CGeomFuncBodyDriver::_get_cart_buffers_def(const SG3Terms&    cterms,
                                             const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned Cartesian integrals");
    
    size_t tcomps= 0;
    
    for (const auto& term : cterms)
    {
        tcomps += term.second.components<T1CPair, T2CPair>().size();
    }
    
    std::string label = "CSimdArray<double> cbuffer";
    
    label += "(" + std::to_string(tcomps) + ", 1);";
    
    vstr.push_back(label);
    
    return vstr;
}

std::vector<std::string>
T3CGeomFuncBodyDriver::_get_half_spher_buffers_def(const SG3Terms&    skterms,
                                                   const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (skterms.size() == 0) return vstr;
    
    vstr.push_back("// allocate aligned half transformed integrals");
    
    size_t tcomps = 0;
    
    for (const auto& term : skterms)
    {
        const auto tint = term.second;
                
        auto icomps = t2c::number_of_spherical_components(std::array<int, 1>({tint[0],}));
            
        auto angpair = std::array<int, 2>({tint[1], tint[2]});
                
        icomps *= t2c::number_of_cartesian_components(angpair);
        
        for (const auto& prefix : tint.prefixes())
        {
            icomps *= prefix.components().size();
        }
        
        tcomps += icomps;
    }
    
    std::string label = "CSimdArray<double> ";
            
    label += "skbuffer(" + std::to_string(tcomps) + ", 1);";
            
    vstr.push_back(label);
        
    return vstr;
}

std::vector<std::string>
T3CGeomFuncBodyDriver::_get_spher_buffers_def(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned spherical integrals");
                    
    auto tcomps = t2c::number_of_spherical_components(std::array<int, 1>({integral[0], }));
                    
    auto angpair = std::array<int, 2>({integral[1], integral[2]});
                    
    tcomps *= t2c::number_of_spherical_components(angpair);
    
    for (const auto& prefix : integral.prefixes())
    {
        tcomps *= prefix.components().size();
    }
    
    std::string label = "CSimdArray<double> ";
                    
    label += "sbuffer(" + std::to_string(tcomps) + ", 1);";
                    
    vstr.push_back(label);
   
    return vstr;
}

std::vector<std::string>
T3CGeomFuncBodyDriver::_get_boys_function_def(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto order = integral[0] + integral[1] + integral[2];
    
    for (const auto& porder : integral.prefixes_order())
    {
        order += porder;
    }
        
    vstr.push_back("// setup Boys fuction data");
        
    vstr.push_back("const CBoysFunc<" + std::to_string(order) + "> bf_table;");

    vstr.push_back("CSimdArray<double> bf_data(" + std::to_string(order + 2) + ", ket_npgtos);");
       
    return vstr;
}

void
T3CGeomFuncBodyDriver::_add_loop_start(      VCodeLines&    lines,
                                       const I3CIntegral&   integral) const
{
    lines.push_back({1, 0, 2, "// set up ket partitioning"});

    lines.push_back({1, 0, 2, "const auto ket_dim = ket_gto_pair_block.number_of_contracted_pairs();"});

    lines.push_back({1, 0, 2, "const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());"});

    lines.push_back({1, 0, 1, "for (size_t i = 0; i < ket_blocks; i++)"});
                    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), size_t{0});"});

    lines.push_back({2, 0, 2, "pfactors.load(c_vec_exps, ket_range, 0, ket_npgtos);"});
    
    lines.push_back({2, 0, 2, "pfactors.load(d_vec_exps, ket_range, 1, ket_npgtos);"});
 
    lines.push_back({2, 0, 2, "pfactors.load(cd_vec_ovls, ket_range, 2, ket_npgtos);"});
    
    lines.push_back({2, 0, 2, "pfactors.load(cd_vec_norms, ket_range, 3, ket_npgtos);"});
    
    lines.push_back({2, 0, 2, "pfactors.replicate_points(c_coords, ket_range, 4, ket_npgtos);"});
    
    lines.push_back({2, 0, 2, "pfactors.replicate_points(d_coords, ket_range, 7, ket_npgtos);"});
    
    if (_need_hrr(integral))
    {
        lines.push_back({2, 0, 2, "cfactors.replicate_points(c_coords, ket_range, 0, 1);"});
        
        lines.push_back({2, 0, 2, "cfactors.replicate_points(d_coords, ket_range, 3, 1);"});
        
        lines.push_back({2, 0, 2, "t4cfunc::comp_distances_cd(cfactors, 6, 0, 3);"});
    }
   
    lines.push_back({2, 0, 2, "// set up active SIMD width"});
    
    lines.push_back({2, 0, 2, "const auto ket_width = ket_range.second - ket_range.first;"});
    
    lines.push_back({2, 0, 2, "pbuffer.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "cbuffer.set_active_width(ket_width);"});
    
    if (_need_hrr(integral) || (integral[0] > 0))
    {
        lines.push_back({2, 0, 2, "skbuffer.set_active_width(ket_width);"});
    }
    
    lines.push_back({2, 0, 2, "sbuffer.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "bf_data.set_active_width(ket_width);"});
      
    lines.push_back({2, 0, 2, "// loop over basis function pairs on bra side"});

    lines.push_back({2, 0, 1, "for (auto j = bra_range.first; j < bra_range.second; j++)"});

    lines.push_back({2, 0, 1, "{"});
    
    lines.push_back({3, 0, 2, "// zero integral buffers"});
    
    lines.push_back({3, 0, 2, "cbuffer.zero();"});
    
    if (_need_hrr(integral) || (integral[0] > 0))
    {
        lines.push_back({3, 0, 2, "skbuffer.zero();"});
    }
    
    lines.push_back({3, 0, 2, "sbuffer.zero();"});

    lines.push_back({3, 0, 2, "// set up coordinates on bra side"});

    lines.push_back({3, 0, 2, "const auto r_a = bra_gto_coords[j];"});
}

void
T3CGeomFuncBodyDriver::_add_loop_end(      VCodeLines&  lines,
                                     const I3CIntegral& integral) const
{
    lines.push_back({2, 0, 1, "}"});
   
    lines.push_back({1, 0, 2, "}"});
}

void
T3CGeomFuncBodyDriver::_add_ket_loop_start(      VCodeLines&  lines,
                                           const I3CIntegral& integral) const
{
    lines.push_back({3, 0, 1, "for (int k = 0; k < bra_npgtos; k++)"});
   
    lines.push_back({3, 0, 1, "{"});
    
    lines.push_back({4, 0, 2, "const auto a_exp = bra_gto_exps[k * bra_ncgtos + j];"});
            
    lines.push_back({4, 0, 2, "const auto a_norm = bra_gto_norms[k * bra_ncgtos + j];"});
  
    lines.push_back({4, 0, 2, "t4cfunc::comp_coordinates_q(pfactors, 10, 4, 7);"});
    
    lines.push_back({4, 0, 2, "t3cfunc::comp_distances_aq(pfactors, 13, 10, r_a);"});
    
    if (_need_center_w(integral))
    {
        auto label_w = std::to_string(_get_index_w(integral));
        
        lines.push_back({4, 0, 2, "t3cfunc::comp_coordinates_w(pfactors, " + label_w + ", 10, r_a, a_exp);"});
    }
    
    if (_need_distances_qd(integral))
    {
        auto label_qd = std::to_string(_get_index_qd(integral));
        
        lines.push_back({4, 0, 2, "t4cfunc::comp_distances_qd(pfactors, " + label_qd + ", 10, 7);"});
    }
    
    if (_need_distances_wq(integral))
    {
        auto label_w = std::to_string(_get_index_w(integral));
        
        auto label_wq = std::to_string(_get_index_wq(integral));
        
        lines.push_back({4, 0, 2, "t4cfunc::comp_distances_wq(pfactors, " + label_wq + ", " + label_w + ", 10);"});
    }
     
    if (_need_distances_wa(integral))
    {
        auto label_w = std::to_string(_get_index_w(integral));
        
        auto label_wa = std::to_string(_get_index_wa(integral));
        
        lines.push_back({4, 0, 2, "t4cfunc::comp_distances_wp(pfactors, " + label_wa + ", " + label_w + ", r_a);"});
    }
    
    const auto gorders = integral.prefixes_order();
    
    auto border = integral[0] + integral[1] + integral[2] + 1;
    
    border += gorders[0] + gorders[1] + gorders[2];
    
    lines.push_back({4, 0, 2, "t3cfunc::comp_boys_args(bf_data, " + std::to_string(border) + ", pfactors, 13, a_exp);"});
    
    lines.push_back({4, 0, 2, "bf_table.compute(bf_data, 0, " + std::to_string(border) + ");"});
    
    lines.push_back({4, 0, 2, "t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);"});
}

void
T3CGeomFuncBodyDriver::_add_ket_loop_end(      VCodeLines&    lines,
                                         const SG3Terms&      cterms,
                                         const SI3CIntegrals& vrr_integrals,
                                         const I3CIntegral&   integral) const
{
    // non-scaled integrals
    
    for (const auto& term : cterms)
    {
        if (!term.second.prefixes().empty()) continue;
        
        if (term.first == std::array<int, 3>({0, 0, 0}))
        {
            const auto tint = term.second;
            
            std::string label = "t2cfunc::reduce(cbuffer, ";
                
            label +=  std::to_string(_get_index(term, cterms)) + ", ";
                
            label += "pbuffer, ";
                
            label += std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
                
            label += std::to_string(tint.components<T1CPair, T2CPair>().size()) + ", ";
                
            label += "ket_width, ket_npgtos);";
                
            lines.push_back({4, 0, 2, label});
        }
    }
    
    // scales integrals on center A

    for (const auto& term : cterms)
    {
        if (!term.second.prefixes().empty()) continue;
        
        if (term.first == std::array<int, 3>({1, 0, 0}))
        {
            const auto tint = term.second;
            
            std::string label = "pbuffer.scale(2.0 * a_exp, {";
           
            label +=  std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
           
            label +=  std::to_string(_get_index(0, tint, vrr_integrals) + tint.components<T1CPair, T2CPair>().size()) + "});";
           
            lines.push_back({4, 0, 2, label});
        }
    }
    
    for (const auto& term : cterms)
    {
        if (!term.second.prefixes().empty()) continue;
        
        if (term.first == std::array<int, 3>({1, 0, 0}))
        {
            const auto tint = term.second;
            
            std::string label = "t2cfunc::reduce(cbuffer, ";
                
            label +=  std::to_string(_get_index(term, cterms)) + ", ";
                
            label += "pbuffer, ";
                
            label += std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
                
            label += std::to_string(tint.components<T1CPair, T2CPair>().size()) + ", ";
                
            label += "ket_width, ket_npgtos);";
                
            lines.push_back({4, 0, 2, label});
        }
    }
                
    lines.push_back({3, 0, 2, "}"});
}

size_t
T3CGeomFuncBodyDriver::_get_index_w(const I3CIntegral& integral) const
{
    return 17;
}

size_t
T3CGeomFuncBodyDriver::_get_index_qd(const I3CIntegral& integral) const
{
    auto index = _get_index_w(integral);
    
    if (_need_center_w(integral)) index += 3;
    
    return index;
}

size_t
T3CGeomFuncBodyDriver::_get_index_wq(const I3CIntegral& integral) const
{
    auto index = _get_index_qd(integral);
    
    if (_need_distances_qd(integral)) index += 3;
    
    return index;
}

size_t
T3CGeomFuncBodyDriver::_get_index_wa(const I3CIntegral& integral) const
{
    auto index = _get_index_wq(integral);
    
    if (_need_distances_wq(integral)) index += 3;
    
    return index;
}

void
T3CGeomFuncBodyDriver::_add_auxilary_integrals(      VCodeLines&    lines,
                                               const SI3CIntegrals& integrals,
                                               const I3CIntegral&   integral,
                                               const size_t         spacer) const
{
    for (const auto& tint : integrals)
    {
        if ((tint[0] + tint[1] + tint[2]) == 0)
        {
            const auto blabel = std::to_string(tint.order());
            
            const auto ilabel = std::to_string(_get_index(0, tint, integrals));
                    
            lines.push_back({spacer, 0, 2, "t3ceri::comp_prim_electron_repulsion_sss(pbuffer, " + ilabel + ", pfactors, 16, bf_data, " + blabel + ");"});
        }
    }
}

size_t
T3CGeomFuncBodyDriver::_get_index(const size_t         start,
                                  const I3CIntegral&   integral,
                                  const SI3CIntegrals& integrals) const
{
    size_t index = start;
    
    for (const auto& tint : integrals)
    {
        if (tint == integral) return index;
        
        index += tint.components<T1CPair, T2CPair>().size();
    }
    
    return 0;
}

void
T3CGeomFuncBodyDriver::_add_vrr_call_tree(      VCodeLines&  lines,
                                          const SI3CIntegrals& integrals,
                                          const I3CIntegral&   integral,
                                          const size_t         spacer) const
{
    for (const auto& tint : integrals)
    {
        if ((tint[1] == 0) && ((tint[0] + tint[2]) > 0))
        {
            const auto name = t3c::prim_compute_func_name(tint);
            
            auto label = t3c::namespace_label(tint) + "::" + name + "(pbuffer, ";
            
            label += _get_vrr_arguments(0, integrals, tint);
            
            label += "pfactors, ";
            
            if (_need_distances_wa(tint))
            {
                label += std::to_string(_get_index_wa(integral)) + ", ";
            }
            else
            {
                label += std::to_string(_get_index_qd(integral)) + ", ";
                
                label += std::to_string(_get_index_wq(integral)) + ", ";
            }
            
            if ((tint[0] + tint[2]) > 1)
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

std::string
T3CGeomFuncBodyDriver::_get_vrr_arguments(const size_t start,
                                          const SI3CIntegrals& integrals,
                                          const I3CIntegral&  integral) const
{
    std::string label = std::to_string(_get_index(start, integral, integrals)) + ", ";
    
    for (const auto& tint : t3c::get_vrr_integrals(integral))
    {
        label += std::to_string(_get_index(start, tint, integrals)) + ", ";
    }
    
    return label;
}

size_t
T3CGeomFuncBodyDriver::_get_index(const G3Term&   term,
                                  const SG3Terms& terms) const
{
    size_t index = 0;
    
    for (const auto& cterm : terms)
    {
        if (term == cterm) return index;
        
        index += cterm.second.components<T1CPair, T2CPair>().size();
    }
    
    return 0;
}

void
T3CGeomFuncBodyDriver::_add_ket_trafo_call_tree(      VCodeLines&  lines,
                                                const SG3Terms&    skterms,
                                                const I3CIntegral& integral) const
{
    size_t gcomps = 1;
    
    for (const auto& prefix : integral.prefixes())
    {
        gcomps *= prefix.components().size();
    }
    
    auto angpair = std::array<int, 2>({integral[1], integral[2]});
    
    auto kccomps = t2c::number_of_cartesian_components(angpair);
    
    auto kscomps = t2c::number_of_spherical_components(angpair);

    auto bscomps = t2c::number_of_spherical_components(std::array<int, 1>({integral[0], }));
    
    auto gterm = t3c::prune_term(G3Term({std::array<int, 3>({0, 0, 0}), integral}));
    
    const auto gindex = _get_half_spher_index(gterm, skterms);
    
    for (size_t i = 0; i < gcomps; i++)
    {
        std::string label = "t3cfunc::ket_transform<" + std::to_string(integral[1]) + ", " + std::to_string(integral[2]) + ">";
            
        if (_need_hrr(integral) || (integral[0] > 0))
        {
            label += "(sbuffer, " + std::to_string(i * bscomps * kscomps) + ", skbuffer, ";
        }
        else
        {
            label += "(sbuffer, " + std::to_string(i * bscomps * kscomps) + ", cbuffer, ";
        }
        
        label += std::to_string(gindex + i * kccomps * bscomps) + ", ";
        
        label += std::to_string(integral[0]) + ");";
            
        lines.push_back({3, 0, 2, label});
    }
    
    std::string label = "distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, ";
    
    label += std::to_string(integral[0]) + ", ";
    
    label += std::to_string(integral[1]) + ", ";
    
    label += std::to_string(integral[2]) + ", ";
    
    label += "j, ket_range);";
    
    lines.push_back({3, 0, 1, label});
}

size_t
T3CGeomFuncBodyDriver::_get_half_spher_index(const G3Term&   term,
                                             const SG3Terms& terms) const
{
    size_t index = 0;
    
    for (const auto& cterm : terms)
    {
        if (term == cterm) return index;
        
        const auto tint = cterm.second;
        
        auto icomps = t2c::number_of_spherical_components(std::array<int, 1>({tint[0], }));
            
        auto angpair = std::array<int, 2>({tint[1], tint[2]});
                
        icomps *= t2c::number_of_cartesian_components(angpair);
        
        for (const auto& prefix : tint.prefixes())
        {
            icomps *= prefix.components().size();
        }
        
        index += icomps;
    }
    
    return 0;
}

void
T3CGeomFuncBodyDriver::_add_bra_trafo_call_tree(      VCodeLines&  lines,
                                                const SG3Terms&    cterms,
                                                const SG3Terms&    skterms,
                                                const I3CIntegral& integral) const
{
    if ((integral[0] + integral[1]) > 0)
    {
        for (const auto& skterm : skterms)
        {
            const auto tint = skterm.second;
            
            const auto angpair = std::array<int, 2>({tint[1], tint[2]});
            
            const auto ket_comps = t2c::number_of_cartesian_components(angpair);
            
            const auto bra_ccomps = t2c::number_of_cartesian_components(std::array<int, 1>{tint[0],});
            
            const auto bra_scomps = t2c::number_of_spherical_components(std::array<int, 1>{tint[0],});
            
            if ((integral[0] == 0) && (tint[0] == 1) && tint[1] == 0)
            {
                for (int i = 0; i < 3; i++)
                {
                    std::string label = "t3cfunc::bra_transform<" + std::to_string(integral[0]) + ">";
                    
                    label += "(skbuffer, " + std::to_string(_get_half_spher_index(skterm, skterms) + i * ket_comps) + ", ";
                    
                    label += "cbuffer, " + std::to_string(_get_index(skterm, cterms) + i * ket_comps)  + ", ";
                    
                    label += std::to_string(tint[1]) + ", " + std::to_string(tint[2]) + ");";
                    
                    lines.push_back({3, 0, 2, label});
                }
            }
            else
            {
                if (tint[1] == 0)
                {
                    for (int i = 0; i < 3; i++)
                    {
                        std::string label = "t3cfunc::bra_transform<" + std::to_string(integral[0]) + ">";
                        
                        label += "(skbuffer, " + std::to_string(_get_half_spher_index(skterm, skterms) + i * bra_scomps * ket_comps) + ", ";
                        
                        label += "cbuffer, " + std::to_string(_get_index(skterm, cterms) + i * bra_ccomps * ket_comps)  + ", ";
                        
                        label += std::to_string(tint[1]) + ", " + std::to_string(tint[2]) + ");";
                        
                        lines.push_back({3, 0, 2, label});
                    }
                }
            }
        }
    }
}

void
T3CGeomFuncBodyDriver::_add_hrr_call_tree(      VCodeLines&  lines,
                                          const SG3Terms&    skterms,
                                          const I3CIntegral& integral) const
{
    for (const auto& term : skterms)
    {
        const auto tint = term.second;
        
        if (tint[1] > 0)
        {
            const auto gorders = tint.prefixes_order();
            
            if (gorders == std::vector<int>({0, 0, 0}))
            {
                const auto name = t3c::hrr_compute_func_name(tint);

                auto label = t3c::namespace_label(tint) + "::" + name + "(skbuffer, ";
                
                label += std::to_string(_get_half_spher_index(term, skterms)) + ", ";
                
                label += _get_hrr_arguments(skterms, term);
                            
                label += "cfactors, 6, ";
                
                label += std::to_string(tint[0]);
                
                label += ");";
                
                lines.push_back({3, 0, 2, label});
            }
            
            if (gorders == std::vector<int>({1, 0, 0}))
            {
                const auto bra_comps = t2c::number_of_spherical_components(std::array<int, 1>{tint[0],});
                
                const auto ket_comps = t2c::number_of_cartesian_components(std::array<int, 2>{tint[1], tint[2],});
                
                for (int i = 0; i < 3; i++)
                {
                    const auto name = t3c::hrr_compute_func_name(tint);

                    auto label = t3c::namespace_label(tint) + "::" + name + "(skbuffer, ";
                    
                    label += std::to_string(_get_half_spher_index(term, skterms) + i * bra_comps * ket_comps) + ", ";
                    
                    label += _get_hrr_arguments(skterms, term, i);
                                
                    label += "cfactors, 6, ";
                    
                    label += std::to_string(tint[0]);
                    
                    label += ");";
                    
                    lines.push_back({3, 0, 2, label});
                }
            }
        }
    }
}

std::string
T3CGeomFuncBodyDriver::_get_hrr_arguments(const SG3Terms& skterms,
                                          const G3Term&   term) const
{
    std::string label;
    
    for (const auto& tint : t3c::get_hrr_integrals(term.second))
    {
        label += std::to_string(_get_half_spher_index(G3Term({term.first, tint}), skterms)) + ", ";
    }
    
    return label;
}

std::string
T3CGeomFuncBodyDriver::_get_hrr_arguments(const SG3Terms& skterms,
                                          const G3Term&   term,
                                          const int       icomponent) const
{
    std::string label;
    
    for (const auto& tint : t3c::get_hrr_integrals(term.second.base()))
    {
        std::cout << "XXX"  <<  tint.prefix_label() << " " << tint.label() << std::endl;
        
        const auto bra_comps = t2c::number_of_spherical_components(std::array<int, 1>{tint[0], });
        
        const auto ket_comps = t2c::number_of_cartesian_components(std::array<int, 2>{tint[1], tint[2],});
        
        auto ctint = tint;
        
        ctint.set_prefixes(term.second.prefixes());
        
        label += std::to_string(_get_half_spher_index(G3Term({term.first, ctint}), skterms) + icomponent * bra_comps * ket_comps) + ", ";
    }
    
    return label;
}

void
T3CGeomFuncBodyDriver::_add_bra_geom_call_tree(      VCodeLines&  lines,
                                               const SG3Terms&    cterms,
                                               const I3CIntegral& integral) const
{
    
    for (const auto& term : cterms)
    {
        const auto tint = term.second;
        
        if (tint.prefixes_order() == std::vector<int>({1, 0, 0}))
        {
            const auto name = t3c::bra_geom_compute_func_name(tint);
           
            auto label = t3c::namespace_label(tint) + "::" + name + "(cbuffer, ";
            
            label += std::to_string(_get_index(term, cterms)) + ", ";
            
            label += _get_bra_geom_arguments(term, cterms);
            
            label += std::to_string(tint[1]) + ", " + std::to_string(tint[2]);
            
            label += ");";
            
            lines.push_back({3, 0, 2, label});
        }
    }
}

std::string
T3CGeomFuncBodyDriver::_get_bra_geom_arguments(const G3Term&   term,
                                               const SG3Terms& cterms) const
{
    std::string label;
    
    const auto tint = term.second;
    
    if (tint.prefixes_order() == std::vector<int>({1, 0, 0}))
    {
        for (const auto& rtint : t3c::get_bra_geom_integrals(tint))
        {
            const auto rterm = (tint[0] > rtint[0]) ? G3Term(std::array<int, 3>({0, 0, 0}), rtint)
                                                    : G3Term(std::array<int, 3>({1, 0, 0}), rtint);
            
            label += std::to_string(_get_index(rterm, cterms))  + ", ";
        }
    }
    
    return label;
}

