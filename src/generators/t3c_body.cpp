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

#include "t3c_body.hpp"

#include "t2c_utils.hpp"
#include "t3c_utils.hpp"

void
T3CFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                   const SI3CIntegrals& hrr_integrals,
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
    
    for (const auto& label : _get_cart_buffers_def(hrr_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_half_spher_buffers_def(hrr_integrals, integral))
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
    
    _add_loop_start(lines, hrr_integrals, integral);

    _add_ket_loop_start(lines, integral);

    _add_auxilary_integrals(lines, vrr_integrals, integral, 4);

    _add_vrr_call_tree(lines, vrr_integrals, integral, 4);

    _add_ket_loop_end(lines, vrr_integrals, hrr_integrals, integral);
    
    _add_bra_trafo_call_tree(lines, hrr_integrals, integral);

    _add_hrr_call_tree(lines, hrr_integrals, integral);
    
    _add_ket_trafo_call_tree(lines, hrr_integrals, integral);

    _add_loop_end(lines, integral);
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T3CFuncBodyDriver::_get_gto_pairs_def() const
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
T3CFuncBodyDriver::_get_ket_variables_def(const I3CIntegral& integral) const
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

std::vector<std::string>
T3CFuncBodyDriver::_get_cart_buffers_def(const SI3CIntegrals& integrals,
                                         const I3CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned Cartesian integrals");
    
    auto tcomps = _get_all_components(_get_cart_buffer_integrals(integrals));
    
    std::string label = "CSimdArray<double> cbuffer";
    
    label += "(" + std::to_string(tcomps) + ", 1);";
    
    vstr.push_back(label);
    
    return vstr;
}

bool
T3CFuncBodyDriver::_need_center_w(const I3CIntegral& integral) const
{
    return (integral[0] + integral[1] + integral[2]) > 0;
}

bool
T3CFuncBodyDriver::_need_distances_qd(const I3CIntegral& integral) const
{
    return (integral[1] + integral[2]) > 0;
}

bool
T3CFuncBodyDriver::_need_distances_wq(const I3CIntegral& integral) const
{
    return (integral[1] + integral[2]) > 0;
}

bool
T3CFuncBodyDriver::_need_distances_wa(const I3CIntegral& integral) const
{
    return integral[0]  > 0;
}

bool
T3CFuncBodyDriver::_need_hrr(const I3CIntegral& integral) const
{
    return integral[1] > 0;
}
std::vector<std::string>
T3CFuncBodyDriver::_get_prim_buffers_def(const SI3CIntegrals& integrals,
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

size_t
T3CFuncBodyDriver::_get_all_components(const SI3CIntegrals& integrals) const
{
    size_t tcomps= 0;
    
    for (const auto& tint : integrals)
    {
        tcomps += tint.components<T1CPair, T2CPair>().size();
    }
    
    return tcomps;
}

SI3CIntegrals
T3CFuncBodyDriver::_get_cart_buffer_integrals(const SI3CIntegrals& integrals) const
{
    SI3CIntegrals tints;
    
    for (const auto& tint : integrals)
    {
        if (tint[1]  == 0)
        {
            tints.insert(tint);
        }
    }
    
    return tints;
}

std::vector<std::string>
T3CFuncBodyDriver::_get_half_spher_buffers_def(const SI3CIntegrals& integrals,
                                               const I3CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    if (_need_hrr(integral) || (integral[0] > 0))
    {
        vstr.push_back("// allocate aligned half transformed integrals");
        
        auto tcomps = _get_all_half_spher_components(_get_half_spher_buffers_integrals(integrals, integral));
        
        std::string label = "CSimdArray<double> ";
                
        label += "skbuffer(" + std::to_string(tcomps) + ", 1);";
                
        vstr.push_back(label);
    }
    
    return vstr;
}

SI3CIntegrals
T3CFuncBodyDriver::_get_half_spher_buffers_integrals(const SI3CIntegrals& integrals,
                                                     const I3CIntegral&   integral) const
{
    SI3CIntegrals tints;
    
    for (const auto& tint : integrals)
    {
        if (tint[0] == integral[0])
        {
            tints.insert(tint);
        }
    }
    
    tints.insert(integral);
    
    return tints;
}

size_t
T3CFuncBodyDriver::_get_all_half_spher_components(const SI3CIntegrals& integrals) const
{
    size_t tcomps = 0;
    
    for (const auto& tint : integrals)
    {
        auto angpair = std::array<int, 2>({0, tint[0]});
                
        auto icomps = t2c::number_of_spherical_components(angpair);
            
        angpair = std::array<int, 2>({tint[1], tint[2]});
                
        icomps *= t2c::number_of_cartesian_components(angpair);
        
        tcomps += icomps;
    }
    
    return tcomps;
}

std::vector<std::string>
T3CFuncBodyDriver::_get_spher_buffers_def(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned spherical integrals");
                    
    const auto tcomps = _get_all_spher_components(integral);
    
    std::string label = "CSimdArray<double> ";
                    
    label += "sbuffer(" + std::to_string(tcomps) + ", 1);";
                    
    vstr.push_back(label);
   
    return vstr;
}

size_t
T3CFuncBodyDriver::_get_all_spher_components(const I3CIntegral& integral) const
{
    auto angpair = std::array<int, 2>({integral[1], integral[2]});
                    
    auto tcomps = t2c::number_of_spherical_components(angpair);
                    
    angpair = std::array<int, 2>({0, integral[0]});
                    
    tcomps *= t2c::number_of_spherical_components(angpair);
    
    return tcomps;
}

std::vector<std::string>
T3CFuncBodyDriver::_get_boys_function_def(const I3CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto order = integral[0] + integral[1] + integral[2];
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        for (const auto& prefix : prefixes)
        {
            order += prefix.shape().order();
        }
    }
        
    vstr.push_back("// setup Boys fuction data");
        
    vstr.push_back("const CBoysFunc<" + std::to_string(order) + "> bf_table;");

    vstr.push_back("CSimdArray<double> bf_data(" + std::to_string(order + 2) + ", ket_npgtos);");
       
    return vstr;
}

void
T3CFuncBodyDriver::_add_loop_start(      VCodeLines&    lines,
                                   const SI3CIntegrals& integrals,
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
        //lines.push_back({2, 0, 2, "ckbuffer.set_active_width(ket_width);"});
        
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
        // lines.push_back({3, 0, 2, "ckbuffer.zero();"});
        
        lines.push_back({3, 0, 2, "skbuffer.zero();"});
    }
    
    lines.push_back({3, 0, 2, "sbuffer.zero();"});

    lines.push_back({3, 0, 2, "// set up coordinates on bra side"});

    lines.push_back({3, 0, 2, "const auto r_a = bra_gto_coords[j];"});
    
    //lines.push_back({3, 0, 2, "const auto a_xyz = r_a.coordinates();"});
}

void
T3CFuncBodyDriver::_add_loop_end(      VCodeLines&  lines,
                                 const I3CIntegral& integral) const
{
    lines.push_back({2, 0, 1, "}"});
   
    lines.push_back({1, 0, 2, "}"});
}

void
T3CFuncBodyDriver::_add_ket_loop_start(      VCodeLines&  lines,
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
    
    auto border = integral[0] + integral[1] + integral[2] + 1;
    
    lines.push_back({4, 0, 2, "t3cfunc::comp_boys_args(bf_data, " + std::to_string(border) + ", pfactors, 13, a_exp);"});
    
    lines.push_back({4, 0, 2, "bf_table.compute(bf_data, 0, " + std::to_string(border) + ");"});
    
    lines.push_back({4, 0, 2, "t3cfunc::comp_ovl_factors(pfactors, 16, 2, 3, a_norm, a_exp);"});
}

size_t
T3CFuncBodyDriver::_get_index_w(const I3CIntegral& integral) const
{
    return 17;
}

size_t
T3CFuncBodyDriver::_get_index_qd(const I3CIntegral& integral) const
{
    auto index = _get_index_w(integral);
    
    if (_need_center_w(integral)) index += 3;
    
    return index;
}

size_t
T3CFuncBodyDriver::_get_index_wq(const I3CIntegral& integral) const
{
    auto index = _get_index_qd(integral);
    
    if (_need_distances_qd(integral)) index += 3;
    
    return index;
}

size_t
T3CFuncBodyDriver::_get_index_wa(const I3CIntegral& integral) const
{
    auto index = _get_index_wq(integral);
    
    if (_need_distances_wq(integral)) index += 3;
    
    return index;
}

void
T3CFuncBodyDriver::_add_auxilary_integrals(      VCodeLines&    lines,
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
T3CFuncBodyDriver::_get_index(const size_t         start,
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
T3CFuncBodyDriver::_add_vrr_call_tree(      VCodeLines&  lines,
                                      const SI3CIntegrals& integrals,
                                      const I3CIntegral&   integral,
                                      const size_t         spacer) const
{
    size_t nterms = 0;
    
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
            
            nterms++;
        }
    }
}

std::string
T3CFuncBodyDriver::_get_vrr_arguments(const size_t start,
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

void
T3CFuncBodyDriver::_add_ket_loop_end(      VCodeLines&  lines,
                                     const SI3CIntegrals& vrr_integrals,
                                     const SI3CIntegrals& hrr_integrals,
                                     const I3CIntegral& integral) const
{
    const auto cints = _get_cart_buffer_integrals(hrr_integrals);
    
    for (const auto& tint : cints)
    {
        if (tint[1]  == 0)
        {
            std::string label = "t2cfunc::reduce(cbuffer, ";
            
            label +=  std::to_string(_get_index(0, tint, cints)) + ", ";
            
            label += "pbuffer, ";
            
            label += std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
            
            label += std::to_string(tint.components<T1CPair, T2CPair>().size()) + ", ";
            
            label += "ket_width, ket_npgtos);";
            
            lines.push_back({4, 0, 2, label});
        }
    }
    
    lines.push_back({3, 0, 2, "}"});
}

void
T3CFuncBodyDriver::_add_bra_trafo_call_tree(      VCodeLines&    lines,
                                            const SI3CIntegrals& integrals,
                                            const I3CIntegral&   integral) const
{
    const auto skints = _get_half_spher_buffers_integrals(integrals, integral);
    
    const auto ckints = _get_cart_buffer_integrals(integrals);
    
    if ((integral[0] + integral[1]) > 0) 
    {
        for (const auto& tint : ckints)
        {
            if (tint[0] == integral[0])
            {
                std::string label = "t3cfunc::bra_transform<" + std::to_string(tint[0]) + ">";
                
                label += "(skbuffer, " + std::to_string(_get_half_spher_index(0, tint, skints)) + ", ";
                
                label += "cbuffer, " + std::to_string(_get_index(0, tint, ckints))  + ", ";
                
                label += std::to_string(tint[1]) + ", " + std::to_string(tint[2]) + ");";
                
                lines.push_back({3, 0, 2, label});
            }
        }
    }
}

void
T3CFuncBodyDriver::_add_hrr_call_tree(      VCodeLines&  lines,
                                      const SI3CIntegrals& integrals,
                                      const I3CIntegral&   integral) const
{
    const auto skints = _get_half_spher_buffers_integrals(integrals, integral);
    
    for (const auto& tint : skints)
    {
        if (tint[1] > 0)
        {
            const auto name = t3c::hrr_compute_func_name(tint);

            auto label = t3c::namespace_label(tint) + "::" + name + "(skbuffer, ";
            
            label += std::to_string(_get_half_spher_index(0, tint, skints)) + ", ";
            
            label += _get_hrr_arguments(0, tint, integrals);
                        
            label += "cfactors, 6, ";
            
            label += std::to_string(tint[0]);
            
            label += ");";
            
            lines.push_back({3, 0, 2, label});
        }
    }
}

size_t
T3CFuncBodyDriver::_get_half_spher_index(const size_t         start,
                                         const I3CIntegral&   integral,
                                         const SI3CIntegrals& integrals) const
{
    size_t index = start;
    
    for (const auto& tint : integrals)
    {
        if (tint == integral) return index;
        
        auto angpair = std::array<int, 2>({0, tint[0]});
                
        auto icomps = t2c::number_of_spherical_components(angpair);
            
        angpair = std::array<int, 2>({tint[1], tint[2]});
                
        icomps *= t2c::number_of_cartesian_components(angpair);
        
        index += icomps;
    }
    
    return index;
}

std::string
T3CFuncBodyDriver::_get_hrr_arguments(const size_t       start,
                                      const I3CIntegral& integral,
                                      const SI3CIntegrals& integrals) const
{
    std::string label;
    
    const auto skints = _get_half_spher_buffers_integrals(integrals, integral);
    
    for (const auto& tint : t3c::get_hrr_integrals(integral))
    {
        label += std::to_string(_get_half_spher_index(start, tint, skints)) + ", ";
    }
    
    return label;
}

void
T3CFuncBodyDriver::_add_ket_trafo_call_tree(      VCodeLines&    lines,
                                            const SI3CIntegrals& integrals,
                                            const I3CIntegral&   integral) const
{
    const auto skints = _get_half_spher_buffers_integrals(integrals, integral);
    
    std::string label = "t3cfunc::ket_transform<" + std::to_string(integral[1]) + ", " + std::to_string(integral[2]) + ">";
        
    if (_need_hrr(integral) || (integral[0] > 0))
    {
        label += "(sbuffer, 0, skbuffer, ";
    }
    else
    {
        label += "(sbuffer, 0, cbuffer, ";
    }
    
    label += std::to_string(_get_half_spher_index(0, integral, skints)) + ", ";
    
    label += std::to_string(integral[0]) + ");";
        
    lines.push_back({3, 0, 2, label});
   
    label = "distributor.distribute(sbuffer, 0, bra_gto_indices, c_indices, d_indices, ";
    
    label += std::to_string(integral[0]) + ", ";
    
    label += std::to_string(integral[1]) + ", ";
    
    label += std::to_string(integral[2]) + ", ";
    
    label += "j, ket_range);";
    
    lines.push_back({3, 0, 1, label});
}
