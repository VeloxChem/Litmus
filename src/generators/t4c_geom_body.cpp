#include "t4c_geom_body.hpp"

#include <algorithm>

#include "t4c_utils.hpp"
#include "t2c_utils.hpp"
#include "t4c_vrr_eri_driver.hpp"
#include "t4c_center_driver.hpp"

void
T4CGeomFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                       const SI4CIntegrals& geom_integrals,
                                       const SI4CIntegrals& bra_base_integrals,
                                       const SI4CIntegrals& bra_rec_base_integrals,
                                       const SI4CIntegrals& ket_base_integrals,
                                       const SI4CIntegrals& ket_rec_base_integrals,
                                       const SI4CIntegrals& vrr_integrals,
                                       const I4CIntegral&   integral) const
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
    
    for (const auto& label : _get_cart_buffers_def(bra_base_integrals, bra_rec_base_integrals, ket_base_integrals, ket_rec_base_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_contr_buffers_def(ket_base_integrals, ket_rec_base_integrals))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_half_spher_buffers_def(geom_integrals, bra_base_integrals, bra_rec_base_integrals, ket_base_integrals, ket_rec_base_integrals, integral))
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
    
    _add_ket_loop_end(lines, bra_base_integrals, bra_rec_base_integrals, ket_base_integrals, ket_rec_base_integrals, vrr_integrals, integral);
    
    _add_ket_hrr_call_tree(lines, bra_base_integrals, bra_rec_base_integrals, ket_base_integrals, ket_rec_base_integrals, integral, 3);
    
    _add_ket_trafo_call_tree(lines, bra_base_integrals, bra_rec_base_integrals, ket_base_integrals, ket_rec_base_integrals, integral, 3);
    
    _add_bra_hrr_call_tree(lines, bra_base_integrals, bra_rec_base_integrals, ket_base_integrals, ket_rec_base_integrals, integral, 3);
    
    _add_bra_geom_hrr_call_tree(lines, geom_integrals, bra_base_integrals, bra_rec_base_integrals, ket_base_integrals, ket_rec_base_integrals, integral, 3);
    
    _add_bra_trafo_call_tree(lines, geom_integrals, bra_base_integrals, bra_rec_base_integrals, ket_base_integrals, ket_rec_base_integrals, integral);
    
    _add_loop_end(lines, integral);
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CGeomFuncBodyDriver::_get_gto_pairs_def() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// intialize GTOs pair data on bra side");

    vstr.push_back("const auto a_coords = bra_gto_pair_block.bra_coordinates();");
       
    vstr.push_back("const auto b_coords = bra_gto_pair_block.ket_coordinates();");
    
    vstr.push_back("const auto a_vec_exps = bra_gto_pair_block.bra_exponents();");
        
    vstr.push_back("const auto b_vec_exps = bra_gto_pair_block.ket_exponents();");
        
    vstr.push_back("const auto ab_vec_norms = bra_gto_pair_block.normalization_factors();");
        
    vstr.push_back("const auto ab_vec_ovls = bra_gto_pair_block.overlap_factors();");
        
    vstr.push_back("const auto a_indices = bra_gto_pair_block.bra_orbital_indices();");
        
    vstr.push_back("const auto b_indices = bra_gto_pair_block.ket_orbital_indices();");
      
    vstr.push_back("const auto bra_ncgtos = bra_gto_pair_block.number_of_contracted_pairs();");
        
    vstr.push_back("const auto bra_npgtos = bra_gto_pair_block.number_of_primitive_pairs();");
        
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
T4CGeomFuncBodyDriver::_get_ket_variables_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned 2D arrays for ket side");
  
    // c_exps, d_exps, cd_ovls, cd_norms, c_coords, d_coords, q_coords, pq_coords, f_ss
    
    size_t nelems = 17;
    
    if (_need_center_w(integral)) nelems += 3;
    
    if (_need_distances_qd(integral)) nelems += 3;
    
    if (_need_distances_wq(integral)) nelems += 3;
    
    if (_need_distances_wp(integral)) nelems += 3;
        
    vstr.push_back("CSimdArray<double> pfactors(" + std::to_string(nelems) +  ", ket_npgtos);");
  
    if (_need_hrr_for_ket(integral))
    {
        vstr.push_back("CSimdArray<double> cfactors(9, 1);");
    }
    
    return vstr;
}

bool
T4CGeomFuncBodyDriver::_need_center_w(const I4CIntegral& integral) const
{
    const auto orders = integral.prefixes_order();
    
    if (orders.empty())
    {
        return (integral[0] + integral[1] + integral[2] + integral[3]) > 0;
    }
    else
    {
        return (integral[0] + integral[1] + integral[2] + integral[3] + orders[0] + orders[1] + orders[2] + orders[3]) > 0;
    }
}

bool
T4CGeomFuncBodyDriver::_need_distances_qd(const I4CIntegral& integral) const
{
    const auto orders = integral.prefixes_order();
    
    if (orders.empty())
    {
        return (integral[2] + integral[3]) > 0;
    }
    else
    {
        return (integral[2] + integral[3] + orders[2] + orders[3]) > 0;
    }
}

bool
T4CGeomFuncBodyDriver::_need_distances_wq(const I4CIntegral& integral) const
{
    const auto orders = integral.prefixes_order();
    
    if (orders.empty())
    {
        return (integral[2] + integral[3]) > 0;
    }
    else
    {
        return (integral[2] + integral[3] + orders[2] + orders[3]) > 0;
    }
}

bool
T4CGeomFuncBodyDriver::_need_distances_wp(const I4CIntegral& integral) const
{
    const auto orders = integral.prefixes_order();
    
    if (orders.empty())
    {
        return (integral[0] + integral[1]) > 0;
    }
    else
    {
        return (integral[0] + integral[1] + orders[0] + orders[1]) > 0;
    }
}

bool
T4CGeomFuncBodyDriver::_need_hrr_for_ket(const I4CIntegral& integral) const
{
    const auto orders = integral.prefixes_order();
    
    if (orders.empty())
    {
        return integral[2] > 0;
    }
    else
    {
        return (integral[2] + orders[2]) > 0;
    }
}

bool
T4CGeomFuncBodyDriver::_need_hrr_for_bra(const I4CIntegral& integral) const
{
    const auto orders = integral.prefixes_order();
    
    if (orders.empty())
    {
        return integral[0] > 0;
    }
    else
    {
        return (integral[0] + orders[0]) > 0;
    }
}

size_t
T4CGeomFuncBodyDriver::_get_index_w(const I4CIntegral& integral) const
{
    return 17;
}

size_t
T4CGeomFuncBodyDriver::_get_index_qd(const I4CIntegral& integral) const
{
    auto index = _get_index_w(integral);
    
    if (_need_center_w(integral)) index += 3;
    
    return index;
}

size_t
T4CGeomFuncBodyDriver::_get_index_wq(const I4CIntegral& integral) const
{
    auto index = _get_index_qd(integral);
    
    if (_need_distances_qd(integral)) index += 3;
    
    return index;
}

size_t
T4CGeomFuncBodyDriver::_get_index_wp(const I4CIntegral& integral) const
{
    auto index = _get_index_wq(integral);
    
    if (_need_distances_wq(integral)) index += 3;
    
    return index;
}

size_t
T4CGeomFuncBodyDriver::_get_index(const size_t         start,
                                  const I4CIntegral&   integral,
                                  const SI4CIntegrals& integrals) const
{
    size_t index = start;
    
    for (const auto& tint : integrals)
    {
        if (tint == integral) return index;
        
        index += tint.components<T2CPair, T2CPair>().size();
    }
    
    return 0;
}

size_t
T4CGeomFuncBodyDriver::_get_half_spher_index(const size_t         start,
                                             const I4CIntegral&   integral,
                                             const SI4CIntegrals& integrals) const
{
    size_t index = start;
    
    for (const auto& tint : integrals)
    {
        if (tint == integral) return index;
        
        auto angpair = std::array<int, 2>({tint[2], tint[3]});
                
        auto icomps = t2c::number_of_spherical_components(angpair);
            
        angpair = std::array<int, 2>({tint[0], tint[1]});
                
        icomps *= t2c::number_of_cartesian_components(angpair);
        
        index += icomps;
    }
    
    return index;
}

size_t
T4CGeomFuncBodyDriver::_get_geom_half_spher_index(const size_t         start,
                                                  const I4CIntegral&   integral,
                                                  const SI4CIntegrals& integrals) const
{
    size_t index = start;
    
    for (const auto& tint : integrals)
    {
        if (tint == integral) return index;
        
        auto angpair = std::array<int, 2>({tint[2], tint[3]});
                
        auto icomps = t2c::number_of_spherical_components(angpair);
            
        angpair = std::array<int, 2>({tint[0], tint[1]});
                
        icomps *= t2c::number_of_cartesian_components(angpair);
        
        for (const auto& prefix : tint.prefixes())
        {
            icomps *= prefix.components().size();
        }
        
        index += icomps;
    }
    
    return index;
}

R4Group
T4CGeomFuncBodyDriver::_generate_integral_group(const VT4CIntegrals& components,
                                                const I4CIntegral&   integral) const
{
    R4Group rgroup;
        
    T4CCenterDriver t4c_geom_drv;
        
    return t4c_geom_drv.create_recursion(components);
}

std::vector<std::string>
T4CGeomFuncBodyDriver::_get_buffers_str(const SI4CIntegrals& geom_integrals,
                                        const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    for (const auto& tint : geom_integrals)
    {
        auto label = t4c::get_geom_buffer_label(tint);
        
        vstr.push_back("/// Set up components of auxilary buffer : " + label);
        
        const auto tlabel = _get_tensor_label(tint);
        
        int index = 0;
        
        for (const auto& tcomp : tint.components<T2CPair, T2CPair>())
        {
            const auto line = "auto " + _get_component_label(tcomp) + " = " + label;
                
            vstr.push_back(line + "[" + std::to_string(index) + "];");
            
            index++;
        }
    }
    
    auto label = t4c::get_geom_buffer_label(integral);
    
    vstr.push_back("/// Set up components of integrals buffer : " + label);
    
    int index = 0;
    
    for (const auto& tcomp : integral.components<T2CPair, T2CPair>())
    {
        const auto line = "auto " + _get_component_label(tcomp) + " = " + label;
            
        vstr.push_back(line + "[" + std::to_string(index) + "];");
        
        index++;
    }
    
    return vstr;
}

std::string
T4CGeomFuncBodyDriver::_get_tensor_label(const I4CIntegral& integral) const
{
    std::string label = "g";

    return label;
}

std::string
T4CGeomFuncBodyDriver::_get_tensor_label(const T4CIntegral& integral) const
{
    std::string label = "g";
    
    return label;
}

void
T4CGeomFuncBodyDriver::_add_recursion_loop(      VCodeLines&         lines,
                                           const R4Group&            rgroup,
                                           const I4CIntegral&        integral,
                                           const std::array<int, 2>& rec_range) const
{
    const auto var_str = _get_pragma_str(rgroup, integral, rec_range);
    
    lines.push_back({1, 0, 2, "// integrals block (" + std::to_string(rec_range[0]) + "-" +  std::to_string(rec_range[1]) + ")"});
    
    lines.push_back({1, 0, 1, "#pragma omp simd aligned(" + var_str + " : 64)"});
    
    lines.push_back({1, 0, 1, "for (size_t i = 0; i < ndims; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    // _get_factor_lines(lines, rec_dists);
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        if (i < (rec_range[1] - 1))
        {
            lines.push_back({2, 0, 2, _get_code_line(rgroup[i])});
        }
        else
        {
            lines.push_back({2, 0, 1, _get_code_line(rgroup[i])});
        }
    }
    
    lines.push_back({1, 0, 1, "}"});
}

std::string
T4CGeomFuncBodyDriver::_get_pragma_str(const R4Group&            rgroup,
                                       const I4CIntegral&        integral,
                                       const std::array<int, 2>& rec_range) const
{
    std::set<std::string> tlabels;
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        tlabels.insert(_get_component_label(rgroup[i].root().integral()));
        
        for (size_t j = 0; j < rgroup[i].terms(); j++)
        {
            tlabels.insert(_get_component_label(rgroup[i][j].integral().base()));
        }
    }
        
    std::string label;
    
    for (const auto& tlabel : tlabels)
    {
        label += tlabel + ", ";
    }
    
    if (const auto prefixes = integral.prefixes(); !prefixes.empty())
    {
        if (prefixes[2].shape().order() > 0) label += "c_exps, ";
        
        if (prefixes[3].shape().order() > 0) label += "d_exps";
    }
    
    if (label[label.size() - 2] == ',') label.erase(label.end() - 2);
    
    return label;
}

std::string
T4CGeomFuncBodyDriver::_get_code_line(const R4CDist& rec_distribution) const
{
    auto tint = rec_distribution.root().integral();
    
    std::string line = _get_component_label(tint) + "[i] = ";
    
    for (size_t i = 0; i < rec_distribution.terms(); i++)
    {
        line += _get_rterm_code(rec_distribution[i], i == 0);
    }
    
    return line + ";";
}

std::string
T4CGeomFuncBodyDriver::_get_rterm_code(const R4CTerm& rec_term,
                                       const bool     is_first) const
{
    const auto pre_fact = rec_term.prefactor();
        
    auto plabel = pre_fact.label();
        
    if (plabel == "1.0")  plabel = "";
        
    if (plabel == "-1.0") plabel = "-";
    
    if (plabel.size() > 1) plabel += " * ";
    
    auto tint = rec_term.integral().base();
    
    plabel += _get_component_label(tint) + "[i]";
        
    for (const auto& fact : rec_term.factors())
    {
        for (int i = 0; i < rec_term.factor_order(fact); i++)
        {
            plabel+= " * " + fact.label();
            
            plabel.erase(plabel.end() - 2, plabel.end());
            
            if (fact.label() == "c_exps_0") plabel += "[i]";
            
            if (fact.label() == "d_exps_0") plabel += "[i]";
        }
    }
    
    if (!is_first)
    {
        if (plabel[0] == '-')
        {
            plabel.insert(plabel.begin() + 1, 1, ' ');
            
            plabel = " " + plabel;
        }
        else
        {
            plabel = " + " + plabel;
        }
    }
        
    return plabel;
}

std::string
T4CGeomFuncBodyDriver::_get_component_label(const T4CIntegral& integral) const
{
    std::string label = _get_tensor_label(integral) + "_" + integral.label();
    
    return label;
}

std::vector<std::string>
T4CGeomFuncBodyDriver::_get_prim_buffers_def(const SI4CIntegrals& integrals,
                                             const I4CIntegral&   integral) const
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
T4CGeomFuncBodyDriver::_get_cart_buffers_def(const SI4CIntegrals& bra_base_integrals,
                                             const SI4CIntegrals& bra_rec_base_integrals,
                                             const SI4CIntegrals& ket_base_integrals,
                                             const SI4CIntegrals& ket_rec_base_integrals,
                                             const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned Cartesian integrals");
    
    auto tcomps = _get_all_components(_get_cart_buffer_integrals(bra_base_integrals, ket_base_integrals));
    
    tcomps += _get_all_components(_get_cart_buffer_integrals(bra_rec_base_integrals, ket_rec_base_integrals));
    
    const auto geom_orders = integral.prefixes_order();
    
    if (geom_orders == std::vector<int>({2, 0, 0, 0}))
    {
        tcomps += _get_geom20_cart_2a_size(bra_rec_base_integrals, ket_rec_base_integrals, integral);
    }
    
    std::string label = "CSimdArray<double> cbuffer";
    
    label += "(" + std::to_string(tcomps) + ", 1);";
    
    vstr.push_back(label);
    
    return vstr;
}

std::vector<std::string>
T4CGeomFuncBodyDriver::_get_contr_buffers_def(const SI4CIntegrals& ket_base_integrals,
                                              const SI4CIntegrals& ket_rec_base_integrals) const
{
    std::vector<std::string> vstr;
    
    auto tcomps = _get_all_components(_get_contr_buffers_integrals(ket_base_integrals));
    
    tcomps += _get_all_components(_get_contr_buffers_integrals(ket_rec_base_integrals));
    
    if (tcomps > 0)
    {
        vstr.insert(vstr.begin(), "// allocate aligned contracted integrals");
        
        std::string label = "CSimdArray<double> ";
        
        label += "ckbuffer(" + std::to_string(tcomps) + ", 1);";
        
        vstr.push_back(label);
    }
    
    return vstr;
}

SI4CIntegrals
T4CGeomFuncBodyDriver::_get_cart_buffer_integrals(const SI4CIntegrals& bra_integrals,
                                                  const SI4CIntegrals& ket_integrals) const
{
    SI4CIntegrals tints;
    
    for (const auto& tint : bra_integrals)
    {
        if ((tint[0] + tint[2]) == 0)
        {
            tints.insert(tint);
        }
    }
    
    for (const auto& tint : ket_integrals)
    {
        if ((tint[0] + tint[2]) == 0)
        {
            tints.insert(tint);
        }
    }
    
    return tints;
}

SI4CIntegrals
T4CGeomFuncBodyDriver::_get_contr_buffers_integrals(const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;

    for (const auto& tint : integrals)
    {
        if ((tint[0] == 0) && (tint[2] > 0))
        {
            tints.insert(tint);
        }
    }
    
    return tints;
}

std::vector<std::string>
T4CGeomFuncBodyDriver::_get_half_spher_buffers_def(const SI4CIntegrals& geom_integrals, 
                                                   const SI4CIntegrals& bra_base_integrals,
                                                   const SI4CIntegrals& bra_rec_base_integrals,
                                                   const SI4CIntegrals& ket_base_integrals,
                                                   const SI4CIntegrals& ket_rec_base_integrals,
                                                   const I4CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned half transformed integrals");
    
    const auto geom_orders = integral.prefixes_order();
    
    size_t tcomps = 0;
    
    if (geom_orders == std::vector<int>({1, 0, 0, 0}))
    {
        tcomps += _get_all_half_spher_components(_get_half_spher_buffers_integrals(bra_base_integrals, ket_base_integrals, integral));
        
        tcomps += _get_all_half_spher_components(_get_half_spher_buffers_integrals(bra_rec_base_integrals, ket_rec_base_integrals, integral));
        
        if (integral[0] > 0)
        {
            tcomps += _get_all_geom_half_spher_components(_get_geom_half_spher_buffers_integrals(geom_integrals, integral));
        }
    }
    
    if (geom_orders == std::vector<int>({2, 0, 0, 0}))
    {
        tcomps += _get_all_half_spher_components(_get_half_spher_buffers_integrals(bra_rec_base_integrals, ket_rec_base_integrals, integral));
        
        tcomps += _get_geom20_half_spher_2a_size(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        if (integral[0] == 0)
        {
            tcomps += _get_geom20_half_spher_size(integral);
        }
    }
    
    std::string label = "CSimdArray<double> ";
            
    label += "skbuffer(" + std::to_string(tcomps) + ", 1);";
            
    vstr.push_back(label);
        
    return vstr;
}

std::vector<std::string>
T4CGeomFuncBodyDriver::_get_spher_buffers_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned spherical integrals");
                    
    const auto tcomps = _get_all_spher_components(integral);
    
    std::string label = "CSimdArray<double> ";
                    
    label += "sbuffer(" + std::to_string(tcomps) + ", 1);";
                    
    vstr.push_back(label);
   
    return vstr;
}

std::vector<std::string>
T4CGeomFuncBodyDriver::_get_boys_function_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    auto order = integral[0] + integral[1] + integral[2] + integral[3];
    
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
T4CGeomFuncBodyDriver::_add_loop_start(      VCodeLines&    lines,
                                       const I4CIntegral&   integral) const
{
    lines.push_back({1, 0, 2, "// set up ket partitioning"});

    lines.push_back({1, 0, 2, "const auto ket_dim = ket_indices.second - ket_indices.first;"});

    lines.push_back({1, 0, 2, "const auto ket_blocks = batch::number_of_batches(ket_dim, simd::width<double>());"});

    lines.push_back({1, 0, 1, "for (size_t i = 0; i < ket_blocks; i++)"});
                    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "auto ket_range = batch::batch_range(i, ket_dim, simd::width<double>(), ket_indices.first);"});

    lines.push_back({2, 0, 2, "pfactors.load(c_vec_exps, ket_range, 0, ket_npgtos);"});
    
    lines.push_back({2, 0, 2, "pfactors.load(d_vec_exps, ket_range, 1, ket_npgtos);"});
 
    lines.push_back({2, 0, 2, "pfactors.load(cd_vec_ovls, ket_range, 2, ket_npgtos);"});
    
    lines.push_back({2, 0, 2, "pfactors.load(cd_vec_norms, ket_range, 3, ket_npgtos);"});
    
    lines.push_back({2, 0, 2, "pfactors.replicate_points(c_coords, ket_range, 4, ket_npgtos);"});
    
    lines.push_back({2, 0, 2, "pfactors.replicate_points(d_coords, ket_range, 7, ket_npgtos);"});
    
    if (_need_hrr_for_ket(integral))
    {
        lines.push_back({2, 0, 2, "cfactors.replicate_points(c_coords, ket_range, 0, 1);"});
        
        lines.push_back({2, 0, 2, "cfactors.replicate_points(d_coords, ket_range, 3, 1);"});
        
        lines.push_back({2, 0, 2, "t4cfunc::comp_distances_cd(cfactors, 6, 0, 3);"});
    }
   
    lines.push_back({2, 0, 2, "// set up active SIMD width"});
    
    lines.push_back({2, 0, 2, "const auto ket_width = ket_range.second - ket_range.first;"});
    
    lines.push_back({2, 0, 2, "pbuffer.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "cbuffer.set_active_width(ket_width);"});
    
    if (_need_hrr_for_ket(integral))
    {
        lines.push_back({2, 0, 2, "ckbuffer.set_active_width(ket_width);"});
    }
    
    lines.push_back({2, 0, 2, "skbuffer.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "sbuffer.set_active_width(ket_width);"});
    
    lines.push_back({2, 0, 2, "bf_data.set_active_width(ket_width);"});
      
    lines.push_back({2, 0, 2, "// loop over basis function pairs on bra side"});

    lines.push_back({2, 0, 1, "for (auto j = bra_indices.first; j < bra_indices.second; j++)"});

    lines.push_back({2, 0, 1, "{"});
    
    lines.push_back({3, 0, 2, "// zero integral buffers"});
    
    lines.push_back({3, 0, 2, "cbuffer.zero();"});
    
    if (_need_hrr_for_ket(integral))
    {
        lines.push_back({3, 0, 2, "ckbuffer.zero();"});
    }
    
    lines.push_back({3, 0, 2, "skbuffer.zero();"});
    
    lines.push_back({3, 0, 2, "sbuffer.zero();"});

    lines.push_back({3, 0, 2, "// set up coordinates on bra side"});

    lines.push_back({3, 0, 2, "const auto r_a = a_coords[j];"});
    
    lines.push_back({3, 0, 2, "const auto r_b = b_coords[j];"});
    
    lines.push_back({3, 0, 2, "const auto a_xyz = r_a.coordinates();"});
    
    lines.push_back({3, 0, 2, "const auto b_xyz = r_b.coordinates();"});

    if (_need_hrr_for_bra(integral))
    {
        lines.push_back({3, 0, 2, "const auto r_ab = TPoint<double>({a_xyz[0] - b_xyz[0], a_xyz[1] - b_xyz[1], a_xyz[2] - b_xyz[2]});"});
    }
}

void
T4CGeomFuncBodyDriver::_add_loop_end(      VCodeLines&  lines,
                                     const I4CIntegral& integral) const
{
    lines.push_back({2, 0, 1, "}"});
   
    lines.push_back({1, 0, 2, "}"});
}

void
T4CGeomFuncBodyDriver::_add_ket_loop_start(      VCodeLines&  lines,
                                           const I4CIntegral& integral) const
{
    auto geom_orders = integral.prefixes_order();
    
    lines.push_back({3, 0, 1, "for (int k = 0; k < bra_npgtos; k++)"});
   
    lines.push_back({3, 0, 1, "{"});
    
    lines.push_back({4, 0, 2, "const auto a_exp = a_vec_exps[k * bra_ncgtos + j];"});
        
    lines.push_back({4, 0, 2, "const auto b_exp = b_vec_exps[k * bra_ncgtos + j];"});
            
    lines.push_back({4, 0, 2, "const auto ab_norm = ab_vec_norms[k * bra_ncgtos + j];"});
        
    lines.push_back({4, 0, 2, "const auto ab_ovl = ab_vec_ovls[k * bra_ncgtos + j];"});
    
    lines.push_back({4, 0, 2, "const auto p_x = (a_xyz[0] * a_exp + b_xyz[0] * b_exp) / (a_exp + b_exp);"});
    
    lines.push_back({4, 0, 2, "const auto p_y = (a_xyz[1] * a_exp + b_xyz[1] * b_exp) / (a_exp + b_exp);"});
    
    lines.push_back({4, 0, 2, "const auto p_z = (a_xyz[2] * a_exp + b_xyz[2] * b_exp) / (a_exp + b_exp);"});
    
    lines.push_back({4, 0, 2, "const auto r_p = TPoint<double>({p_x, p_y, p_z});"});
    
    if ((integral[0] + integral[1] + geom_orders[0] + geom_orders[1]) > 0)
    {
        lines.push_back({4, 0, 2, "const auto pb_x = p_x - b_xyz[0];"});
        
        lines.push_back({4, 0, 2, "const auto pb_y = p_y - b_xyz[1];"});
        
        lines.push_back({4, 0, 2, "const auto pb_z = p_z - b_xyz[2];"});
        
        lines.push_back({4, 0, 2, "const auto r_pb = TPoint<double>({pb_x, pb_y, pb_z});"});
    }
    
    lines.push_back({4, 0, 2, "t4cfunc::comp_coordinates_q(pfactors, 10, 4, 7);"});
    
    lines.push_back({4, 0, 2, "t4cfunc::comp_distances_pq(pfactors, 13, 10, r_p);"});
    
    if (_need_center_w(integral))
    {
        auto label_w = std::to_string(_get_index_w(integral));
        
        lines.push_back({4, 0, 2, "t4cfunc::comp_coordinates_w(pfactors, " + label_w + ", 10, r_p, a_exp, b_exp);"});
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
     
    if (_need_distances_wp(integral))
    {
        auto label_w = std::to_string(_get_index_w(integral));
        
        auto label_wp = std::to_string(_get_index_wp(integral));
        
        lines.push_back({4, 0, 2, "t4cfunc::comp_distances_wp(pfactors, " + label_wp + ", " + label_w + ", r_p);"});
    }
    
    auto border = integral[0] + integral[1] + integral[2] + integral[3] + 1;
    
    border += geom_orders[0] + geom_orders[1] + geom_orders[2] + geom_orders[3]; 
    
    lines.push_back({4, 0, 2, "t4cfunc::comp_boys_args(bf_data, " + std::to_string(border) + ", pfactors, 13, a_exp, b_exp);"});
    
    lines.push_back({4, 0, 2, "bf_table.compute(bf_data, 0, " + std::to_string(border) + ");"});
    
    lines.push_back({4, 0, 2, "t4cfunc::comp_ovl_factors(pfactors, 16, 2, 3, ab_ovl, ab_norm, a_exp, b_exp);"});
}

void
T4CGeomFuncBodyDriver::_add_ket_loop_end(      VCodeLines&  lines,
                                         const SI4CIntegrals& bra_base_integrals,
                                         const SI4CIntegrals& bra_rec_base_integrals,
                                         const SI4CIntegrals& ket_base_integrals,
                                         const SI4CIntegrals& ket_rec_base_integrals,
                                         const SI4CIntegrals& vrr_integrals,
                                         const I4CIntegral& integral) const
{
    if (integral.prefixes_order() == std::vector<int>{1, 0, 0, 0})
    {
        auto cints = _get_cart_buffer_integrals(bra_base_integrals, ket_base_integrals);
        
        // non scaled integrals
        
        for (const auto& tint : cints)
        {
            if ((tint[0] + tint[2]) == 0)
            {
                std::string label = "t2cfunc::reduce(cbuffer, ";
                
                label +=  std::to_string(_get_index(0, tint, cints)) + ", ";
                
                label += "pbuffer, ";
                
                label += std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
                
                label += std::to_string(tint.components<T2CPair, T2CPair>().size()) + ", ";
                
                label += "ket_width, ket_npgtos);";
                
                lines.push_back({4, 0, 2, label});
            }
        }
        
        auto tcomps = _get_all_components(cints);
        
        // scaling of integrals
        
        cints = _get_cart_buffer_integrals(bra_rec_base_integrals, ket_rec_base_integrals);
        
        for (const auto& tint : cints)
        {
            if ((tint[0] + tint[2]) == 0)
            {
                std::string label = "pbuffer.scale(2.0 * a_exp, {";
                
                label +=  std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
                
                label +=  std::to_string(_get_index(0, tint, vrr_integrals) + tint.components<T2CPair, T2CPair>().size()) + "});";
                
                lines.push_back({4, 0, 2, label});
            }
        }
        
        // scaled integrals
        
        for (const auto& tint : cints)
        {
            if ((tint[0] + tint[2]) == 0)
            {
                std::string label = "t2cfunc::reduce(cbuffer, ";
                
                label +=  std::to_string(_get_index(tcomps, tint, cints)) + ", ";
                
                label += "pbuffer, ";
                
                label += std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
                
                label += std::to_string(tint.components<T2CPair, T2CPair>().size()) + ", ";
                
                label += "ket_width, ket_npgtos);";
                
                lines.push_back({4, 0, 2, label});
            }
        }
        
        lines.push_back({3, 0, 2, "}"});
    }
    
    if (integral.prefixes_order() == std::vector<int>{2, 0, 0, 0})
    {
        auto cints = _get_cart_buffer_integrals(bra_rec_base_integrals, ket_rec_base_integrals);
        
        // scaling of integrals (2 alpha)
        
        for (const auto& tint : cints)
        {
            if (((tint[0] + tint[2]) == 0) && (tint[1] == integral[1]))
            {
                std::string label = "pbuffer.scale(2.0 * a_exp, {";
                
                label +=  std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
                
                label +=  std::to_string(_get_index(0, tint, vrr_integrals) + tint.components<T2CPair, T2CPair>().size()) + "});";
                
                lines.push_back({4, 0, 2, label});
            }
        }
        
        // reduction of scaled integrals (2 alpha)
        
        for (const auto& tint : cints)
        {
            if (((tint[0] + tint[2]) == 0) && (tint[1] == integral[1]))
            {
                std::string label = "t2cfunc::reduce(cbuffer, ";
                
                label +=  std::to_string(_get_index(0, tint, cints)) + ", ";
                
                label += "pbuffer, ";
                
                label += std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
                
                label += std::to_string(tint.components<T2CPair, T2CPair>().size()) + ", ";
                
                label += "ket_width, ket_npgtos);";
                
                lines.push_back({4, 0, 2, label});
            }
        }
        
        // scaling of integrals (4 alpha)
        
        for (const auto& tint : cints)
        {
            if ((tint[0] + tint[2]) == 0)
            {
                std::string label;
                
                if (tint[1] == integral[1])
                {
                    label += "pbuffer.scale(2.0 * a_exp, {";
                }
                else
                {
                    label += "pbuffer.scale(4.0 * a_exp * a_exp, {";
                }
 
                label +=  std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
                
                label +=  std::to_string(_get_index(0, tint, vrr_integrals) + tint.components<T2CPair, T2CPair>().size()) + "});";
                
                lines.push_back({4, 0, 2, label});
            }
        }
        
        // reduction of scaled integrals (4 alpha)
        
        const auto cstart = _get_geom20_cart_2a_size(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        for (const auto& tint : cints)
        {
            if ((tint[0] + tint[2]) == 0)
            {
                std::string label = "t2cfunc::reduce(cbuffer, ";
                
                label +=  std::to_string(_get_index(cstart, tint, cints)) + ", ";
                
                label += "pbuffer, ";
                
                label += std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
                
                label += std::to_string(tint.components<T2CPair, T2CPair>().size()) + ", ";
                
                label += "ket_width, ket_npgtos);";
                
                lines.push_back({4, 0, 2, label});
            }
        }
        
        lines.push_back({3, 0, 2, "}"});
    }
}

void
T4CGeomFuncBodyDriver::_add_auxilary_integrals(      VCodeLines&    lines,
                                               const SI4CIntegrals& integrals,
                                               const I4CIntegral&   integral,
                                               const size_t         spacer) const
{
    for (const auto& tint : integrals)
    {
        if ((tint[0] + tint[1] + tint[2] + tint[3]) == 0)
        {
            const auto blabel = std::to_string(tint.order());
            
            const auto ilabel = std::to_string(_get_index(0, tint, integrals));
                    
            lines.push_back({spacer, 0, 2, "erirec::comp_prim_electron_repulsion_ssss(pbuffer, " + ilabel + ", pfactors, 16, bf_data, " + blabel + ");"});
        }
    }
}

void
T4CGeomFuncBodyDriver::_add_vrr_call_tree(      VCodeLines&  lines,
                                          const SI4CIntegrals& integrals,
                                          const I4CIntegral&   integral,
                                          const size_t         spacer) const
{
    for (const auto& tint : integrals)
    {
        if (((tint[0] + tint[2]) == 0) && ((tint[1] + tint[3]) > 0))
        {
            const auto name = t4c::prim_compute_func_name(tint);
            
            auto label = t4c::namespace_label(tint) + "::" + name + "(pbuffer, ";
            
            label += _get_vrr_arguments(0, integrals, tint);
            
            label += "pfactors, ";
            
            if (_need_distances_wp(tint))
            {
                label += std::to_string(_get_index_wp(integral)) + ", r_pb, ";
            }
            else
            {
                label += std::to_string(_get_index_qd(integral)) + ", ";
                
                label += std::to_string(_get_index_wq(integral)) + ", ";
            }
            
            if ((tint[1] + tint[3]) > 1)
            {
                label += "a_exp, b_exp";
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
T4CGeomFuncBodyDriver::_add_ket_hrr_call_tree(      VCodeLines&  lines,
                                              const SI4CIntegrals& bra_base_integrals,
                                              const SI4CIntegrals& bra_rec_base_integrals,
                                              const SI4CIntegrals& ket_base_integrals,
                                              const SI4CIntegrals& ket_rec_base_integrals,
                                              const I4CIntegral&   integral,
                                              const size_t         spacer) const
{
    const auto geom_orders = integral.prefixes_order();
    
    if (geom_orders == std::vector<int>({1, 0, 0, 0}))
    {
        auto ckints = _get_contr_buffers_integrals(ket_base_integrals);
        
        for (const auto& tint : ckints)
        {
            if ((tint[0] == 0) && (tint[2] > 0))
            {
                const auto name = t4c::ket_hrr_compute_func_name(tint);
                
                auto label = t4c::namespace_label(tint) + "::" + name + "(ckbuffer, ";
                
                label += std::to_string(_get_index(0, tint, ckints)) + ", ";
                
                if (tint[2] == 1)
                {
                    label += "cbuffer, ";
                }
                
                label += _get_ket_hrr_arguments(0, tint, bra_base_integrals, ket_base_integrals);
                
                label += "cfactors, 6, ";
                
                label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]);
                
                label += ");";
                
                lines.push_back({spacer, 0, 2, label});
            }
        }
        
        const auto ccomps = _get_all_components(_get_cart_buffer_integrals(bra_base_integrals, ket_base_integrals));
        
        auto ckcomps = _get_all_components(ckints);
        
        ckints = _get_contr_buffers_integrals(ket_rec_base_integrals);
        
        for (const auto& tint : ckints)
        {
            if ((tint[0] == 0) && (tint[2] > 0))
            {
                const auto name = t4c::ket_hrr_compute_func_name(tint);
                
                auto label = t4c::namespace_label(tint) + "::" + name + "(ckbuffer, ";
                
                label += std::to_string(_get_index(ckcomps, tint, ckints)) + ", ";
                
                if (tint[2] == 1)
                {
                    label += "cbuffer, ";
                    
                    label += _get_ket_hrr_arguments(ccomps, tint, bra_rec_base_integrals, ket_rec_base_integrals);
                }
                else
                {
                    label += _get_ket_hrr_arguments(ckcomps, tint, bra_rec_base_integrals, ket_rec_base_integrals);
                }
                
                label += "cfactors, 6, ";
                
                label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]);
                
                label += ");";
                
                lines.push_back({spacer, 0, 2, label});
            }
        }
    }
    
    if (geom_orders == std::vector<int>({2, 0, 0, 0}))
    {
        auto ckints = _get_contr_buffers_integrals(ket_rec_base_integrals);
        
        for (const auto& tint : ckints)
        {
            if ((tint[0] == 0) && (tint[1] == integral[1]) && (tint[2] > 0))
            {
                const auto name = t4c::ket_hrr_compute_func_name(tint);
                
                auto label = t4c::namespace_label(tint) + "::" + name + "(ckbuffer, ";
                
                label += std::to_string(_get_index(0, tint, ckints)) + ", ";
                
                if (tint[2] == 1)
                {
                    label += "cbuffer, ";
                    
                    label += _get_ket_hrr_arguments(0, tint, bra_rec_base_integrals, ket_rec_base_integrals);
                }
                else
                {
                    label += _get_ket_hrr_arguments(0, tint, bra_rec_base_integrals, ket_rec_base_integrals);
                }
                
                label += "cfactors, 6, ";
                
                label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]);
                
                label += ");";
                
                lines.push_back({spacer, 0, 2, label});
            }
        }
        
        const auto cstart = _get_geom20_cart_2a_size(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        const auto ckstart = _get_geom20_contr_2a_size(ket_rec_base_integrals, integral);
        
        for (const auto& tint : ckints)
        {
            if ((tint[0] == 0) && (tint[2] > 0))
            {
                const auto name = t4c::ket_hrr_compute_func_name(tint);
                
                auto label = t4c::namespace_label(tint) + "::" + name + "(ckbuffer, ";
                
                label += std::to_string(_get_index(ckstart, tint, ckints)) + ", ";
                
                if (tint[2] == 1)
                {
                    label += "cbuffer, ";
                    
                    label += _get_ket_hrr_arguments(cstart, tint, bra_rec_base_integrals, ket_rec_base_integrals);
                }
                else
                {
                    label += _get_ket_hrr_arguments(ckstart, tint, bra_rec_base_integrals, ket_rec_base_integrals);
                }
                
                label += "cfactors, 6, ";
                
                label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]);
                
                label += ");";
                
                lines.push_back({spacer, 0, 2, label});
            }
        }
    }
    
}

void
T4CGeomFuncBodyDriver::_add_ket_trafo_call_tree(      VCodeLines&  lines,
                                                const SI4CIntegrals& bra_base_integrals,
                                                const SI4CIntegrals& bra_rec_base_integrals,
                                                const SI4CIntegrals& ket_base_integrals,
                                                const SI4CIntegrals& ket_rec_base_integrals,
                                                const I4CIntegral&   integral,
                                                const size_t         spacer) const
{
    const auto geom_orders = integral.prefixes_order();
    
    if (geom_orders == std::vector<int>({1, 0, 0, 0}))
    {
        auto skints = _get_half_spher_buffers_integrals(bra_base_integrals, ket_base_integrals, integral);
        
        auto ckints = _get_contr_buffers_integrals(ket_base_integrals);
        
        auto cints = _get_cart_buffer_integrals(bra_base_integrals, ket_base_integrals);
        
        // non scaled integrals
        
        if (integral[2] > 0)
        {
            for (const auto& tint : ket_base_integrals)
            {
                if ((tint[0] == 0) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
                {
                    std::string label = "t4cfunc::ket_transform<" + std::to_string(tint[2]) + ", " + std::to_string(tint[3]) + ">";
                    
                    label += "(skbuffer, "  + std::to_string(_get_half_spher_index(0, tint, skints)) + ", ";
                    
                    label += "ckbuffer, " + std::to_string(_get_index(0, tint, ckints))  + ", ";
                    
                    label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]) + ");";
                    
                    lines.push_back({spacer, 0, 2, label});
                }
            }
        }
        
        if ((integral[0] >= 0) && (integral[2] == 0))
        {
            for (const auto& tint : cints)
            {
                if ((tint[0] == 0) && (tint[2] == 0))
                {
                    std::string label = "t4cfunc::ket_transform<" + std::to_string(tint[2]) + ", " + std::to_string(tint[3]) + ">";
                        
                    label += "(skbuffer, "  +  std::to_string(_get_half_spher_index(0, tint, skints)) + ", ";
                    
                    label += "cbuffer, " + std::to_string(_get_index(0, tint, cints))  + ", ";
                    
                    label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]) + ");";
                        
                    lines.push_back({spacer, 0, 2, label});
                }
            }
        }
        
        // scaled integrals
        
        const auto skcomps = _get_all_half_spher_components(skints);
        
        const auto ckcomps = _get_all_components(ckints);
        
        const auto ccomps = _get_all_components(cints);
        
        skints = _get_half_spher_buffers_integrals(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        ckints = _get_contr_buffers_integrals(ket_rec_base_integrals);
        
        cints = _get_cart_buffer_integrals(bra_rec_base_integrals, ket_rec_base_integrals);
        
        if (integral[2] > 0)
        {
            for (const auto& tint : ket_rec_base_integrals)
            {
                if ((tint[0] == 0) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
                {
                    std::string label = "t4cfunc::ket_transform<" + std::to_string(tint[2]) + ", " + std::to_string(tint[3]) + ">";
                    
                    label += "(skbuffer, "  + std::to_string(_get_half_spher_index(skcomps, tint, skints)) + ", ";
                    
                    label += "ckbuffer, " + std::to_string(_get_index(ckcomps, tint, ckints))  + ", ";
                    
                    label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]) + ");";
                    
                    lines.push_back({spacer, 0, 2, label});
                }
            }
        }
        
        if ((integral[0] >= 0) && (integral[2] == 0))
        {
            for (const auto& tint : bra_rec_base_integrals)
            {
                if ((tint[0] == 0) && (tint[2] == 0))
                {
                    std::string label = "t4cfunc::ket_transform<" + std::to_string(tint[2]) + ", " + std::to_string(tint[3]) + ">";
                        
                    label += "(skbuffer, "  +  std::to_string(_get_half_spher_index(skcomps, tint, skints)) + ", ";
                    
                    label += "cbuffer, " + std::to_string(_get_index(ccomps, tint, cints))  + ", ";
                    
                    label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]) + ");";
                        
                    lines.push_back({spacer, 0, 2, label});
                }
            }
        }
    }
    
    if (geom_orders == std::vector<int>({2, 0, 0, 0}))
    {
        auto cints = _get_cart_buffer_integrals(bra_rec_base_integrals, ket_rec_base_integrals);
        
        auto ckints = _get_contr_buffers_integrals(ket_rec_base_integrals);
        
        auto skints = _get_half_spher_buffers_integrals(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        if (integral[2] > 0)
        {
            // scaled integrals (2 alpha)
            
            for (const auto& tint : ket_rec_base_integrals)
            {
                if ((tint[0] == 0) && (tint[1] == integral[1]) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
                {
                    std::string label = "t4cfunc::ket_transform<" + std::to_string(tint[2]) + ", " + std::to_string(tint[3]) + ">";
                    
                    label += "(skbuffer, "  + std::to_string(_get_half_spher_index(0, tint, skints)) + ", ";
                    
                    label += "ckbuffer, " + std::to_string(_get_index(0, tint, ckints))  + ", ";
                    
                    label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]) + ");";
                    
                    lines.push_back({spacer, 0, 2, label});
                }
            }
            
            // scaled integrals (4 alpha)
            
            const auto ckstart = _get_geom20_contr_2a_size(ket_rec_base_integrals, integral);
            
            const auto skstart = _get_geom20_half_spher_2a_size(bra_rec_base_integrals, ket_rec_base_integrals, integral);
            
            for (const auto& tint : ket_rec_base_integrals)
            {
                if ((tint[0] == 0) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
                {
                    std::string label = "t4cfunc::ket_transform<" + std::to_string(tint[2]) + ", " + std::to_string(tint[3]) + ">";
                    
                    label += "(skbuffer, "  + std::to_string(_get_half_spher_index(skstart, tint, skints)) + ", ";
                    
                    label += "ckbuffer, " + std::to_string(_get_index(ckstart, tint, ckints))  + ", ";
                    
                    label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]) + ");";
                    
                    lines.push_back({spacer, 0, 2, label});
                }
            }
        }
        
        if (integral[2] == 0)
        {
            // scaled integrals (2 alpha)
            
            for (const auto& tint : bra_rec_base_integrals)
            {
                if ((tint[0] == 0) && (tint[1] == integral[1]) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
                {
                    std::string label = "t4cfunc::ket_transform<" + std::to_string(tint[2]) + ", " + std::to_string(tint[3]) + ">";
                        
                    label += "(skbuffer, "  +  std::to_string(_get_half_spher_index(0, tint, skints)) + ", ";
                    
                    label += "cbuffer, " + std::to_string(_get_index(0, tint, cints))  + ", ";
                    
                    label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]) + ");";
                        
                    lines.push_back({spacer, 0, 2, label});
                }
            }
            
            // scaled integrals (4 alpha)
            
            const auto cstart = _get_geom20_cart_2a_size(bra_rec_base_integrals, ket_rec_base_integrals, integral);
            
            const auto skstart = _get_geom20_half_spher_2a_size(bra_rec_base_integrals, ket_rec_base_integrals, integral);
            
            for (const auto& tint : bra_rec_base_integrals)
            {
                if ((tint[0] == 0) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
                {
                    std::string label = "t4cfunc::ket_transform<" + std::to_string(tint[2]) + ", " + std::to_string(tint[3]) + ">";
                        
                    label += "(skbuffer, "  +  std::to_string(_get_half_spher_index(skstart, tint, skints)) + ", ";
                    
                    label += "cbuffer, " + std::to_string(_get_index(cstart, tint, cints))  + ", ";
                    
                    label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]) + ");";
                        
                    lines.push_back({spacer, 0, 2, label});
                }
            }
        }
    }
}

void
T4CGeomFuncBodyDriver::_add_bra_hrr_call_tree(      VCodeLines&  lines,
                                              const SI4CIntegrals& bra_base_integrals,
                                              const SI4CIntegrals& bra_rec_base_integrals,
                                              const SI4CIntegrals& ket_base_integrals,
                                              const SI4CIntegrals& ket_rec_base_integrals,
                                              const I4CIntegral&   integral,
                                              const size_t         spacer) const
{
    const auto geom_orders = integral.prefixes_order();
    
    if (geom_orders == std::vector<int>({1, 0, 0, 0}))
    {
        auto skints = _get_half_spher_buffers_integrals(bra_base_integrals, ket_base_integrals, integral);
        
        for (const auto& tint : bra_base_integrals)
        {
            if (tint[0] > 0)
            {
                const auto name = t4c::bra_hrr_compute_func_name(tint);
                
                auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
                
                label += _get_bra_hrr_arguments(0, tint, skints);
                
                label += "r_ab, ";
                
                label += std::to_string(tint[2]) + ", " + std::to_string(tint[3]);
                
                label += ");";
                
                lines.push_back({spacer, 0, 2, label});
            }
        }
        
        const auto skcomps = _get_all_half_spher_components(skints);
        
        skints = _get_half_spher_buffers_integrals(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        for (const auto& tint : bra_rec_base_integrals)
        {
            if (tint[0] > 0)
            {
                const auto name = t4c::bra_hrr_compute_func_name(tint);
                
                auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
                
                label += _get_bra_hrr_arguments(skcomps, tint, skints);
                
                label += "r_ab, ";
                
                label += std::to_string(tint[2]) + ", " + std::to_string(tint[3]);
                
                label += ");";
                
                lines.push_back({spacer, 0, 2, label});
            }
        }
    }
    
    if (geom_orders == std::vector<int>({2, 0, 0, 0}))
    {
        auto skints = _get_half_spher_buffers_integrals(bra_base_integrals, ket_base_integrals, integral);
        
        
        const auto skstart = _get_geom20_half_spher_2a_size(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        skints = _get_half_spher_buffers_integrals(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        for (const auto& tint : bra_rec_base_integrals)
        {
            if (tint[0] > 0)
            {
                const auto name = t4c::bra_hrr_compute_func_name(tint);
                
                auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
                
                label += _get_bra_hrr_arguments(skstart, tint, skints);
                
                label += "r_ab, ";
                
                label += std::to_string(tint[2]) + ", " + std::to_string(tint[3]);
                
                label += ");";
                
                lines.push_back({spacer, 0, 2, label});
            }
        }
    }
    
}

void
T4CGeomFuncBodyDriver::_add_bra_geom_hrr_call_tree(      VCodeLines&  lines,
                                                   const SI4CIntegrals& geom_integrals,
                                                   const SI4CIntegrals& bra_base_integrals,
                                                   const SI4CIntegrals& bra_rec_base_integrals,
                                                   const SI4CIntegrals& ket_base_integrals,
                                                   const SI4CIntegrals& ket_rec_base_integrals,
                                                   const I4CIntegral&   integral,
                                                   const size_t         spacer) const
{
    const auto geom_orders = integral.prefixes_order();
    
    if (geom_orders == std::vector<int>({1, 0, 0, 0}))
    {
        auto bskints = _get_half_spher_buffers_integrals(bra_base_integrals, ket_base_integrals, integral);
        
        auto btcomps = _get_all_half_spher_components(bskints);
        
        auto rskints = _get_half_spher_buffers_integrals(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        auto rtcomps = _get_all_half_spher_components(rskints);
        
        const auto gints = _get_geom_half_spher_buffers_integrals(geom_integrals, integral);
        
        for (const auto& tint : gints)
        {
            if ((tint[0] > 0) && (!tint.prefixes().empty()))
            {
                const auto name = t4c::bra_geom_hrr_compute_func_name(tint);
                
                auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
                
                label += std::to_string(_get_geom_half_spher_index(btcomps + rtcomps, tint, gints)) + ", ";
                
                for (const auto& cint : t4c::get_bra_geom_hrr_integrals(tint))
                {
                    if (const auto gorders = cint.prefixes_order(); gorders.empty())
                    {
                        label += std::to_string(_get_geom_half_spher_index(0, cint, bskints)) + ", ";
                    }
                    else
                    {
                        if (cint[0] == 0)
                        {
                            const auto rint = *cint.shift(1, 0);
                            
                            label += std::to_string(_get_geom_half_spher_index(btcomps, rint.base(), rskints)) + ", ";
                        }
                        else
                        {
                            label += std::to_string(_get_geom_half_spher_index(btcomps + rtcomps, cint, gints)) + ", ";
                        }
                    }
                }
                
                label += "r_ab, ";
                
                label += std::to_string(tint[2]) + ", " + std::to_string(tint[3]);
                
                label += ");";
                
                lines.push_back({spacer, 0, 2, label});
            }
        }
    }
    
    if (geom_orders == std::vector<int>({2, 0, 0, 0}))
    {
        auto rskints = _get_half_spher_buffers_integrals(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        auto rscomps = _get_all_half_spher_components(rskints);
        
        auto r2acomps = _get_geom20_half_spher_2a_size(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        if (integral[0] == 0)
        {
            auto name = t4c::bra_geom_hrr_compute_func_name(integral);
            
            auto label = t4c::namespace_label(integral) + "::" + name + "(skbuffer, ";
            
            label += std::to_string(rscomps + r2acomps) + ", ";
            
            label += std::to_string(_get_half_spher_index(0, integral.base(), rskints)) + ", ";
            
            const auto rint = *integral.base().shift(2, 0);
            
            label += std::to_string(_get_half_spher_index(r2acomps, rint, bra_rec_base_integrals)) + ", ";
            
            label += std::to_string(integral[1]) + ", "  + std::to_string(integral[2]) + ", " + std::to_string(integral[3]);
            
            label += ");";
            
            lines.push_back({spacer, 0, 2, label});
        }
    }
}

void
T4CGeomFuncBodyDriver::_add_bra_trafo_call_tree(      VCodeLines&  lines,
                                                const SI4CIntegrals& geom_integrals,
                                                const SI4CIntegrals& bra_base_integrals,
                                                const SI4CIntegrals& bra_rec_base_integrals,
                                                const SI4CIntegrals& ket_base_integrals,
                                                const SI4CIntegrals& ket_rec_base_integrals,
                                                const I4CIntegral&   integral) const
{
    const auto geom_orders = integral.prefixes_order();
    
    if (geom_orders == std::vector<int>{1, 0, 0, 0})
    {
        auto skints = _get_half_spher_buffers_integrals(bra_base_integrals, ket_base_integrals, integral);
        
        auto btcomps = _get_all_half_spher_components(skints);
        
        skints = _get_half_spher_buffers_integrals(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        auto rtcomps = _get_all_half_spher_components(skints);
        
        if (integral[0] > 0)
        {
            skints = _get_geom_half_spher_buffers_integrals(geom_integrals, integral);
        }
        
        if (integral.prefixes_order() == std::vector<int>({1, 0, 0, 0}))
        {
            auto angpair = std::array<int, 2>({integral[0], integral[1]});
                            
            auto bcart_comps = t2c::number_of_cartesian_components(angpair);
            
            auto bspher_comps = t2c::number_of_spherical_components(angpair);
                            
            angpair = std::array<int, 2>({integral[2], integral[3]});
                            
            auto kspher_comps = t2c::number_of_spherical_components(angpair);
            
            auto gindex = _get_geom_half_spher_index(btcomps + rtcomps, integral, skints);
            
            if (integral[0] == 0)
            {
                const auto tint = *integral.shift(1, 0);
                
                gindex = _get_geom_half_spher_index(btcomps, tint.base(), skints);
            }
            
            for (int i = 0; i < 3; i++)
            {
                std::string label = "t4cfunc::bra_transform<" + std::to_string(integral[0]) + ", " + std::to_string(integral[1]) + ">";
                    
                label += "(sbuffer, " + std::to_string(i * bspher_comps * kspher_comps) + ", skbuffer, ";
                
                label += std::to_string(gindex + i * bcart_comps * kspher_comps) + ", ";
                
                label += std::to_string(integral[2]) + ", " + std::to_string(integral[3]) + ");";
                
                lines.push_back({3, 0, 2, label});
            }
        }
    }
    
    if (geom_orders == std::vector<int>{2, 0, 0, 0})
    {
        auto angpair = std::array<int, 2>({integral[0], integral[1]});
                        
        auto bcart_comps = t2c::number_of_cartesian_components(angpair);
        
        auto bspher_comps = t2c::number_of_spherical_components(angpair);
                        
        angpair = std::array<int, 2>({integral[2], integral[3]});
                        
        auto kspher_comps = t2c::number_of_spherical_components(angpair);
        
        auto rskints = _get_half_spher_buffers_integrals(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        auto rscomps = _get_all_half_spher_components(rskints);
        
        auto r2acomps = _get_geom20_half_spher_2a_size(bra_rec_base_integrals, ket_rec_base_integrals, integral);
        
        auto gindex = rscomps + r2acomps; 
        
        for (int i = 0; i < 6; i++)
        {
            std::string label = "t4cfunc::bra_transform<" + std::to_string(integral[0]) + ", " + std::to_string(integral[1]) + ">";
                
            label += "(sbuffer, " + std::to_string(i * bspher_comps * kspher_comps) + ", skbuffer, ";
            
            label += std::to_string(gindex + i * bcart_comps * kspher_comps) + ", ";
            
            label += std::to_string(integral[2]) + ", " + std::to_string(integral[3]) + ");";
            
            lines.push_back({3, 0, 2, label});
        }
    }
    
    std::string label = "distributor.distribute(sbuffer, 0, a_indices, b_indices, c_indices, d_indices, ";
    
    label += std::to_string(integral[0]) + ", ";
    
    label += std::to_string(integral[1]) + ", ";
    
    label += std::to_string(integral[2]) + ", ";
    
    label += std::to_string(integral[3]) + ", ";
   
    label += "j, ket_range);";
   
    lines.push_back({3, 0, 1, label});
}

std::string
T4CGeomFuncBodyDriver::_get_vrr_arguments(const size_t start,
                                          const SI4CIntegrals& integrals,
                                          const I4CIntegral&  integral) const
{
    std::string label = std::to_string(_get_index(start, integral, integrals)) + ", ";
    
    for (const auto& tint : t4c::get_vrr_integrals(integral))
    {
        label += std::to_string(_get_index(start, tint, integrals)) + ", ";
    }
    
    return label;
}

std::string
T4CGeomFuncBodyDriver::_get_ket_hrr_arguments(const size_t       start,
                                              const I4CIntegral& integral,
                                              const SI4CIntegrals& bra_integrals,
                                              const SI4CIntegrals& ket_integrals) const
{
    std::string label;
    
    if (integral[2] == 1)
    {
        const auto cints = _get_cart_buffer_integrals(bra_integrals, ket_integrals);
        
        for (const auto& tint : t4c::get_ket_hrr_integrals(integral))
        {
            label += std::to_string(_get_index(start, tint, cints)) + ", ";
        }
    }
    else
    {
        const auto ckints = _get_contr_buffers_integrals(ket_integrals);
        
        for (const auto& tint : t4c::get_ket_hrr_integrals(integral))
        {
            label += std::to_string(_get_index(start, tint, ckints)) + ", ";
        }
    }
    
    return label;
}

std::string
T4CGeomFuncBodyDriver::_get_bra_hrr_arguments(const size_t start,
                                              const I4CIntegral& integral,
                                              const SI4CIntegrals& integrals) const
{
    std::string label = std::to_string(_get_half_spher_index(start, integral, integrals)) + ", ";
    
    for (const auto& tint : t4c::get_bra_hrr_integrals(integral))
    {
        label += std::to_string(_get_half_spher_index(start, tint, integrals))  + ", ";
    }
    
    return label;
}

SI4CIntegrals
T4CGeomFuncBodyDriver::_get_half_spher_buffers_integrals(const SI4CIntegrals& bra_integrals,
                                                         const SI4CIntegrals& ket_integrals,
                                                         const I4CIntegral&   integral) const
{
    SI4CIntegrals tints;
    
    std::vector<std::string> vstr;
    
    for (const auto& tint : ket_integrals)
    {
        if ((tint[2] == integral[2]) && (tint[3] == integral[3]))
        {
            tints.insert(tint);
        }
    }
    
    for (const auto& tint : bra_integrals)
    {
        if ((tint[0] >= 0) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
        {
            tints.insert(tint);
        }
    }
    
    return tints;
}

SI4CIntegrals
T4CGeomFuncBodyDriver::_get_geom_half_spher_buffers_integrals(const SI4CIntegrals& integrals,
                                                              const I4CIntegral&   integral) const
{
    SI4CIntegrals tints;
    
    std::vector<std::string> vstr;
    
    if (integral[0] > 0)
    {
        for (const auto& tint : integrals)
        {
            if ((tint[0] > 0) && (tint[2] == integral[2]) && (tint[3] == integral[3]) && (!tint.prefixes().empty()))
            {
                tints.insert(tint);
            }
        }
    }
    
    tints.insert(integral);
    
    return tints;
}


size_t
T4CGeomFuncBodyDriver::_get_all_components(const SI4CIntegrals& integrals) const
{
    size_t tcomps= 0;
    
    for (const auto& tint : integrals)
    {
        tcomps += tint.components<T2CPair, T2CPair>().size();
    }
    
    return tcomps;
}

size_t
T4CGeomFuncBodyDriver::_get_all_half_spher_components(const SI4CIntegrals& integrals) const
{
    size_t tcomps = 0;
    
    for (const auto& tint : integrals)
    {
        auto angpair = std::array<int, 2>({tint[2], tint[3]});
                
        auto icomps = t2c::number_of_spherical_components(angpair);
            
        angpair = std::array<int, 2>({tint[0], tint[1]});
                
        icomps *= t2c::number_of_cartesian_components(angpair);
        
        tcomps += icomps;
    }
    
    return tcomps;
}

size_t
T4CGeomFuncBodyDriver::_get_all_geom_half_spher_components(const SI4CIntegrals& integrals) const
{
    size_t tcomps = 0;
    
    for (const auto& tint : integrals)
    {
        auto angpair = std::array<int, 2>({tint[2], tint[3]});
                
        auto icomps = t2c::number_of_spherical_components(angpair);
            
        angpair = std::array<int, 2>({tint[0], tint[1]});
                
        icomps *= t2c::number_of_cartesian_components(angpair);
        
        for (const auto& prefix : tint.prefixes())
        {
            icomps *= prefix.components().size();
        }
        
        tcomps += icomps;
    }
    
    return tcomps;
}

size_t
T4CGeomFuncBodyDriver::_get_all_spher_components(const I4CIntegral& integral) const
{
    auto angpair = std::array<int, 2>({integral[2], integral[3]});
                    
    auto tcomps = t2c::number_of_spherical_components(angpair);
                    
    angpair = std::array<int, 2>({integral[0], integral[1]});
                    
    tcomps *= t2c::number_of_spherical_components(angpair);
    
    for (const auto& prefix : integral.prefixes())
    {
        tcomps *= prefix.components().size();
    }
    
    return tcomps;
}

size_t
T4CGeomFuncBodyDriver::_get_geom20_cart_2a_size(const SI4CIntegrals& bra_integrals,
                                                const SI4CIntegrals& ket_integrals,
                                                const I4CIntegral& integral) const
{
    size_t tcomp = 0;
    
    for (const auto& tint : _get_cart_buffer_integrals(bra_integrals, ket_integrals))
    {
        if (((tint[0] + tint[2]) == 0) && (tint[1] == integral[1]))
        {
            tcomp += tint.components<T2CPair, T2CPair>().size();
        }
    }
    
    return tcomp;
}

size_t
T4CGeomFuncBodyDriver::_get_geom20_contr_2a_size(const SI4CIntegrals& integrals,
                                                 const I4CIntegral& integral) const
{
    size_t tcomps = 0;
    
    for (const auto& tint : _get_contr_buffers_integrals(integrals))
    {
        if ((tint[0] == 0) && (tint[1] == integral[1]) && (tint[2] > 0))
        {
            tcomps += tint.components<T2CPair, T2CPair>().size();
        }
    }
    
    return tcomps;
}

size_t
T4CGeomFuncBodyDriver::_get_geom20_half_spher_2a_size(const SI4CIntegrals& bra_integrals,
                                                      const SI4CIntegrals& ket_integrals,
                                                      const I4CIntegral& integral) const
{
    size_t tcomps = 0;
    
    for (const auto& tint : _get_half_spher_buffers_integrals(bra_integrals, ket_integrals, integral))
    {
        if ((tint[0] == 0) && (tint[1] == integral[1]) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
        {
            auto angpair = std::array<int, 2>({tint[2], tint[3]});
                    
            auto icomps = t2c::number_of_spherical_components(angpair);
                
            angpair = std::array<int, 2>({tint[0], tint[1]});
                    
            icomps *= t2c::number_of_cartesian_components(angpair);
            
            tcomps += icomps;
        }
    }
    
    return tcomps;
}

size_t
T4CGeomFuncBodyDriver::_get_geom20_half_spher_size(const I4CIntegral& integral) const
{
    size_t tcomps = 0;
    
    auto angpair = std::array<int, 2>({integral[2], integral[3]});
            
    tcomps += t2c::number_of_spherical_components(angpair);
        
    angpair = std::array<int, 2>({integral[0], integral[1]});
            
    tcomps *= 6 * t2c::number_of_cartesian_components(angpair);
    
    return tcomps;
}
