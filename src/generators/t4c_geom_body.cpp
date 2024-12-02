#include "t4c_geom_body.hpp"

#include <algorithm>

#include "t4c_utils.hpp"
#include "t2c_utils.hpp"
#include "t4c_vrr_eri_driver.hpp"
#include "t4c_center_driver.hpp"

void
T4CGeomFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                       const SG4Terms&      cterms,
                                       const SG4Terms&      ckterms,
                                       const SG4Terms&      skterms,
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
    
    for (const auto& label : _get_cart_buffers_def(cterms, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
   
    for (const auto& label : _get_contr_buffers_def(ckterms, integral))
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

    _add_ket_hrr_call_tree(lines, cterms, ckterms, integral, 3);

    _add_ket_trafo_call_tree(lines, cterms, ckterms, skterms, integral, 3);

    _add_bra_hrr_call_tree(lines, skterms, integral, 3);
    
    _add_bra_geom_hrr_call_tree(lines, skterms, integral, 3);
    
    _add_bra_trafo_call_tree(lines, skterms, integral);
    
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
T4CGeomFuncBodyDriver::_get_index(const G4Term&   term,
                                  const SG4Terms& terms) const
{
    size_t index = 0;
    
    for (const auto& cterm : terms)
    {
        if (term == cterm) return index;
        
        index += cterm.second.components<T2CPair, T2CPair>().size();
    }
    
    return 0;
}

bool
T4CGeomFuncBodyDriver::_find_term(const G4Term&   term,
                                  const SG4Terms& terms) const
{
    for (const auto& cterm : terms)
    {
        if (term == cterm) return true;
    }
    
    return false;
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
T4CGeomFuncBodyDriver::_get_half_spher_index(const G4Term&   term,
                                             const SG4Terms& terms) const
{
    size_t index = 0;
    
    for (const auto& cterm : terms)
    {
        if (term == cterm) return index;
        
        const auto tint = cterm.second;
        
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
    
    return 0;
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
T4CGeomFuncBodyDriver::_get_cart_buffers_def(const SG4Terms&    cterms,
                                             const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned Cartesian integrals");
    
    size_t tcomps= 0;
    
    for (const auto& term : cterms)
    {
        tcomps += term.second.components<T2CPair, T2CPair>().size();
    }
    
    std::string label = "CSimdArray<double> cbuffer";
    
    label += "(" + std::to_string(tcomps) + ", 1);";
    
    vstr.push_back(label);
    
    return vstr;
}

std::vector<std::string>
T4CGeomFuncBodyDriver::_get_contr_buffers_def(const SG4Terms&    ckterms,
                                              const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    size_t tcomps= 0;
    
    for (const auto& term : ckterms)
    {
        tcomps += term.second.components<T2CPair, T2CPair>().size();
    }
    
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
T4CGeomFuncBodyDriver::_get_half_spher_buffers_def(const SG4Terms&    skterms,
                                                   const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned half transformed integrals");
    
    size_t tcomps = 0;
    
    for (const auto& term : skterms)
    {
        const auto tint = term.second;
        
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
                    
    auto angpair = std::array<int, 2>({integral[2], integral[3]});
                    
    auto tcomps = t2c::number_of_spherical_components(angpair);
                    
    angpair = std::array<int, 2>({integral[0], integral[1]});
                    
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
T4CGeomFuncBodyDriver::_add_ket_loop_end(      VCodeLines&    lines,
                                         const SG4Terms&      cterms,
                                         const SI4CIntegrals& vrr_integrals,
                                         const I4CIntegral&   integral) const
{
    // non-scaled integrals
    
    for (const auto& term : cterms)
    {
        if (term.first == std::array<int, 4>({0, 0, 0, 0}))
        {
            const auto tint = term.second;
            
            std::string label = "t2cfunc::reduce(cbuffer, ";
                
            label +=  std::to_string(_get_index(term, cterms)) + ", ";
                
            label += "pbuffer, ";
                
            label += std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
                
            label += std::to_string(tint.components<T2CPair, T2CPair>().size()) + ", ";
                
            label += "ket_width, ket_npgtos);";
                
            lines.push_back({4, 0, 2, label});
        }
    }
    
    // scalec integrals on center B
    
    for (const auto& term : cterms)
    {
        if (term.first == std::array<int, 4>({0, 1, 0, 0}))
        {
            const auto tint = term.second;
            
            std::string label = "pbuffer.scale(2.0 * b_exp, {";
           
            label +=  std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
           
            label +=  std::to_string(_get_index(0, tint, vrr_integrals) + tint.components<T2CPair, T2CPair>().size()) + "});";
           
            lines.push_back({4, 0, 2, label});
        }
    }
    
    for (const auto& term : cterms)
    {
        if (term.first == std::array<int, 4>({0, 1, 0, 0}))
        {
            const auto tint = term.second;
            
            std::string label = "t2cfunc::reduce(cbuffer, ";
                
            label +=  std::to_string(_get_index(term, cterms)) + ", ";
                
            label += "pbuffer, ";
                
            label += std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
                
            label += std::to_string(tint.components<T2CPair, T2CPair>().size()) + ", ";
                
            label += "ket_width, ket_npgtos);";
                
            lines.push_back({4, 0, 2, label});
        }
    }
    
    // scaled integrals on center B
    
    for (const auto& term : cterms)
    {
        if (term.first == std::array<int, 4>({1, 0, 0, 0}))
        {
            const auto tint = term.second;
            
            std::string label;
            
            const auto gterm = t4c::prune_term(G4Term({std::array<int, 4>({0, 1, 0, 0}), tint}));
                
            if (_find_term(gterm, cterms))
            {
                label += "pbuffer.scale(a_exp / b_exp, {";
            }
            else
            {
                label += "pbuffer.scale(2.0 * a_exp, {";
            }
           
            label +=  std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
           
            label +=  std::to_string(_get_index(0, tint, vrr_integrals) + tint.components<T2CPair, T2CPair>().size()) + "});";
           
            lines.push_back({4, 0, 2, label});
        }
    }
    
    for (const auto& term : cterms)
    {
        if (term.first == std::array<int, 4>({1, 0, 0, 0}))
        {
            const auto tint = term.second;
            
            std::string label = "t2cfunc::reduce(cbuffer, ";
                
            label +=  std::to_string(_get_index(term, cterms)) + ", ";
                
            label += "pbuffer, ";
                
            label += std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
                
            label += std::to_string(tint.components<T2CPair, T2CPair>().size()) + ", ";
                
            label += "ket_width, ket_npgtos);";
                
            lines.push_back({4, 0, 2, label});
        }
    }
    
    // scaled integrals on centers A and B
    
    for (const auto& term : cterms)
    {
        if (term.first == std::array<int, 4>({1, 1, 0, 0}))
        {
            const auto tint = term.second;
            
            std::string label;
            
            const auto gterm = t4c::prune_term(G4Term({std::array<int, 4>({1, 0, 0, 0}), tint}));
                
            if (_find_term(gterm, cterms))
            {
                label += "pbuffer.scale(2.0 * b_exp, {";
            }
            else
            {
                label += "pbuffer.scale(4.0 * a_exp * b_exp, {";
            }
                
            label +=  std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
           
            label +=  std::to_string(_get_index(0, tint, vrr_integrals) + tint.components<T2CPair, T2CPair>().size()) + "});";
           
            lines.push_back({4, 0, 2, label});
        }
    }
    
    for (const auto& term : cterms)
    {
        if (term.first == std::array<int, 4>({1, 1, 0, 0}))
        {
            const auto tint = term.second;
            
            std::string label = "t2cfunc::reduce(cbuffer, ";
            
            label +=  std::to_string(_get_index(term, cterms)) + ", ";
            
            label += "pbuffer, ";
            
            label += std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
            
            label += std::to_string(tint.components<T2CPair, T2CPair>().size()) + ", ";
            
            label += "ket_width, ket_npgtos);";
            
            lines.push_back({4, 0, 2, label});
        }
    }
    
    for (const auto& term : cterms)
    {
        if (term.first == std::array<int, 4>({1, 0, 1, 0}))
        {
            const auto tint = term.second;
            
            std::string label;
            
            const auto gterm = t4c::prune_term(G4Term({std::array<int, 4>({1, 0, 0, 0}), tint}));
                
            if (_find_term(gterm, cterms))
            {
                label += "pbuffer.scale(pfactors, 0, 2.0, {";
            }
            else
            {
                label += "pbuffer.scale(pfactors, 0, 4.0 * a_exp, {";
            }
           
            label +=  std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
           
            label +=  std::to_string(_get_index(0, tint, vrr_integrals) + tint.components<T2CPair, T2CPair>().size()) + "});";
           
            lines.push_back({4, 0, 2, label});
        }
    }
    
    for (const auto& term : cterms)
    {
        if (term.first == std::array<int, 4>({1, 0, 1, 0}))
        {
            const auto tint = term.second;
            
            std::string label = "t2cfunc::reduce(cbuffer, ";
            
            label +=  std::to_string(_get_index(term, cterms)) + ", ";
            
            label += "pbuffer, ";
            
            label += std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
            
            label += std::to_string(tint.components<T2CPair, T2CPair>().size()) + ", ";
            
            label += "ket_width, ket_npgtos);";
            
            lines.push_back({4, 0, 2, label});
        }
    }
    
    
    for (const auto& term : cterms)
    {
        if (term.first == std::array<int, 4>({2, 0, 0, 0}))
        {
            const auto tint = term.second;
            
            std::string label;
                
            const auto gterm = t4c::prune_term(G4Term({std::array<int, 4>({1, 0, 0, 0}), tint}));
                
            if (_find_term(gterm, cterms))
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
    
    for (const auto& term : cterms)
    {
        if (term.first == std::array<int, 4>({2, 0, 0, 0}))
        {
            const auto tint = term.second;
            
            std::string label = "t2cfunc::reduce(cbuffer, ";
                
            label +=  std::to_string(_get_index(term, cterms)) + ", ";
                
            label += "pbuffer, ";
                
            label += std::to_string(_get_index(0, tint, vrr_integrals)) + ", ";
                
            label += std::to_string(tint.components<T2CPair, T2CPair>().size()) + ", ";
                
            label += "ket_width, ket_npgtos);";
                
            lines.push_back({4, 0, 2, label});
        }
    }
            
    lines.push_back({3, 0, 2, "}"});
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
                                              const SG4Terms&    cterms,
                                              const SG4Terms&    ckterms,
                                              const I4CIntegral& integral,
                                              const size_t       spacer) const
{
    for (const auto& term : ckterms)
    {
        const auto tint = term.second;
        
        if (!tint.prefixes().empty()) continue;
        
        const auto name = t4c::ket_hrr_compute_func_name(tint);
            
        auto label = t4c::namespace_label(tint) + "::" + name + "(ckbuffer, ";
            
        label += std::to_string(_get_index(term, ckterms)) + ", ";
            
        if (tint[2] == 1)
        {
            label += "cbuffer, ";
        }
            
        label += _get_ket_hrr_arguments(term, cterms, ckterms);
            
        label += "cfactors, 6, ";
            
        label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]);
            
        label += ");";
            
        lines.push_back({spacer, 0, 2, label});
    }
    
    for (const auto& term : ckterms)
    {
        const auto tint = term.second;
        
        if (const auto gorders = tint.prefixes_order(); !gorders.empty())
        {
            if ((gorders[2] + gorders[3]) == 0) continue;
            
            const auto name = t4c::ket_geom_hrr_compute_func_name(tint);
                
            auto label = t4c::namespace_label(tint) + "::" + name + "(ckbuffer, ";
                
            label += std::to_string(_get_index(term, ckterms)) + ", ";
                
            if (tint[2] == 0)
            {
                label += "cbuffer, ";
            }
                
            label += _get_ket_geom_hrr_arguments(term, cterms, ckterms);
                
            label += "cfactors, 6, ";
                
            label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]);
                
            label += ");";
                
            lines.push_back({spacer, 0, 2, label});
        }
    }
}

void
T4CGeomFuncBodyDriver::_add_ket_trafo_call_tree(      VCodeLines&  lines,
                                                const SG4Terms&    cterms,
                                                const SG4Terms&    ckterms,
                                                const SG4Terms&    skterms,
                                                const I4CIntegral& integral,
                                                const size_t       spacer) const
{
    for (const auto& term : skterms)
    {
        std::string label;
        
        const auto tint = term.second;
        
        if (tint[0] == 0 && tint.prefixes().empty())
        {
            if (term.first[2] > 0)
            {
                const auto gcomps = Tensor(term.first[2]).components().size();
               
                const auto ccomps = t2c::number_of_cartesian_components(std::array<int, 2>{tint[0], tint[1]});
                
                const auto scomps = t2c::number_of_spherical_components(std::array<int, 2>{tint[0], tint[1]});
                
                for (size_t i = 0; i < gcomps; i++)
                {
                    label = "t4cfunc::ket_transform<" + std::to_string(tint[2] - term.first[2]) + ", " + std::to_string(tint[3]) + ">";
                    
                    label += "(skbuffer, "  + std::to_string(_get_half_spher_index(term, skterms) + i * scomps ) + ", ";
                    
                    label += "ckbuffer, " + std::to_string(_get_index(term, ckterms) + i * ccomps)  + ", ";
                    
                    label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]) + ");";
                    
                    lines.push_back({spacer, 0, 2, label});
                }
            }
            else
            {
                label = "t4cfunc::ket_transform<" + std::to_string(tint[2]) + ", " + std::to_string(tint[3]) + ">";
                    
                label += "(skbuffer, "  + std::to_string(_get_half_spher_index(term, skterms)) + ", ";
                    
                if (tint[2] == 0)
                {
                    label += "cbuffer, " + std::to_string(_get_index(term, cterms))  + ", ";
                }
                else
                {
                    label += "ckbuffer, " + std::to_string(_get_index(term, ckterms))  + ", ";
                }
                    
                label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]) + ");";
                    
                lines.push_back({spacer, 0, 2, label});
            }
        }
    }
    
    for (const auto& term : skterms)
    {
        std::string label;
        
        const auto tint = term.second;
        
        if (const auto gorders = tint.prefixes_order(); !gorders.empty())
        {
            if ((gorders[0] + gorders[1]) > 0) continue;
            
            if (gorders[2] > 0)
            {
                const auto gcomps = Tensor(gorders[2]).components().size();
               
                const auto bcomps = t2c::number_of_cartesian_components(std::array<int, 2>{tint[0], tint[1]});
                
                const auto kccomps = t2c::number_of_cartesian_components(std::array<int, 2>{tint[2], tint[3]});
                
                const auto kscomps = t2c::number_of_spherical_components(std::array<int, 2>{tint[2], tint[3]});
                
                for (size_t i = 0; i < gcomps; i++)
                {
                    label = "t4cfunc::ket_transform<" + std::to_string(tint[2]) + ", " + std::to_string(tint[3]) + ">";
                    
                    label += "(skbuffer, "  + std::to_string(_get_half_spher_index(term, skterms) + i * bcomps * kscomps ) + ", ";
                    
                    label += "ckbuffer, " + std::to_string(_get_index(term, ckterms) + i * bcomps * kccomps)  + ", ";
                    
                    label += std::to_string(tint[0]) + ", " + std::to_string(tint[1]) + ");";
                    
                    lines.push_back({spacer, 0, 2, label});
                }
            }
        }
    }
}

void
T4CGeomFuncBodyDriver::_add_bra_hrr_call_tree(      VCodeLines&  lines,
                                              const SG4Terms&    skterms,
                                              const I4CIntegral& integral,
                                              const size_t       spacer) const
{
    const auto geom_orders = integral.prefixes_order();
    
    for (const auto& term : skterms)
    {
        const auto tint = term.second;
        
        if ((tint[0] > 0) && tint.prefixes().empty())
        {
            if (term.first[2] > 0)
            {
                const auto gcomps = Tensor(term.first[2]).components().size();
                
                for (size_t i = 0; i < gcomps; i++)
                {
                    const auto name = t4c::bra_hrr_compute_func_name(tint);
                    
                    auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
                    
                    label += _get_bra_hrr_arguments(term, skterms);
                    
                    label += "r_ab, ";
                    
                    label += std::to_string(tint[2] - term.first[2]) + ", " + std::to_string(tint[3]);
                    
                    label += ");";
                    
                    lines.push_back({spacer, 0, 2, label});
                }
            }
            else
            {
                const auto name = t4c::bra_hrr_compute_func_name(tint);
                
                auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
                
                label += _get_bra_hrr_arguments(term, skterms);
                
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
                                                   const SG4Terms&    skterms,
                                                   const I4CIntegral& integral,
                                                   const size_t       spacer) const
{
    for (const auto& term : skterms)
    {
        const auto tint = term.second;
        
        if (tint.prefixes_order() == std::vector<int>({1, 0, 0, 0}))
        {
            const auto name = t4c::bra_geom_hrr_compute_func_name(tint);
           
            auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
            
            label += _get_bra_geom_hrr_arguments(term, skterms);
            
            label += "r_ab, ";
            
            label += std::to_string(tint[2]) + ", " + std::to_string(tint[3]);
            
            label += ");";
            
            lines.push_back({spacer, 0, 2, label});
        }
    }
    
    for (const auto& term : skterms)
    {
        const auto tint = term.second;
        
        if (tint.prefixes_order() == std::vector<int>({0, 1, 0, 0}))
        {
            const auto name = t4c::bra_geom_hrr_compute_func_name(tint);
           
            auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
            
            label += _get_bra_geom_hrr_arguments(term, skterms);
            
            if (tint[0] > 0) label += "r_ab, ";
            
            label += std::to_string(tint[2]) + ", " + std::to_string(tint[3]);
            
            label += ");";
            
            lines.push_back({spacer, 0, 2, label});
        }
    }
    
    for (const auto& term : skterms)
    {
        const auto tint = term.second;
        
        if (tint.prefixes_order() == std::vector<int>({1, 1, 0, 0}))
        {
            const auto name = t4c::bra_geom_hrr_compute_func_name(tint);
           
            auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
            
            label += _get_bra_geom_hrr_arguments(term, skterms);
            
            label += "r_ab, ";
            
            label += std::to_string(tint[2]) + ", " + std::to_string(tint[3]);
            
            label += ");";
            
            lines.push_back({spacer, 0, 2, label});
        }
    }
    
    for (const auto& term : skterms)
    {
        const auto tint = term.second;
        
        if (tint.prefixes_order() == std::vector<int>({2, 0, 0, 0}))
        {
            const auto name = t4c::bra_geom_hrr_compute_func_name(tint);
           
            auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
            
            label += _get_bra_geom_hrr_arguments(term, skterms);
            
            label += "r_ab, ";
            
            if (tint[0] == 0)
            {
                label += std::to_string(tint[1]) + ", ";
            }
            
            label += std::to_string(tint[2]) + ", " + std::to_string(tint[3]);
            
            label += ");";
            
            lines.push_back({spacer, 0, 2, label});
        }
    }
            
    for (const auto& term : skterms)
    {
        const auto tint = term.second;
        
        if (tint.prefixes_order() == std::vector<int>({1, 0, 1, 0}))
        {
            const auto name = t4c::bra_geom_hrr_compute_func_name(tint);
           
            auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
            
            label += _get_bra_geom_hrr_arguments(term, skterms);
            
            label += "r_ab, " ;
            
            label += std::to_string(tint[2]) + ", " + std::to_string(tint[3]);
            
            label += ");";
            
            lines.push_back({spacer, 0, 2, label});
        }
    }
            
//    const auto geom_orders = integral.prefixes_order();
//    
//    if (geom_orders == std::vector<int>({1, 0, 0, 0}))
//    {
//        auto bskints = _get_half_spher_buffers_integrals(bra_base_integrals, ket_base_integrals, integral);
//        
//        auto btcomps = _get_all_half_spher_components(bskints);
//        
//        auto rskints = _get_half_spher_buffers_integrals(bra_rec_base_integrals, ket_rec_base_integrals, integral);
//        
//        auto rtcomps = _get_all_half_spher_components(rskints);
//        
//        const auto gints = _get_geom_half_spher_buffers_integrals(geom_integrals, integral);
//        
//        for (const auto& tint : gints)
//        {
//            if ((tint[0] > 0) && (!tint.prefixes().empty()))
//            {
//                const auto name = t4c::bra_geom_hrr_compute_func_name(tint);
//                
//                auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
//                
//                label += std::to_string(_get_geom_half_spher_index(btcomps + rtcomps, tint, gints)) + ", ";
//                
//                for (const auto& cint : t4c::get_bra_geom_hrr_integrals(tint))
//                {
//                    if (const auto gorders = cint.prefixes_order(); gorders.empty())
//                    {
//                        label += std::to_string(_get_geom_half_spher_index(0, cint, bskints)) + ", ";
//                    }
//                    else
//                    {
//                        if (cint[0] == 0)
//                        {
//                            const auto rint = *cint.shift(1, 0);
//                            
//                            label += std::to_string(_get_geom_half_spher_index(btcomps, rint.base(), rskints)) + ", ";
//                        }
//                        else
//                        {
//                            label += std::to_string(_get_geom_half_spher_index(btcomps + rtcomps, cint, gints)) + ", ";
//                        }
//                    }
//                }
//                
//                label += "r_ab, ";
//                
//                label += std::to_string(tint[2]) + ", " + std::to_string(tint[3]);
//                
//                label += ");";
//                
//                lines.push_back({spacer, 0, 2, label});
//            }
//        }
//    }
//    
//    if (geom_orders == std::vector<int>({2, 0, 0, 0}))
//    {
//        auto rskints = _get_half_spher_buffers_integrals(bra_rec_base_integrals, ket_rec_base_integrals, integral);
//        
//        auto rscomps = _get_all_half_spher_components(rskints);
//        
//        auto r2acomps = _get_geom20_half_spher_2a_size(bra_rec_base_integrals, ket_rec_base_integrals, integral);
//        
//        const auto gints =_get_geom_half_spher_buffers_integrals(geom_integrals, integral); 
//        
//        for (const auto& tint : gints)
//        {
//            if (tint[0] > 0)
//            {
//                if (tint.prefixes_order() == std::vector<int>({1, 0, 0, 0}))
//                {
//                    // FIX ME :
//                }
//            }
//            
//            if (tint[0] == 0)
//            {
//                if (tint.prefixes_order() == std::vector<int>({2, 0, 0, 0}))
//                {
//                    auto name = t4c::bra_geom_hrr_compute_func_name(tint);
//                    
//                    auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
//                    
//                    label += std::to_string(_get_geom_half_spher_index(rscomps + r2acomps, tint, geom_integrals)) + ", ";
//                    
//                    label += std::to_string(_get_half_spher_index(0, tint.base(), rskints)) + ", ";
//                    
//                    const auto rint = *tint.base().shift(2, 0);
//                    
//                    label += std::to_string(_get_half_spher_index(r2acomps, rint, bra_rec_base_integrals)) + ", ";
//                    
//                    label += std::to_string(tint[1]) + ", "  + std::to_string(tint[2]) + ", " + std::to_string(tint[3]);
//                    
//                    label += ");";
//                    
//                    lines.push_back({spacer, 0, 2, label});
//                }
//            }
//            
//            if (tint[0] > 0)
//            {
//                if (tint.prefixes_order() == std::vector<int>({2, 0, 0, 0}))
//                {
//                    auto name = t4c::bra_geom_hrr_compute_func_name(tint);
//                    
//                    auto label = t4c::namespace_label(tint) + "::" + name + "(skbuffer, ";
//                    
//                    label += std::to_string(_get_geom_half_spher_index(rscomps + r2acomps, tint, geom_integrals)) + ", ";
//                    
//                    //label += std::to_string(_get_half_spher_index(0, tint.base(), rskints)) + ", ";
//                    
//                    //const auto rint = *tint.base().shift(2, 0);
//                    
//                    //label += std::to_string(_get_half_spher_index(r2acomps, rint, bra_rec_base_integrals)) + ", ";
//                    
//                    label += std::to_string(tint[1]) + ", "  + std::to_string(tint[2]) + ", " + std::to_string(tint[3]);
//                    
//                    label += ");";
//                    
//                    lines.push_back({spacer, 0, 2, label});
//                }
//            }
//        }
//    }
}

void
T4CGeomFuncBodyDriver::_add_bra_trafo_call_tree(      VCodeLines&  lines,
                                                const SG4Terms&    skterms,
                                                const I4CIntegral& integral) const
{
    size_t gcomps = 1;
    
    for (const auto& prefix : integral.prefixes())
    {
        gcomps *= prefix.components().size();
    }
    
    auto angpair = std::array<int, 2>({integral[0], integral[1]});
    
    auto bccomps = t2c::number_of_cartesian_components(angpair);
    
    auto bscomps = t2c::number_of_spherical_components(angpair);
    
    angpair = std::array<int, 2>({integral[2], integral[3]});
    
    auto kscomps = t2c::number_of_spherical_components(angpair);
    
    auto gterm = t4c::prune_term(G4Term({std::array<int, 4>({0, 0, 0, 0}), integral}));
    
    const auto gindex = _get_half_spher_index(gterm, skterms);
    
    for (size_t i = 0; i < gcomps; i++)
    {
        std::string label = "t4cfunc::bra_transform<" + std::to_string(integral[0]) + ", " + std::to_string(integral[1]) + ">";
       
        label += "(sbuffer, " + std::to_string(i * bscomps * kscomps) + ", skbuffer, ";
       
        label += std::to_string(gindex + i * bccomps * kscomps) + ", ";
       
        label += std::to_string(integral[2]) + ", " + std::to_string(integral[3]) + ");";
       
        lines.push_back({3, 0, 2, label});
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
T4CGeomFuncBodyDriver::_get_ket_hrr_arguments(const G4Term&  term,
                                              const SG4Terms& cterms,
                                              const SG4Terms& ckterms) const
{
    std::string label;
    
    for (const auto& tint : t4c::get_ket_hrr_integrals(term.second))
    {
        if (term.second[2] == 1)
        {
            label += std::to_string(_get_index(G4Term({term.first, tint}), cterms)) + ", ";
        }
        else
        {
            label += std::to_string(_get_index(G4Term({term.first, tint}), ckterms)) + ", ";
        }
    }
    
    return label;
}

std::string
T4CGeomFuncBodyDriver::_get_ket_geom_hrr_arguments(const G4Term&  term,
                                                   const SG4Terms& cterms,
                                                   const SG4Terms& ckterms) const
{
    std::string label;
    
    for (const auto& tint : t4c::get_ket_geom_hrr_integrals(term.second))
    {
        auto efacts = term.first;
        
        if (term.first == std::array<int, 4>({1, 0, 0, 0}))
        {
            efacts = std::array<int, 4>({1, 0, 1, 0});
        }
        
        if (term.second[2] == 0)
        {
            label += std::to_string(_get_index(G4Term({efacts, tint}), cterms)) + ", ";
        }
        else
        {
            label += std::to_string(_get_index(G4Term({efacts, tint}), ckterms)) + ", ";
        }
    }
    
    return label;
}

std::string
T4CGeomFuncBodyDriver::_get_bra_hrr_arguments(const G4Term&  term,
                                              const SG4Terms& skterms) const
{
    std::string label = std::to_string(_get_half_spher_index(term, skterms))  + ", ";
    
    for (const auto& tint : t4c::get_bra_hrr_integrals(term.second))
    {
        label += std::to_string(_get_half_spher_index(G4Term({term.first, tint}), skterms))  + ", ";
    }
    
    return label;
}

std::string
T4CGeomFuncBodyDriver::_get_bra_geom_hrr_arguments(const G4Term&  term,
                                                   const SG4Terms& skterms) const
{
    std::string label = std::to_string(_get_half_spher_index(term, skterms))  + ", ";
    
    const auto tint = term.second;
    
    if (tint[0] == 0)
    {
        if (tint.prefixes_order() == std::vector<int>({0, 1, 0, 0}))
        {
            for (const auto& rtint : t4c::get_aux_geom_hrr_integrals(tint))
            {
                if (rtint[1] > tint[1])
                {
                    auto rterm = G4Term({std::array<int, 4>{1, 1, 0, 0}, rtint});
                    
                    label += std::to_string(_get_half_spher_index(rterm, skterms))  + ", ";
                }
                else
                {
                    auto rterm = G4Term({std::array<int, 4>{1, 0, 0, 0}, rtint});
                    
                    label += std::to_string(_get_half_spher_index(rterm, skterms))  + ", ";
                }
            }
        }
        
        if (tint.prefixes_order() == std::vector<int>({1, 1, 0, 0}))
        {
            for (const auto& rtint : t4c::get_aux_geom_hrr_integrals(tint))
            {
                auto rterm = G4Term({std::array<int, 4>{1, 0, 0, 0}, rtint});
                
                label += std::to_string(_get_half_spher_index(rterm, skterms))  + ", ";
            }
        }
                
        if (tint.prefixes_order() == std::vector<int>({2, 0, 0, 0}))
        {
            auto rterm = G4Term({std::array<int, 4>{1, 0, 0, 0}, tint.base()});
            
            label += std::to_string(_get_half_spher_index(rterm, skterms))  + ", ";
            
            rterm = G4Term({std::array<int, 4>{2, 0, 0, 0}, tint.shift(2, 0)->base()});
                                    
            label += std::to_string(_get_half_spher_index(rterm, skterms))  + ", ";
        }
        
        if (tint.prefixes_order() == std::vector<int>({1, 0, 1, 0}))
        {
            G4Term rterm;
            
            auto cint = tint.shift_prefix(-1, 0, false);
            
            if (tint[2] > 0)
            {
                rterm = G4Term({std::array<int, 4>{0, 0, 0, 0}, *cint});
            }
            else
            {
                rterm = G4Term({std::array<int, 4>{1, 0, 0, 0}, *cint});
            }
            
            label += std::to_string(_get_half_spher_index(rterm, skterms))  + ", ";
            
            if (tint[2] > 0)
            {
                rterm = G4Term({std::array<int, 4>{0, 0, 0, 0}, *(cint->shift(1, 1))});
            }
            else
            {
                rterm = G4Term({std::array<int, 4>{1, 0, 0, 0}, *(cint->shift(1, 1))});
            }
                                    
            label += std::to_string(_get_half_spher_index(rterm, skterms))  + ", ";
        }
    }
    else
    {
        for (const auto& rtint : t4c::get_bra_geom_hrr_integrals(tint))
        {
            auto rterm = t4c::prune_term(G4Term({term.first, rtint}));
            
            label += std::to_string(_get_half_spher_index(rterm, skterms))  + ", ";
        }
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
    
    const auto geom_orders = integral.prefixes_order();
    
    if (geom_orders == std::vector<int>({1, 0, 0, 0}))
    {
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
    }
    
    if (geom_orders == std::vector<int>({2, 0, 0, 0}))
    {
     
        for (const auto& tint : integrals)
        {
            if ((tint[2] == integral[2]) && (tint[3] == integral[3]) && (!tint.prefixes().empty()))
            {
                if ((tint.prefixes_order() == std::vector<int>({1, 0, 0, 0})) && (tint[0] == 0)) continue;
                
                tints.insert(tint);
            }
        }
        
        tints.insert(integral);
    }
    
    
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
        if (((tint[0] + tint[2]) == 0) && (tint[1] >= integral[1]) && (tint[1] <= (integral[0] + integral[1])))
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
        if ((tint[0] == 0) && (tint[1] >= integral[1]) && (tint[1] <= (integral[0] + integral[1])) && (tint[2] > 0))
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
        if ((tint[0] == 0) && (tint[1] >= integral[1]) && (tint[1] <= (integral[0] + integral[1])) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
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
