#include "t4c_geom_body.hpp"

#include <algorithm>

#include "t4c_utils.hpp"
#include "t2c_utils.hpp"
#include "t4c_vrr_eri_driver.hpp"
#include "t4c_center_driver.hpp"

void
T4CGeomFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                       const SI4CIntegrals& geom_integrals,
                                       const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    lines.push_back({1, 0, 2, "const auto ndims = " +
                              t4c::get_geom_buffer_label(integral) +
                              ".number_of_columns();"});
    
    for (const auto& label : _get_buffers_str(geom_integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    // generate integrals group
    
    const auto components = integral.components<T2CPair, T2CPair>();
    
    const auto rgroup = _generate_integral_group(components, integral);
    
    // set up dimensions needed for splitting recursion into subblocks
    
    const auto ncomps = static_cast<int>(components.size());
    
    const auto acomps = t2c::number_of_cartesian_components(integral[0]);
    
    const auto bcomps = t2c::number_of_cartesian_components(integral[1]);
    
    const auto ccomps = t2c::number_of_cartesian_components(integral[2]);
    
    const auto dcomps = t2c::number_of_cartesian_components(integral[3]);
    
    const int ocomps = ncomps / (acomps * bcomps * ccomps * dcomps);
    
    if ((acomps * bcomps * ccomps * dcomps) == 1)
    {
        _add_recursion_loop(lines, rgroup, integral, {0, ncomps});
    }
    else
    {
        int rcomps = dcomps;
        
        if (rcomps == 1) rcomps = ccomps;
        
        if (rcomps == 1) rcomps = bcomps;
        
        if (rcomps == 1) rcomps = acomps;
        
        int blkoff = 0;
        
        for (int i = 0; i < ocomps; i++)
        {
            for (int j = 0; j < ncomps / (ocomps * rcomps); j++)
            {
                _add_recursion_loop(lines, rgroup, integral, {blkoff, blkoff + rcomps});
                
                blkoff += rcomps;
            }
        }
    }
        
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
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
