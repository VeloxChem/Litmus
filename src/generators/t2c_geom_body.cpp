#include "t2c_geom_body.hpp"

#include "t2c_center_driver.hpp"
#include "t2c_utils.hpp"

void
T2CGeomFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                       const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    lines.push_back({1, 0, 2, "const auto nelems = pbuffer.number_of_active_elements();"});
    
    for (const auto& label : _get_factors_str(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    const auto components = integral.components<T1CPair, T1CPair>();
    
    const auto ncomps = static_cast<int>(components.size());
    
    auto bcomps = static_cast<int>(Tensor(integral[0]).components().size());
    
    auto kcomps = static_cast<int>(Tensor(integral[1]).components().size());
    
    std::vector<R2CDist> rec_dists;

    for (const auto& component : components)
    {
        rec_dists.push_back(_get_geom_recursion(component));
    }
    
    const auto prefixes = integral.prefixes();
    
    lines.push_back({1, 0, 1, "for (size_t i = 0; i < op_comps; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    if (prefixes[1].shape().order() == 0)
    {
        lines.push_back({2, 0, 1, "for (size_t j = 0; j < ket_comps; j++)"});
        
        lines.push_back({2, 0, 1, "{"});
    }
    
    auto spacer = (prefixes[1].shape().order() == 0) ? 3 : 2;
    
    for (const auto& label : _get_buffers_str(rec_dists, integral))
    {
        lines.push_back({spacer, 0, 2, label});
    }
    
    if ((integral[0] == 0) && (integral[1] == 0))
    {
        const std::array<int, 2> rec_range({0, ncomps});

        for (const auto& label : _get_buffers_str(integral, components, rec_range))
        {
            lines.push_back({spacer, 0, 2, label});
        }
        
        _add_recursion_loop(lines, integral, components, rec_range);
    }
    else
    {
        if (prefixes[1].shape().order() == 0)
        {
            const auto nblocks = ncomps / bcomps;
            
            for (int i = 0; i < nblocks; i++)
            {
                const std::array<int, 2> rec_range({i * bcomps, (i + 1) * bcomps});
                
                for (const auto& label : _get_buffers_str(integral, components, rec_range))
                {
                    lines.push_back({spacer, 0, 2, label});
                }
                
                _add_recursion_loop(lines, integral, components, {i * bcomps, (i + 1) * bcomps});
                
                if (i < (ncomps - 1))  lines.push_back({0, 0, 1, ""});;
            }
        }
        else
        {
            kcomps = (kcomps == 1) ? bcomps : kcomps; 
            
            auto nblocks = ncomps / kcomps;
            
            for (int i = 0; i < nblocks; i++)
            {
                const std::array<int, 2> rec_range({i * kcomps, (i + 1) * kcomps});
                
                for (const auto& label : _get_buffers_str(integral, components, rec_range))
                {
                    lines.push_back({spacer, 0, 2, label});
                }
                
                _add_recursion_loop(lines, integral, components, {i * kcomps, (i + 1) * kcomps});
                
                if (i < (ncomps - 1))  lines.push_back({0, 0, 1, ""});;
            }
        }
    }
        
    if (prefixes[1].shape().order() == 0)
    {
        lines.push_back({2, 0, 1, "}"});
    }
    
    lines.push_back({1, 0, 2, "}"});
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CGeomFuncBodyDriver::_get_factors_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (integral.prefixes()[1].shape().order() > 0)
    {
        vstr.push_back("// Set up exponents");

        vstr.push_back("auto b_exps = factors.data(0);");
    }
    
    return vstr;
}

R2CDist
T2CGeomFuncBodyDriver::_get_geom_recursion(const T2CIntegral& integral) const
{
    R2CDist rdist{R2CTerm(integral)};

    if (!integral.prefixes().empty())
    {
        T2CCenterDriver geom_drv;

        geom_drv.apply_recursion(rdist);
    }
    
    rdist.simplify();
    
    return rdist;
}

std::vector<std::string>
T2CGeomFuncBodyDriver::_get_buffers_str(const std::vector<R2CDist>& rec_dists,
                                        const I2CIntegral&          integral) const
{
    std::vector<std::string> vstr;
    
    for (const auto& tint : t2c::get_geom_integrals(integral))
    {
        vstr.push_back("// Set up components of auxiliary buffer : " + tint.label());

        const auto tlabel = _get_tensor_label(tint);
        
        int index = 0;
        
        for (const auto& tcomp : tint.components<T1CPair, T1CPair>())
        {
            auto line = "auto " + _get_component_label(tcomp) + " = pbuffer.data(" + t2c::get_index_label(tint);
            
            const auto bcomps = Tensor(tint[0]).components().size();
            
            const auto kcomps = Tensor(tint[1]).components().size();
            
            if (integral.prefixes()[1].shape().order() == 0)
            {
                line += " + i * " + std::to_string(bcomps) + " * ket_comps";
                
                line += " + "  + std::to_string(index) + " * ket_comps + j";
            }
            else
            {
                line += " + i * " + std::to_string(bcomps * kcomps);
                
                line += " + " + std::to_string(index);
            }
                        
            vstr.push_back(line + ");");

            index++;
        }
    }
    
    return vstr;
}

std::vector<std::string>
T2CGeomFuncBodyDriver::_get_buffers_str(const I2CIntegral&        integral,
                                        const VT2CIntegrals&      components,
                                        const std::array<int, 2>& rec_range) const
{
    std::vector<std::string> vstr;
    
    auto label = t2c::get_buffer_label(integral, "prim");
    
    if ((rec_range[1] - rec_range[0]) == static_cast<int>(components.size()))
    {
        vstr.push_back("// Set up components of targeted buffer : " + integral.label());
    }
    else
    {
        vstr.push_back("// Set up " + std::to_string(rec_range[0]) + "-" + std::to_string(rec_range[1]) +
                       " components of targeted buffer : " + integral.label());
    }
    
    const auto bcomps = static_cast<int>(Tensor(integral[0]).components().size());
    
    const auto kcomps = static_cast<int>(Tensor(integral[1]).components().size());
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        auto line = "auto " + _get_component_label(components[i]) + " = pbuffer.data(" + t2c::get_index_label(integral);
       
        if (integral.prefixes()[1].shape().order() == 0)
        {
            line += " + " + std::to_string(i / bcomps) + " * op_comps * " + std::to_string(bcomps) + " * ket_comps";
            
            line += " + i * " + std::to_string(bcomps) + " * ket_comps";
            
            line += " + "  + std::to_string(i % bcomps) + " * ket_comps + j";
        }
        else
        {
            line += " + " + std::to_string(i / (bcomps * kcomps)) + " * op_comps * " + std::to_string(bcomps * kcomps);
            
            line += " + i * " + std::to_string(bcomps * kcomps)+ " + " + std::to_string(i % (bcomps * kcomps));
        }
        
        vstr.push_back(line + ");");
    }
    
    return vstr;
}

std::string
T2CGeomFuncBodyDriver::_get_tensor_label(const I2CIntegral& integral) const
{
    std::string label;
    
    if (integral.integrand().name() == "R") label = "to";

    return label;
}

std::string
T2CGeomFuncBodyDriver::_get_tensor_label(const T2CIntegral& integral) const
{
    std::string label;
    
    if (integral.integrand().name() == "R") label = "to";

    return label;
}

std::string
T2CGeomFuncBodyDriver::_get_component_label(const T2CIntegral& integral) const
{
    std::string label = _get_tensor_label(integral) + "_" + integral.label();
    
    if (integral.integrand().name() == "A")
    {
        label += "_" + std::to_string(integral.order());
    }

    if (integral.integrand().name() == "AG")
    {
        label += "_" + std::to_string(integral.order());
    }
    return label;
}

bool
T2CGeomFuncBodyDriver::_find_integral(const std::vector<R2CDist>& rec_dists,
                                      const T2CIntegral&          integral) const
{
    for (const auto& rdist : rec_dists)
    {
        for (const auto& tint : rdist.unique_integrals())
        {
            if (integral == tint) return true;
        }
    }
    
    return false;
}

void
T2CGeomFuncBodyDriver::_add_recursion_loop(      VCodeLines&         lines,
                                           const I2CIntegral&        integral,
                                           const VT2CIntegrals&      components,
                                           const std::array<int, 2>& rec_range) const
{
    auto spacer = (integral.prefixes()[1].shape().order() == 0) ? 3 : 2;
    
    std::vector<R2CDist> rec_dists;
    
    for (int i = rec_range[0]; i < rec_range[1]; i++)
    {
        rec_dists.push_back(_get_geom_recursion(components[i]));
    }
    
    // set up recursion loop
    
    const auto var_str = _get_pragma_str(integral, rec_dists);
    
    lines.push_back({spacer, 0, 1, "#pragma omp simd aligned(" + var_str + " : 64)"});
    
    lines.push_back({spacer, 0, 1, "for (size_t k = 0; k < nelems; k++)"});
    
    lines.push_back({spacer, 0, 1, "{"});
    
    _get_factor_lines(lines, integral, rec_dists);
    
    for (size_t i = 0; i < rec_dists.size(); i++)
    {

        if (i < (rec_dists.size() - 1))
        {
            lines.push_back({spacer + 1, 0, 2, _get_code_line(rec_dists[i])});
        }
        else
        {
            lines.push_back({spacer + 1, 0, 1, _get_code_line(rec_dists[i])});
        }
    }
    
    lines.push_back({spacer, 0, 1, "}"});
}

std::string
T2CGeomFuncBodyDriver::_get_pragma_str(const I2CIntegral&          integral,
                                       const std::vector<R2CDist>& rec_distributions) const
{
    std::set<std::string> tlabels;
    
    for (const auto& rdist : rec_distributions)
    {
        auto tint = rdist.root().integral();
        
        tlabels.insert(_get_component_label(tint));
        
        for (size_t i = 0; i < rdist.terms(); i++)
        {
            auto tint = rdist[i].integral().base();
            
            tlabels.insert(_get_component_label(tint));
            
            for (const auto& fact : rdist[i].factors())
            {
                if (fact.order() > 0) tlabels.insert(fact.label());
            }
        }
    }
    
    std::string label;
    
    for (const auto& tlabel : tlabels)
    {
        label += tlabel + ", ";
    }
    
    if (integral.prefixes()[1].shape().order() > 0) label += "b_exps";
    
    if (label[label.size() - 2] == ',') label.erase(label.end() - 2);
    
    return label;
}

void
T2CGeomFuncBodyDriver::_get_factor_lines(                VCodeLines& lines,
                                         const I2CIntegral&          integral,
                                         const std::vector<R2CDist>& rec_distributions) const
{
    std::set<std::string> tlabels;
    
    for (const auto& rdist : rec_distributions)
    {
        auto tint = rdist.root().integral();
        
        tlabels.insert(_get_tensor_label(tint) + "_" + tint.label());
        
        for (size_t i = 0; i < rdist.terms(); i++)
        {
            for (const auto& fact : rdist[i].factors())
            {
                if (fact.order() == 0) tlabels.insert(fact.label());
            }
        }
    }
    
    auto spacer = (integral.prefixes()[1].shape().order() == 0) ? 4 : 3;
    
    if (std::find(tlabels.begin(), tlabels.end(), "tbe_0") !=  tlabels.end())
    {
        lines.push_back({spacer, 0, 2, "const double tbe_0 = a_exp;"});
    }
    
    if (std::find(tlabels.begin(), tlabels.end(), "tke_0") !=  tlabels.end())
    {
        lines.push_back({spacer, 0, 2, "const double tke_0 = b_exps[k];"});
    }
}

std::string
T2CGeomFuncBodyDriver::_get_code_line(const R2CDist& rec_distribution) const
{
    auto tint = rec_distribution.root().integral();
    
    std::string line = _get_component_label(tint) + "[k] = ";

    for (size_t i = 0; i < rec_distribution.terms(); i++)
    {
        line += _get_rterm_code(rec_distribution[i], i == 0);
    }

    return line + ";";
}

std::string
T2CGeomFuncBodyDriver::_get_rterm_code(const R2CTerm& rec_term,
                                       const bool     is_first) const
{
    const auto pre_fact = rec_term.prefactor();
        
    auto plabel = pre_fact.label();
        
    if (plabel == "1.0")  plabel = "";
        
    if (plabel == "-1.0") plabel = "-";
    
    if (plabel.size() > 1) plabel += " * ";
    
    auto tint = rec_term.integral().base();
    
    plabel += _get_component_label(tint) + "[k]";

    for (const auto& [fact, nrep] : rec_term.map_of_factors())
    {
        for (int i = 0; i < nrep; i++) plabel += " * " + fact.label();
        
        if (fact.order() > 0) plabel += "[k]";
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
