#include "t3c_geom_cpu_generators.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "file_stream.hpp"

#include "v3i_geom100_eri_driver.hpp"
#include "v3i_eri_driver.hpp"
#include "t3c_utils.hpp"
#include "t3c_geom_docs.hpp"
#include "t3c_geom_decl.hpp"
#include "t3c_geom_body.hpp"

void
T3CGeomCPUGenerator::generate(const std::string&        label,
                              const int                 max_ang_mom,
                              const int                 max_aux_ang_mom,
                              const std::array<int, 3>& geom_drvs) const
{
    if (_is_available(label))
    {
        for (int i = 0; i <= max_aux_ang_mom; i++)
        {
            for (int j = 0; j <= max_ang_mom; j++)
            {
                const auto kstart = ((geom_drvs[1] + geom_drvs[2]) > 0) ? 0 : j;
                
                for (int k = kstart; k <= max_ang_mom; k++)
                {
                    const auto integral = _get_integral(label, {i, j, k}, geom_drvs);
                    
                    const auto geom_integrals = _generate_geom_integral_group(integral);
                    
                    auto geom_terms = _generate_geom_terms_group(geom_integrals, integral);
                    
                    //_add_ket_hrr_terms_group(geom_terms);
                    
                    const auto cterms = _filter_cbuffer_terms(geom_terms);
                    
                    const auto skterms = _filter_skbuffer_terms(integral, geom_terms);
                    
                    const auto vrr_integrals = _generate_vrr_integral_group(geom_terms);
                    
                    _write_cpp_header(cterms, skterms, vrr_integrals, integral);
                        
                    std::cout << " *** REFERENCE: " << integral.prefix_label() << " | " << integral.label() << std::endl;
                    
                    std::cout << " --- GEOM INTEGRALS. --- " << std::endl;
                    
                    for (const auto& tint : geom_integrals)
                    {
                        std::cout << " <>" << tint.prefix_label() << " | " << tint.label() << std::endl;
                    }
                    
                    std::cout << " --- GEOM TERMS. --- " << std::endl;
                
                    for (const auto& term : geom_terms)
                    {
                        std::cout << " * ";
                        
                        for (int t = 0; t < 3; t++) std::cout << term.first[t] << ",";
                        
                        std::cout << " * <>" << term.second.prefix_label() << " | " << term.second.label() << std::endl;
                    }
                    
                    std::cout << " --- CBUFFER TERMS. --- " << std::endl;
            
                    for (const auto& term : cterms)
                    {
                        std::cout << " * ";
                    
                        for (int t = 0; t < 3; t++) std::cout << term.first[t] << ",";
                    
                        std::cout << " * <>" << term.second.prefix_label() << " | " << term.second.label() << std::endl;
                    }
                    
                    std::cout << " --- SKBUFFER TERMS. --- " << std::endl;
    
                    for (const auto& term : skterms)
                    {
                        std::cout << " * ";
            
                        for (int t = 0; t < 3; t++) std::cout << term.first[t] << ",";
            
                        std::cout << " * <>" << term.second.prefix_label() << " | " << term.second.label() << std::endl;
                    }
                    
                    std::cout << " --- VRR INTEGRALS --- " << std::endl;
                
                    for (const auto& tint : vrr_integrals)
                    {
                        std::cout << " <>" << tint.prefix_label() << " | " << tint.label() << "_"  << tint.order() << std::endl;
                    }
                }
            }
        }
    }
    else
    {
        std::cerr << "*** ERROR *** Unsupported type of three-center integral: ";
        
        std::cerr << label << " !!!" << std::endl;
        
        std::exit(EXIT_FAILURE);
    }
}

bool
T3CGeomCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "electron repulsion") return true;
    
    return false;
}

I3CIntegral
T3CGeomCPUGenerator::_get_integral(const std::string&        label,
                                   const std::array<int, 3>& ang_moms,
                                   const std::array<int, 3>& geom_drvs) const
{
    // bra and ket sides

    const auto bpair = I1CPair("GA", ang_moms[0]);

    const auto kpair = I2CPair("GC", ang_moms[1], "GD", ang_moms[2]);

    VOperators prefixes;

    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[0])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[1])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[2])));
    
    // electron repulsion integrals

    if (fstr::lowercase(label) == "electron repulsion")
    {
        return I3CIntegral(bpair, kpair, Operator("1/|r-r'|"), 0, prefixes);
    }
    
    return I3CIntegral();
}

SI3CIntegrals
T3CGeomCPUGenerator::_generate_geom_integral_group(const I3CIntegral& integral) const
{
    const auto geom_order = integral.prefixes_order();
    
    SI3CIntegrals tints;
    
    if (geom_order == std::vector<int>({1, 0, 0}))
    {
        if (integral[1] > 0)
        {
            V3IElectronRepulsionDriver eri_drv;
            
            for (auto& tint : eri_drv.create_ket_hrr_recursion({integral.base(), }))
            {
                auto ctint = tint;
                
                ctint.set_prefixes(integral.prefixes()); 
                
                tints.insert(ctint);
            }
        }
        
        V3IGeom100ElectronRepulsionDriver geom_drv;
    
        auto new_tints = tints;
        
        if (new_tints.empty()) new_tints.insert(integral);
        
        for (const auto& tint : new_tints)
        {
            for (const auto& ctint : geom_drv.apply_bra_hrr_recursion(tint))
            {
                tints.insert(ctint);
            }
        }
    }
  
    tints.insert(integral);
    
    return tints;
}

SG3Terms
T3CGeomCPUGenerator::_generate_geom_terms_group(const SI3CIntegrals& integrals,
                                                const I3CIntegral&   integral) const
{
    SG3Terms terms;
    
    for (const auto& tint : integrals)
    {
        if (integral.prefixes_order() == std::vector<int>({1, 0, 0}))
        {
            if (tint.prefixes_order() == std::vector<int>({1, 0, 0}))
            {
                terms.insert({std::array<int, 3>{0, 0, 0}, tint});
                
                if (tint[0] == 0)
                {
                    terms.insert({std::array<int, 3>{1, 0, 0}, tint.shift(1, 0)->base()});
                }
            }
            else
            {
                if (tint[0] == (integral[0] + 1))
                {
                    terms.insert({std::array<int, 3>{1, 0, 0}, tint});
                }
                else
                {
                    terms.insert({std::array<int, 3>{0, 0, 0}, tint});
                }
            }
        }
    }
    
    return terms;
}

void
T3CGeomCPUGenerator::_add_ket_hrr_terms_group(SG3Terms& terms) const
{
    SG3Terms new_terms;
    
    for (const auto& term : terms)
    {
        if (term.second[1] > 0)
        {
            V3IElectronRepulsionDriver eri_drv;
            
            if (term.second.prefixes().empty())
            {
                for (const auto& tint : eri_drv.create_ket_hrr_recursion({term.second, }))
                {
                    new_terms.insert({term.first, tint});
                }
            }
            else
            {
                new_terms.insert(term);
            }
        }
        else
        {
            new_terms.insert(term);
        }
    }
    
    terms = new_terms;
}

SG3Terms
T3CGeomCPUGenerator::_filter_cbuffer_terms(const SG3Terms& terms) const
{
    SG3Terms new_terms;
    
    for (const auto& term : terms)
    {
        if ((term.second[1] == 0) && term.second.prefixes().empty())
        {
            new_terms.insert(term);
        }
    }
    
    return new_terms;
}

SI3CIntegrals
T3CGeomCPUGenerator::_generate_vrr_integral_group(const SG3Terms& terms) const
{
    SI3CIntegrals tints;
       
    V3IElectronRepulsionDriver eri_drv;
        
    for (const auto& term : terms)
    {
        if ((term.second[1] == 0) && term.second.prefixes().empty())
        {
            const auto ctints = eri_drv.create_vrr_recursion({term.second, });
                
            tints.insert(ctints.cbegin(), ctints.cend());
        }
    }
    
    return tints;
}

SG3Terms
T3CGeomCPUGenerator::_filter_skbuffer_terms(const I3CIntegral& integral,
                                            const SG3Terms& terms) const
{
    SG3Terms new_terms;
    
    const auto gorders = integral.prefixes_order();
    
    for (const auto& term : terms)
    {
        if (gorders == std::vector<int>({1, 0, 0}))
        {
            if (integral[0] == 0)
            {
                if ((term.second[0] == 1) && term.second.prefixes().empty())
                {
                    new_terms.insert(term);
                }
            }
            else
            {
                if ((term.second[0] == integral[0]) && (!term.second.prefixes().empty()))
                {
                    new_terms.insert(term);
                }
            }
        }
    }
    
    return new_terms;
}

void
T3CGeomCPUGenerator::_write_cpp_header(const SG3Terms&      cterms,
                                       const SG3Terms&      skterms,
                                       const SI3CIntegrals& vrr_integrals,
                                       const I3CIntegral& integral) const
{
    auto fname = _file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, true);
    
    _write_hpp_includes(fstream, skterms, vrr_integrals, integral);
    
    _write_namespace(fstream, integral, true);
    
    T3CGeomDocuDriver docs_drv;
    
    T3CGeomDeclDriver decl_drv;

    T3CGeomFuncBodyDriver func_drv;
    
    docs_drv.write_doc_str(fstream, integral);
    
    decl_drv.write_func_decl(fstream, integral, false);
    
    func_drv.write_func_body(fstream, cterms, skterms, vrr_integrals, integral);
    
    fstream << std::endl;

    _write_namespace(fstream, integral, false);
        
    _write_hpp_defines(fstream, integral, false);
    
    fstream.close();
}

std::string
T3CGeomCPUGenerator::_file_name(const I3CIntegral& integral) const
{
    std::string label = "Rec" + integral.label();
    
    return t3c::integral_label(integral) + label;
}

void
T3CGeomCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                        const I3CIntegral&   integral,
                                        const bool           start) const
{
    auto fname = _file_name(integral) + "_hpp";
    
    auto lines = VCodeLines();
 
    if (start)
    {
        lines.push_back({0, 0, 1, "#ifndef " + fname});
        
        lines.push_back({0, 0, 2, "#define " + fname});
    }
    else
    {
        lines.push_back({0, 0, 1, "#endif /* " + fname + " */"});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T3CGeomCPUGenerator::_write_hpp_includes(      std::ofstream& fstream,
                                         const SG3Terms&      skterms,
                                         const SI3CIntegrals& vrr_integrals,
                                         const I3CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "#include <array>"});
    
    lines.push_back({0, 0, 1, "#include <cstddef>"});
        
    lines.push_back({0, 0, 2, "#include <utility>"});
    
    std::set<std::string> labels;
        
    for (const auto& tint : vrr_integrals)
    {
        labels.insert(t3c::prim_file_name(tint));
    }
        
    for (const auto& term : skterms)
    {
        const auto tint = term.second;
        
        if (tint[2] >= integral[2])
        {
            if (tint[1] > 0)
            {
                if (tint.prefixes().empty())
                {
                    labels.insert(t3c::hrr_file_name(tint));
                }
                else
                {
                   // labels.insert(t4c::bra_geom_hrr_file_name(tint));
                }
            }
            else
            {
                // FIX ME !!!
            }
        }
    }
    
    // Add remaining stuff
    
    for (const auto& label : labels)
    {
        lines.push_back({0, 0, 1, "#include \"" + label + ".hpp\""});
    }
    
    //lines.push_back({0, 0, 1, "#include \"" + t4c::geom_file_name(integral) + ".hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"SimdArray.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"BoysFunc.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T4CUtils.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T3CUtils.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T2CUtils.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"BatchFunc.hpp\""});

    lines.push_back({0, 0, 1, "#include \"GtoPairBlock.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"GtoBlock.hpp\""});
   
    ost::write_code_lines(fstream, lines);
}

void
T3CGeomCPUGenerator::_write_namespace(      std::ofstream& fstream,
                                      const I3CIntegral&   integral,
                                      const bool           start) const
{
    const auto label = t3c::namespace_label(integral);
    
    auto lines = VCodeLines();
    
    if (start)
    {
        lines.push_back({0, 0, 2, "namespace " + label + " { // " + label + " namespace"});
    }
    else
    {
        lines.push_back({0, 0, 2, "} // " + label + " namespace"});
    }
    
    ost::write_code_lines(fstream, lines);
}
