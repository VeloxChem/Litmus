#include "t3c_geom_cpu_generators.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "file_stream.hpp"

#include "v3i_geom100_eri_driver.hpp"
#include "v3i_eri_driver.hpp"

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
                    
                    _add_ket_hrr_terms_group(geom_terms);
                        
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
        V3IGeom100ElectronRepulsionDriver geom_drv;
        
        tints = geom_drv.apply_bra_hrr_recursion(integral);
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
