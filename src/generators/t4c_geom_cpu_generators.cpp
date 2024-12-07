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

#include "t4c_geom_cpu_generators.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "file_stream.hpp"

#include "t4c_utils.hpp"
#include "t4c_geom_docs.hpp"
#include "t4c_geom_decl.hpp"
#include "t4c_geom_body.hpp"
#include "v4i_center_driver.hpp"
#include "v4i_geom10_eri_driver.hpp"
#include "v4i_geom20_eri_driver.hpp"
#include "v4i_geom11_eri_driver.hpp"
#include "v4i_geom1010_eri_driver.hpp"
#include "v4i_eri_driver.hpp"

void
T4CGeomCPUGenerator::generate(const std::string&        label,
                              const int                 max_ang_mom,
                              const std::array<int, 5>& geom_drvs) const
{
    if (_is_available(label))
    {
        for (int i = 0; i <= max_ang_mom; i++)
        {
            for (int j = 0; j <= max_ang_mom; j++)
            {
                for (int k = 0; k <= max_ang_mom; k++)
                {
                    const auto lstart = ((geom_drvs[3] + geom_drvs[4]) > 0) ? 0 : k;
                    
                    for (int l = lstart; l <= max_ang_mom; l++)
                    {
                        const auto integral = _get_integral(label, {i, j, k, l}, geom_drvs);
                        
                        const auto geom_integrals = _generate_geom_integral_group(integral);
                        
                        auto geom_terms = _generate_geom_terms_group(geom_integrals);
                        
                        _prune_terms_group(geom_terms);
                        
                        _add_bra_hrr_terms_group(geom_terms);
                        
                        _add_ket_hrr_terms_group(geom_terms);
                        
                        const auto cterms = _filter_cbuffer_terms(geom_terms);
                        
                        const auto ckterms = _filter_ckbuffer_terms(geom_terms);
                        
                        const auto skterms = _filter_skbuffer_terms(integral, geom_terms);
                        
                        const auto vrr_integrals = _generate_vrr_integral_group(geom_terms);
                                                                    
                        _write_cpp_header(cterms, ckterms, skterms, vrr_integrals, integral);
                        
//                       if ((i == 2) && (j == 2) && (k == 2) && (l == 2))
//                        {
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
                                
                                for (int t = 0; t < 4; t++) std::cout << term.first[t] << ",";
                                
                                std::cout << " * <>" << term.second.prefix_label() << " | " << term.second.label() << std::endl;
                            }
                        
                            std::cout << " --- CBUFFER TERMS. --- " << std::endl;
                    
                            for (const auto& term : cterms)
                            {
                                std::cout << " * ";
                            
                                for (int t = 0; t < 4; t++) std::cout << term.first[t] << ",";
                            
                                std::cout << " * <>" << term.second.prefix_label() << " | " << term.second.label() << std::endl;
                            }
                        
                            std::cout << " --- CKBUFFER TERMS. --- " << std::endl;
                
                            for (const auto& term : ckterms)
                            {
                                std::cout << " * ";
                        
                                for (int t = 0; t < 4; t++) std::cout << term.first[t] << ",";
                        
                                std::cout << " * <>" << term.second.prefix_label() << " | " << term.second.label() << std::endl;
                            }
                        
                            std::cout << " --- SKBUFFER TERMS. --- " << std::endl;
            
                            for (const auto& term : skterms)
                            {
                                std::cout << " * ";
                    
                                for (int t = 0; t < 4; t++) std::cout << term.first[t] << ",";
                    
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
    }
    else
    {
        std::cerr << "*** ERROR *** Unsupported type of four-center integral: ";
        
        std::cerr << label << " !!!" << std::endl;
        
        std::exit(EXIT_FAILURE);
    }
}

bool
T4CGeomCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "electron repulsion") return true;
    
    return false;
}

I4CIntegral
T4CGeomCPUGenerator::_get_integral(const std::string&        label,
                                   const std::array<int, 4>& ang_moms,
                                   const std::array<int, 5>& geom_drvs) const
{
    // bra and ket sides

    const auto bpair = I2CPair("GA", ang_moms[0], "GB", ang_moms[1]);

    const auto kpair = I2CPair("GC", ang_moms[2], "GD", ang_moms[3]);

    VOperators prefixes;

    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[0])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[1])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[3])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[4])));

    // electron repulsion integrals

    if (fstr::lowercase(label) == "electron repulsion")
    {
        return I4CIntegral(bpair, kpair, Operator("1/|r-r'|"), 0, prefixes);
    }
    
    return I4CIntegral();
}

SI4CIntegrals
T4CGeomCPUGenerator::_generate_geom_integral_group(const I4CIntegral& integral) const
{
    const auto geom_order = integral.prefixes_order();
    
    SI4CIntegrals tints;
    
    if (geom_order == std::vector<int>({1, 0, 0, 0}))
    {
        V4IGeom10ElectronRepulsionDriver geom_drv;
        
        tints = geom_drv.apply_bra_hrr_recursion(integral);
    }
    
    if (geom_order == std::vector<int>({2, 0, 0, 0}))
    {
        V4IGeom20ElectronRepulsionDriver geom_drv;
        
        V4IGeom10ElectronRepulsionDriver grad_drv;
        
        for (const auto& tint : geom_drv.apply_bra_hrr_recursion(integral))
        {
            tints.insert(tint);
            
            if (tint.prefixes_order() == std::vector<int>({1, 0, 0, 0}))
            {
                for (const auto& cint : grad_drv.apply_bra_hrr_recursion(tint))
                {
                    tints.insert(cint);
                }
            }
        }
    }
    
    if (geom_order == std::vector<int>({1, 1, 0, 0}))
    {
        V4IGeom11ElectronRepulsionDriver geom_drv;
        
        V4IGeom10ElectronRepulsionDriver grad_drv;
        
        for (const auto& tint : geom_drv.apply_bra_hrr_recursion(integral))
        {
            tints.insert(tint);
            
            if (tint.prefixes_order() == std::vector<int>({1, 0, 0, 0}))
            {
                for (const auto& cint : grad_drv.apply_bra_hrr_recursion(tint))
                {
                    tints.insert(cint);
                }
            }
        }
    }
    
    if (geom_order == std::vector<int>({1, 0, 1, 0}))
    {
        V4IGeom10ElectronRepulsionDriver grad_drv;
        
        if ((integral[0] == 0) && (integral[2] == 0))
        {
            tints.insert(integral);
            
            return tints; 
        }
        
        if ((integral[0] == 0) && (integral[2] > 0))
        {
            tints.insert(integral);
            
            for (const auto& cint : grad_drv.apply_ket_hrr_recursion(integral))
            {
                tints.insert(cint);
            }
            
            return tints;
        }
        
        if ((integral[0] > 0) && (integral[2] == 0))
        {
            tints.insert(integral);
            
            for (const auto& cint : grad_drv.apply_bra_hrr_recursion(integral))
            {
                tints.insert(cint);
            }
            
            return tints;
        }
        
        tints.insert(integral);
        
        auto cints = grad_drv.apply_bra_hrr_recursion(integral);

        for (auto cint : cints)
        {
            tints.insert(cint);
            
            for (auto rint : grad_drv.apply_ket_hrr_recursion(cint))
            {
                tints.insert(rint);
            }
        }
    }
    
    tints.insert(integral); 
    
    return tints;
}

SG4Terms
T4CGeomCPUGenerator::_generate_geom_terms_group(const SI4CIntegrals& integrals) const
{
    SG4Terms terms;
    
    for (const auto& tint : integrals)
    {
        terms.insert({std::array<int, 4>{0, 0, 0, 0}, tint});
    
        if (tint.prefixes_order() == std::vector<int>({1, 0, 0, 0}))
        {
            if (tint[0] == 0)
            {
                terms.insert({std::array<int, 4>{1, 0, 0, 0}, tint.shift(1, 0)->base()});
            }
        }
        
        if (tint.prefixes_order() == std::vector<int>({0, 0, 1, 0}))
        {
            if (tint[2] == 0)
            {
                terms.insert({std::array<int, 4>{0, 0, 1, 0}, tint.base()});
                
                terms.insert({std::array<int, 4>{0, 0, 1, 0}, tint.shift(1, 3)->base()});
            }
        }
        
        if (tint.prefixes_order() == std::vector<int>({0, 1, 0, 0}))
        {
            if (tint[0] == 0)
            {
                terms.insert({std::array<int, 4>{0, 1, 0, 0}, tint.shift(1, 1)->base()});
                
                if (tint[1] > 0)
                {
                    terms.insert({std::array<int, 4>{0, 0, 0, 0}, tint.shift(-1, 1)->base()});
                }
            }
        }
        
        if (tint.prefixes_order() == std::vector<int>({1, 0, 1, 0}))
        {
            if ((tint[0] == 0) && (tint[2] == 0))
            {
                const auto btint = tint.shift_prefix(-1, 0, false);
                
                const auto bptint = btint->shift(1, 1);
                
                terms.insert({std::array<int, 4>{1, 0, 0, 0}, *btint});
                
                terms.insert({std::array<int, 4>{1, 0, 0, 0}, *bptint});
                
                terms.insert({std::array<int, 4>{1, 0, 1, 0}, btint->base()});
                
                terms.insert({std::array<int, 4>{1, 0, 1, 0}, btint->shift(1, 3)->base()});
                
                terms.insert({std::array<int, 4>{1, 0, 1, 0}, bptint->base()});
                
                terms.insert({std::array<int, 4>{1, 0, 1, 0}, bptint->shift(1, 3)->base()});
            }
            
            if ((tint[0] == 0) && (tint[2] > 0))
            {
                const auto btint = tint.shift_prefix(-1, 0, false);
                
                const auto bptint = btint->shift(1, 1);
                
                terms.insert({std::array<int, 4>{1, 0, 0, 0}, *btint});
                
                terms.insert({std::array<int, 4>{1, 0, 0, 0}, *bptint});
            }
        }
        
        if (tint.prefixes_order() == std::vector<int>({2, 0, 0, 0}))
        {
            if (tint[0] == 0)
            {
                terms.insert({std::array<int, 4>{2, 0, 0, 0}, tint.shift(2, 0)->base()});
                
                terms.insert({std::array<int, 4>{1, 0, 0, 0}, tint.base()});
            }
        }
        
        if (tint.prefixes_order() == std::vector<int>({1, 1, 0, 0}))
        {
            if (tint[0] == 0)
            {
                if (tint[1] == 0)
                {
                    const auto rtint = tint.shift(1, 0);
                    
                    const auto ctint = rtint->shift(1, 1);
                    
                    terms.insert({std::array<int, 4>{1, 1, 0, 0}, ctint->base()});
                }
                else
                {
                    const auto rtint = tint.shift_prefix(-1, 0, false);
                    
                    // up rec. terms
                    
                    if (const auto ctint = rtint->shift(1, 1))
                    {
                        terms.insert({std::array<int, 4>{1, 0, 0, 0}, *ctint});
                        
                        if (const auto ptint = ctint->shift(1, 1))
                        {
                            terms.insert({std::array<int, 4>{1, 1, 0, 0}, ptint->base()});
                        }
                        
                        if (const auto ptint = ctint->shift(-1, 1))
                        {
                            terms.insert({std::array<int, 4>{1, 0, 0, 0}, ptint->base()});
                        }
                    }
                    
                    // down rec. terms
                    
                    terms.insert({std::array<int, 4>{1, 0, 0, 0}, *rtint});
                    
                    terms.insert({std::array<int, 4>{1, 0, 0, 0}, rtint->base()});
                    
                    if (const auto ptint = rtint->shift(1, 1))
                    {
                        terms.insert({std::array<int, 4>{1, 1, 0, 0}, ptint->base()});
                    }
                    
                    if (const auto ptint = rtint->shift(-1, 1))
                    {
                        terms.insert({std::array<int, 4>{1, 0, 0, 0}, ptint->base()});
                    }
                }
            }
        }
    }
    
    return terms;
}

void
T4CGeomCPUGenerator::_add_bra_hrr_terms_group(SG4Terms& terms) const
{
    SG4Terms new_terms;
    
    for (const auto& term : terms)
    {
        if (term.second[0] > 0)
        {
            V4IElectronRepulsionDriver eri_drv;
            
            if (term.second.prefixes().empty())
            {
                for (const auto& tint : eri_drv.create_bra_hrr_recursion({term.second, }))
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

void
T4CGeomCPUGenerator::_add_ket_hrr_terms_group(SG4Terms& terms) const
{
    SG4Terms new_terms;
    
    for (const auto& term : terms)
    {
        if ((term.second[0] == 0) && (term.second[2] > 0))
        {
            V4IElectronRepulsionDriver eri_drv;
            
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

SG4Terms
T4CGeomCPUGenerator::_filter_cbuffer_terms(const SG4Terms& terms) const
{
    SG4Terms new_terms;
    
    for (const auto& term : terms)
    {
        if ((term.second[0] == 0) && (term.second[2] == 0) && term.second.prefixes().empty())
        {
            new_terms.insert(term);
        }
    }
    
    return new_terms;
}

SG4Terms
T4CGeomCPUGenerator::_filter_ckbuffer_terms(const SG4Terms& terms) const
{
    SG4Terms new_terms;
    
    for (const auto& term : terms)
    {
        int korder = 0;
        
        const auto gorders = term.second.prefixes_order();
        
        if (!gorders.empty())
        {
            korder = gorders[0]; 
        }
        
        if ((term.second[0] == 0) && (term.second[2] > 0) && (korder == 0))
        {
            new_terms.insert(term);
        }
        
        if ((term.second[0] == 0) && (gorders == std::vector<int>({0, 0, 1, 0})))
        {
            new_terms.insert(term);
        }
    }
    
    return new_terms;
}

SG4Terms
T4CGeomCPUGenerator::_filter_skbuffer_terms(const I4CIntegral& integral,
                                            const SG4Terms& terms) const
{
    SG4Terms new_terms;
    
    const auto gorders = integral.prefixes_order();
    
    for (const auto& term : terms)
    {
        if ((gorders[2] + gorders[3]) > 0)
        {
            if (const auto corders = term.second.prefixes_order(); !corders.empty())
            {
                if ((corders[2] == gorders[2]) && (corders[3] == gorders[3]) &&
                    (term.second[2] == integral[2]) && (term.second[3] == integral[3]))
                {
                    new_terms.insert(term);
                }
            }
//            else
//            {
//                if (gorders == std::vector<int>({1, 0, 1, 0}))
//                {
//                    if ((integral[2] == (term.second[2] + 1)) &&
//                        (term.second[3] == integral[3])       &&
//                        (term.first == std::array<int, 4>({1, 0, 0, 0})))
//                    {
//                        std::cout << " I am here..." << std::endl;
//                            
//                        new_terms.insert(term);
//                    }
//                }
//            }
        }
        else
        {
            if ((term.second[2] == integral[2]) && (term.second[3] == integral[3]))
            {
                new_terms.insert(term);
            }
        }
    }
    
    return new_terms;
}

SI4CIntegrals
T4CGeomCPUGenerator::_generate_vrr_integral_group(const SG4Terms& terms) const
{
    SI4CIntegrals tints;
       
    V4IElectronRepulsionDriver eri_drv;
        
    for (const auto& term : terms)
    {
        if ((term.second[0] == 0) && (term.second[2] == 0) && term.second.prefixes().empty())
        {
            const auto ctints = eri_drv.create_vrr_recursion({term.second, });
                
            tints.insert(ctints.cbegin(), ctints.cend());
        }
    }
    
    return tints;
}

SI4CIntegrals
T4CGeomCPUGenerator::_generate_geom_base_integral_group(const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    for (const auto& tint : integrals)
    {
        if (tint.prefixes().empty())
        {
            tints.insert(tint);
        }
    }
    
    return tints;
}

SI4CIntegrals
T4CGeomCPUGenerator::_generate_geom_rec_integral_group(const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    for (const auto& tint : integrals)
    {
        if (!tint.prefixes().empty())
        {
            tints.insert(tint);
        }
    }
    
    return tints;
}

SI4CIntegrals
T4CGeomCPUGenerator::_generate_bra_hrr_integral_group(const I4CIntegral&   integral,
                                                      const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    if (integral.prefixes_order() == std::vector<int>{1, 0, 0, 0})
    {
        if (integral[0] == 0) tints.insert((integral.shift(1, 0))->base());
        
        for (const auto& tint : integrals)
        {
            if (tint.prefixes_order() == std::vector<int>{1, 0, 0, 0})
            {
                if (tint[0] == 0) tints.insert((tint.shift(1, 0))->base());
            }
        }
    }
    
    if (integral.prefixes_order() == std::vector<int>{2, 0, 0, 0})
    {
        if (integral[0] == 0)
        {
            tints.insert((integral.shift(2, 0))->base());
            
            tints.insert(integral.base());
        }

        for (const auto& tint : integrals)
        {
            if (tint.prefixes_order() == std::vector<int>{2, 0, 0, 0})
            {
                if (tint[0] == 0)
                {
                    tints.insert((tint.shift(2, 0))->base());
                    
                    tints.insert(tint.base());
                }
            }
            
            if (tint.prefixes_order() == std::vector<int>{1, 0, 0, 0})
            {
                if (tint[0] == 0)
                {
                    tints.insert((tint.shift(2, 0))->base());
                    
                    tints.insert(integral.base());
                }
            }
        }
    }
    
    return tints;
}

SI4CIntegrals
T4CGeomCPUGenerator::_generate_bra_base_hrr_integral_group(const I4CIntegral&   integral,
                                                           const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    if (integral.prefixes_order() == std::vector<int>{1, 0, 0, 0})
    {
        // Electron repulsion integrals
        
        if (integral.integrand() == Operator("1/|r-r'|"))
        {
            V4IElectronRepulsionDriver eri_drv;
            
            tints = eri_drv.create_bra_hrr_recursion(integrals);
        }
    }
    
    if (integral.prefixes_order() == std::vector<int>{2, 0, 0, 0})
    {
        // Electron repulsion integrals
        
        if (integral.integrand() == Operator("1/|r-r'|"))
        {
            V4IElectronRepulsionDriver eri_drv;
            
            tints = eri_drv.create_bra_hrr_recursion(integrals);
        }
    }
    
    
    return tints;
}

SI4CIntegrals
T4CGeomCPUGenerator::_generate_ket_hrr_integral_group(const I4CIntegral&   integral,
                                                      const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    if (integral.prefixes_order() == std::vector<int>{1, 0, 0, 0})
    {
        for (const auto& tint : integrals)
        {
            if (tint[2] > 0) tints.insert(tint);
        }
    }
    
    if (integral.prefixes_order() == std::vector<int>{2, 0, 0, 0})
    {
        for (const auto& tint : integrals)
        {
            if (tint[2] > 0) tints.insert(tint);
        }
    }
    
    return tints;
}

SI4CIntegrals
T4CGeomCPUGenerator::_generate_ket_base_hrr_integral_group(const I4CIntegral&   integral,
                                                           const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    if (integral.prefixes_order() == std::vector<int>{1, 0, 0, 0})
    {
        // Electron repulsion integrals
        
        if (integral.integrand() == Operator("1/|r-r'|"))
        {
            V4IElectronRepulsionDriver eri_drv;
            
            for (const auto& tint : integrals)
            {
                if ((tint[0] == 0) && (tint[2] > 0))
                {
                    const auto ctints = eri_drv.create_ket_hrr_recursion({tint, });
                    
                    tints.insert(ctints.cbegin(), ctints.cend());
                }
            }
        }
    }
    
    if (integral.prefixes_order() == std::vector<int>{2, 0, 0, 0})
    {
        // Electron repulsion integrals
        
        if (integral.integrand() == Operator("1/|r-r'|"))
        {
            V4IElectronRepulsionDriver eri_drv;
            
            for (const auto& tint : integrals)
            {
                if ((tint[0] == 0) && (tint[2] > 0))
                {
                    const auto ctints = eri_drv.create_ket_hrr_recursion({tint, });
                    
                    tints.insert(ctints.cbegin(), ctints.cend());
                }
            }
        }
    }
    
    return tints;
}

SI4CIntegrals
T4CGeomCPUGenerator::_generate_vrr_integral_group(const I4CIntegral& integral,
                                                  const SI4CIntegrals& bra_base_integrals,
                                                  const SI4CIntegrals& bra_rec_base_integrals,
                                                  const SI4CIntegrals& ket_base_integrals,
                                                  const SI4CIntegrals& ket_rec_base_integrals) const
{
    SI4CIntegrals tints;
    
    // Electron repulsion integrals
    
    if (integral.integrand() == Operator("1/|r-r'|"))
    {
        V4IElectronRepulsionDriver eri_drv;
        
        for (const auto& tint : bra_base_integrals)
        {
            if ((tint[0] == 0) && (tint[2] == 0))
            {
                const auto ctints = eri_drv.create_vrr_recursion({tint, });
                
                tints.insert(ctints.cbegin(), ctints.cend());
            }
        }
        
        for (const auto& tint : bra_rec_base_integrals)
        {
            if ((tint[0] == 0) && (tint[2] == 0))
            {
                const auto ctints = eri_drv.create_vrr_recursion({tint, });
                
                tints.insert(ctints.cbegin(), ctints.cend());
            }
        }
        
        for (const auto& tint : ket_base_integrals)
        {
            if ((tint[0] == 0) && (tint[2] == 0))
            {
                const auto ctints = eri_drv.create_vrr_recursion({tint, });
                
                tints.insert(ctints.cbegin(), ctints.cend());
            }
        }
        
        for (const auto& tint : ket_rec_base_integrals)
        {
            if ((tint[0] == 0) && (tint[2] == 0))
            {
                const auto ctints = eri_drv.create_vrr_recursion({tint, });
                
                tints.insert(ctints.cbegin(), ctints.cend());
            }
        }
    }
    
    return tints;
}

std::string
T4CGeomCPUGenerator::_file_name(const I4CIntegral& integral) const
{
    std::string label = "Rec" + integral.label();
    
    return t4c::integral_label(integral) + label;
}

void
T4CGeomCPUGenerator::_write_cpp_header(const SG4Terms&      cterms,
                                       const SG4Terms&      ckterms,
                                       const SG4Terms&      skterms,
                                       const SI4CIntegrals& vrr_integrals,
                                       const I4CIntegral& integral) const
{
    auto fname = _file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, true);
    
    _write_hpp_includes(fstream, ckterms, skterms, vrr_integrals, integral);
    
    _write_namespace(fstream, integral, true);
    
    T4CGeomDocuDriver docs_drv;
    
    T4CGeomDeclDriver decl_drv;
    
    T4CGeomFuncBodyDriver func_drv;

    docs_drv.write_doc_str(fstream, integral);
    
    decl_drv.write_func_decl(fstream, integral, false);
    
    func_drv.write_func_body(fstream, cterms, ckterms, skterms, vrr_integrals, integral);
    
    fstream << std::endl;

    _write_namespace(fstream, integral, false);
        
    _write_hpp_defines(fstream, integral, false);
    
    fstream.close();
}

void
T4CGeomCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                        const I4CIntegral&   integral,
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
T4CGeomCPUGenerator::_write_hpp_includes(      std::ofstream& fstream,
                                         const SG4Terms&      ckterms,
                                         const SG4Terms&      skterms,
                                         const SI4CIntegrals& vrr_integrals,
                                         const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "#include <array>"});
    
    lines.push_back({0, 0, 1, "#include <cstddef>"});
        
    lines.push_back({0, 0, 2, "#include <utility>"});
    
    std::set<std::string> labels;
        
    for (const auto& tint : vrr_integrals)
    {
        labels.insert(t4c::prim_file_name(tint));
    }
    
    for (const auto& term : ckterms)
    {
        const auto tint = term.second;
        
        if ((tint[0] == 0) && (tint[2] > 0) && (tint.prefixes().empty()))
        {
            labels.insert(t4c::ket_hrr_file_name(tint));
        }
    }
    
    for (const auto& term : skterms)
    {
        const auto tint = term.second;
        
        if ((tint[2] == integral[2]) && (tint[3] == integral[3]))
        {
            if (tint[0] > 0)
            {
                if (tint.prefixes().empty())
                {
                    labels.insert(t4c::bra_hrr_file_name(tint));
                }
                else
                {
                    labels.insert(t4c::bra_geom_hrr_file_name(tint));
                }
            }
            else
            {
                if (tint.prefixes_order() == std::vector<int>({2, 0, 0, 0}))
                {
                    labels.insert("ElectronRepulsionGeom2000ContrRecSXXX");
                }
                
                if (tint.prefixes_order() == std::vector<int>({0, 1, 0, 0}))
                {
                    labels.insert(t4c::bra_geom_hrr_file_name(tint));
                }
                
                if (tint.prefixes_order() == std::vector<int>({1, 1, 0, 0}))
                {
                    labels.insert(t4c::bra_geom_hrr_file_name(tint));
                }
                
                if (tint.prefixes_order() == std::vector<int>({0, 0, 1, 0}))
                {
                    labels.insert(t4c::ket_geom_hrr_file_name(tint));
                }
                
                if (tint.prefixes_order() == std::vector<int>({1, 0, 1, 0}))
                {
                    labels.insert(t4c::bra_geom_hrr_file_name(tint));
                }
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
    
    lines.push_back({0, 0, 1, "#include \"T2CUtils.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"BatchFunc.hpp\""});

    lines.push_back({0, 0, 2, "#include \"GtoPairBlock.hpp\""});
   
    ost::write_code_lines(fstream, lines);
}

void
T4CGeomCPUGenerator::_write_namespace(      std::ofstream& fstream,
                                      const I4CIntegral&   integral,
                                      const bool           start) const
{
    const auto label = t4c::namespace_label(integral);
    
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

void
T4CGeomCPUGenerator::_prune_terms_group(SG4Terms& terms) const
{
    SG4Terms new_terms;
    
    for (const auto& term : terms)
    {
        new_terms.insert(t4c::prune_term(term));
    }
    
    terms = new_terms;
}
