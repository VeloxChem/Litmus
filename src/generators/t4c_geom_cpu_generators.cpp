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
#include "t4c_docs.hpp"
#include "t4c_decl.hpp"
#include "t4c_body.hpp"
#include "v4i_center_driver.hpp"
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
            for (int j = i; j <= max_ang_mom; j++)
            {
                for (int k = 0; k <= max_ang_mom; k++)
                {
                    for (int l = k; l <= max_ang_mom; l++)
                    {
                        const auto integral = _get_integral(label, {i, j, k, l}, geom_drvs);
                        
                        const auto geom_integrals = _generate_geom_integral_group(integral);
                        
                        const auto vrr_integrals = _generate_vrr_integral_group(integral, geom_integrals);
                        
                        _write_cpp_header(geom_integrals, vrr_integrals, integral);
                        
                        std::cout << " *** REFERENCE: " << integral.prefix_label() << " | " << integral.label() << std::endl;
                        
                        for (const auto& tint : geom_integrals)
                        {
                            std::cout << " <>" << tint.prefix_label() << " | " << tint.label() << std::endl;
                        }
                        
                        std::cout << " --- VRR --- " << std::endl;
                        
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
    V4ICenterDriver geom_drv;

    return geom_drv.apply_bra_ket_vrr(integral);
}

SI4CIntegrals
T4CGeomCPUGenerator::_generate_vrr_integral_group(const I4CIntegral& integral,
                                                  const SI4CIntegrals& integrals) const
{
    SI4CIntegrals ref_tints;
    
    for (const auto& tint : integrals)
    {
        ref_tints.insert(tint.base());
    }
    
    SI4CIntegrals tints;
    
    // Electron repulsion integrals
    
    if (integral.integrand() == Operator("1/|r-r'|"))
    {
        V4IElectronRepulsionDriver eri_drv;
        
        tints = eri_drv.create_full_vrr_recursion(ref_tints);
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
T4CGeomCPUGenerator::_write_cpp_header(const SI4CIntegrals& geom_integrals,
                                       const SI4CIntegrals& vrr_integrals,
                                       const I4CIntegral&   integral) const
{
    auto fname = _file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, true);
    
    _write_hpp_includes(fstream, geom_integrals, vrr_integrals, integral);
    
    _write_namespace(fstream, integral, true);
    
    T4CDocuDriver docs_drv;
    
    T4CDeclDriver decl_drv;
    
    T4CFuncBodyDriver func_drv;

    docs_drv.write_doc_str(fstream, integral, false);
    
    decl_drv.write_func_decl(fstream, integral, false, false);
    
    func_drv.write_geom_func_body(fstream, geom_integrals, vrr_integrals, integral);
    
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
                                         const SI4CIntegrals& geom_integrals,
                                         const SI4CIntegrals& vrr_integrals,
                                         const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include <array>"});
    
    std::set<std::string> labels;
    
    for (const auto& tint : vrr_integrals)
    {
        if ((tint[0] + tint[2]) == 0)
        {
            labels.insert(t4c::prim_file_name(tint));
        }
    }
    
//    for (const auto& tint : ket_integrals)
//    {
//        if ((tint[0] == 0) && (tint[2] > 0))
//        {
//            labels.insert(t4c::ket_hrr_file_name(tint));
//        }
//    }
//
//    for (const auto& tint : bra_integrals)
//    {
//        if ((tint[0] > 0) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
//        {
//            labels.insert(t4c::bra_hrr_file_name(tint));
//        }
//    }
    
    for (const auto& label : labels)
    {
        lines.push_back({0, 0, 1, "#include \"" + label + ".hpp\""});
    }
    
    lines.push_back({0, 0, 1, "#include \"SimdArray.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"BoysFunc.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T4CUtils.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T2CUtils.hpp\""});
    
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
