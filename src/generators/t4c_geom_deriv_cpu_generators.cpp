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

#include "t4c_geom_deriv_cpu_generators.hpp"

#include <iostream>

#include "file_stream.hpp"

#include "t4c_utils.hpp"
#include "t4c_geom_docs.hpp"
#include "t4c_geom_decl.hpp"
#include "t4c_geom_body.hpp"
#include "v4i_center_driver.hpp"

void
T4CGeomDerivCPUGenerator::generate(const int                 max_ang_mom,
                                   const std::array<int, 4>& geom_drvs) const
{
    
    for (int i = 0; i <= max_ang_mom; i++)
    {
        for (int j = i; j <= max_ang_mom; j++)
        {
            for (int k = 0; k <= max_ang_mom; k++)
            {
                for (int l = k; l <= max_ang_mom; l++)
                {
                    const auto integral = _get_integral({i, j, k, l}, geom_drvs);
                                        
                    const auto geom_integrals = t4c::get_geom_integrals(integral);
                                        
                    _write_cpp_header(geom_integrals, integral);
                                        
                    _write_cpp_file(geom_integrals, integral);
                }
            }
        }
    }
}

I4CIntegral
T4CGeomDerivCPUGenerator::_get_integral(const std::array<int, 4>& ang_moms,
                                        const std::array<int, 4>& geom_drvs) const
{
    // bra and ket sides

    const auto bpair = I2CPair("GA", ang_moms[0], "GB", ang_moms[1]);

    const auto kpair = I2CPair("GC", ang_moms[2], "GD", ang_moms[3]);

    VOperators prefixes;

    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[0])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[1])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[2])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[3])));

    return I4CIntegral(bpair, kpair, Operator("1"), 0, prefixes);
}

SI4CIntegrals
T4CGeomDerivCPUGenerator::_generate_geom_integral_group(const I4CIntegral& integral) const
{
    V4ICenterDriver geom_drv;
    
    SI4CIntegrals ref_tints;
    
    for (const auto& tint : geom_drv.apply_bra_ket_vrr(integral))
    {
        ref_tints.insert(tint.base());
    }
    
    return ref_tints;
}

void
T4CGeomDerivCPUGenerator::_write_cpp_header(const SI4CIntegrals& geom_integrals,
                                            const I4CIntegral&   integral) const
{
    auto fname = t4c::geom_file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, true);
    
    _write_hpp_includes(fstream,  integral);
    
    _write_namespace(fstream, true);
    
    T4CGeomDocuDriver docs_drv;

    T4CGeomDeclDriver decl_drv;

    docs_drv.write_doc_str(fstream, geom_integrals, integral);

    decl_drv.write_func_decl(fstream, geom_integrals, integral, true);

    fstream << std::endl;
    
    _write_namespace(fstream, false);

    _write_hpp_defines(fstream, integral, false);
    
    fstream.close();
}

void
T4CGeomDerivCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                             const I4CIntegral&   integral,
                                             const bool           start) const
{
    auto fname = t4c::geom_file_name(integral) + "_hpp";
    
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
T4CGeomDerivCPUGenerator::_write_hpp_includes(      std::ofstream& fstream,
                                              const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"SimdArray.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T4CGeomDerivCPUGenerator::_write_namespace(      std::ofstream& fstream,
                                           const bool           start) const
{
    const auto label = t4c::geom_namespace_label();
    
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
T4CGeomDerivCPUGenerator::_write_cpp_file(const SI4CIntegrals& geom_integrals,
                                          const I4CIntegral&   integral) const
{
    auto fname = t4c::geom_file_name(integral) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_cpp_includes(fstream, integral);

    _write_namespace(fstream, true);

    T4CGeomDeclDriver decl_drv;
    
    decl_drv.write_func_decl(fstream, geom_integrals, integral, false);

    T4CGeomFuncBodyDriver func_drv;

    func_drv.write_func_body(fstream, geom_integrals, integral);
    
    fstream << std::endl;
    
    _write_namespace(fstream, false);
        
    fstream.close();
}

void
T4CGeomDerivCPUGenerator::_write_cpp_includes(      std::ofstream& fstream,
                                              const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + t4c::geom_file_name(integral) +  ".hpp\""});
    
    ost::write_code_lines(fstream, lines);
}
