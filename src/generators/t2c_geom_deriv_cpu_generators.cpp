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

#include "t2c_geom_deriv_cpu_generators.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "file_stream.hpp"

#include "t2c_utils.hpp"
#include "t2c_geom_docs.hpp"

void
T2CGeomDerivCPUGenerator::generate(const int                 max_ang_mom,
                                   const std::array<int, 3>& geom_drvs) const
{
    if (geom_drvs[2] == 0)
    {
        for (int i = 0; i <= max_ang_mom; i++)
        {
            const auto integral = _get_integral({i, 0}, geom_drvs);
            
            const auto geom_integrals = t2c::get_geom_integrals(integral);
            
            _write_cpp_header(geom_integrals, integral, geom_drvs);
                                
            _write_cpp_file(geom_integrals, integral, geom_drvs);
            
            std::cout << " *** REFERENCE: " << integral.prefix_label() << " | " << integral.label() << " : " << geom_integrals.size() << std::endl;
           
            for (const auto& tint : geom_integrals)
            {
                std::cout << " <>" << tint.prefix_label() << " | " << tint.label() << std::endl;
            }
        }
    }
    else
    {
        for (int i = 0; i <= max_ang_mom; i++)
        {
            for (int j = 0; j <= max_ang_mom; j++)
            {
                const auto integral = _get_integral({i, j}, geom_drvs);
                
                const auto geom_integrals = t2c::get_geom_integrals(integral);
                
                _write_cpp_header(geom_integrals, integral, geom_drvs);
                                    
                _write_cpp_file(geom_integrals, integral, geom_drvs);
                
                std::cout << " *** REFERENCE: " << integral.prefix_label() << " | " << integral.label() << " : " << geom_integrals.size() << std::endl;
               
                for (const auto& tint : geom_integrals)
                {
                    std::cout << " <>" << tint.prefix_label() << " | " << tint.label() << std::endl;
                }
            }
        }
    }
}

I2CIntegral
T2CGeomDerivCPUGenerator::_get_integral(const std::array<int, 2>& ang_moms,
                                        const std::array<int, 3>& geom_drvs) const
{
    // bra and ket sides
    
    const auto bra = I1CPair("GA", ang_moms[0]);
    
    const auto ket = I1CPair("GB", ang_moms[1]);

    VOperators prefixes;

    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[0])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[2])));

    return I2CIntegral(bra, ket, Operator("R"), 0, prefixes);
}

void
T2CGeomDerivCPUGenerator::_write_cpp_header(const SI2CIntegrals&      geom_integrals,
                                            const I2CIntegral&        integral,
                                            const std::array<int, 3>& geom_drvs) const
{
    auto fname = t2c::geom_file_name(integral, geom_drvs) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, geom_drvs, true);

    _write_hpp_includes(fstream, integral, geom_drvs);

    _write_namespace(fstream, true);

    T2CGeomDocuDriver docs_drv;

    //T2CGeomDeclDriver decl_drv;

    docs_drv.write_doc_str(fstream, geom_integrals, integral, geom_drvs);

    //decl_drv.write_func_decl(fstream, geom_integrals, integral, geom_drvs, true);

    fstream << std::endl;

    _write_namespace(fstream, false);

    _write_hpp_defines(fstream, integral, geom_drvs, false);
    
    fstream.close();
}

void
T2CGeomDerivCPUGenerator::_write_hpp_defines(      std::ofstream&      fstream,
                                             const I2CIntegral&        integral,
                                             const std::array<int, 3>& geom_drvs,
                                             const bool                start) const
{
    auto fname = t2c::geom_file_name(integral, geom_drvs) + "_hpp";
    
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
T2CGeomDerivCPUGenerator::_write_hpp_includes(      std::ofstream&      fstream,
                                              const I2CIntegral&        integral,
                                              const std::array<int, 3>& geom_drvs) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"SimdArray.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CGeomDerivCPUGenerator::_write_namespace(      std::ofstream& fstream,
                                           const bool           start) const
{
    const auto label = t2c::geom_namespace_label();
    
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
T2CGeomDerivCPUGenerator::_write_cpp_file(const SI2CIntegrals&      geom_integrals,
                                          const I2CIntegral&        integral,
                                          const std::array<int, 3>& geom_drvs) const
{
    auto fname = t2c::geom_file_name(integral, geom_drvs) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
//    _write_cpp_includes(fstream, integral);
//
//    _write_namespace(fstream, true);
//
//    T4CGeomDeclDriver decl_drv;
//
//    decl_drv.write_func_decl(fstream, geom_integrals, integral, false);
//
//    T4CGeomFuncBodyDriver func_drv;
//
//    func_drv.write_func_body(fstream, geom_integrals, integral);
//
//    fstream << std::endl;
//
//    _write_namespace(fstream, false);
        
    fstream.close();
}
