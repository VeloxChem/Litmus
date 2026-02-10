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

#include "t2c_hrr_cpu_generators.hpp"

#include "t2c_utils.hpp"
#include "string_formater.hpp"
#include "file_stream.hpp"
#include "t2c_hrr_docs.hpp"
#include "t2c_hrr_decl.hpp"
#include "t2c_hrr_body.hpp"

void
T2CHRRCPUGenerator::generate(const int max_ang_mom) const
{
    #pragma omp parallel
    {
        #pragma omp single nowait
        {
            for (int i = 1; i <= 2 * max_ang_mom; i++)
            {
                for (int j = 1; j <= 2 * max_ang_mom; j++)
                {
                    if ((i + j) > 2 * max_ang_mom) continue;
                 
                    #pragma omp task firstprivate(i,j)
                    {
                        const auto integral = _get_integral({i, j});

                        _write_hrr_cpp_header(integral);
                                
                        _write_hrr_cpp_file(integral);
                    }
                }
            }
        }
    }
}

I2CIntegral
T2CHRRCPUGenerator::_get_integral(const std::array<int, 2>& ang_moms) const
{
    // bra and ket sides
    
    const auto bra = I1CPair("GA", ang_moms[0]);
    
    const auto ket = I1CPair("GB", ang_moms[1]);

    return I2CIntegral(bra, ket, Operator(), 0, {});
}

void
T2CHRRCPUGenerator::_write_hrr_cpp_header(const I2CIntegral& integral) const
{
    auto fname = t2c::hrr_file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, true);
    
    _write_hrr_hpp_includes(fstream, integral);
    
    _write_namespace(fstream, integral, true);
    
    T2CHRRDocuDriver docs_drv;
    
    docs_drv.write_doc_str(fstream, integral);
    
    T2CHRRDeclDriver decl_drv;
    
    decl_drv.write_func_decl(fstream, integral, true);
    
    _write_namespace(fstream, integral, false);
    
    _write_hpp_defines(fstream, integral, false);
    
    fstream << std::endl;
    
    fstream.close();
}

void
T2CHRRCPUGenerator::_write_hrr_cpp_file(const I2CIntegral& integral) const
{
    auto fname = t2c::hrr_file_name(integral) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_hrr_cpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T2CHRRDeclDriver decl_drv;
    
    decl_drv.write_func_decl(fstream, integral, false);

    T2CHRRFuncBodyDriver func_drv;

    func_drv.write_func_body(fstream, integral);
    
    fstream << std::endl;
    
    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
T2CHRRCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                       const I2CIntegral&   integral,
                                       const bool           start) const
{
    auto fname = t2c::hrr_file_name(integral) + "_hpp";
    
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
T2CHRRCPUGenerator::_write_hrr_hpp_includes(      std::ofstream& fstream,
                                            const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"SimdArray.hpp\""});
        
    ost::write_code_lines(fstream, lines);
}

void
T2CHRRCPUGenerator::_write_hrr_cpp_includes(      std::ofstream& fstream,
                                            const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + t2c::hrr_file_name(integral) +  ".hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CHRRCPUGenerator::_write_namespace(      std::ofstream& fstream,
                                     const I2CIntegral&   integral,
                                     const bool           start) const
{
    const std::string label = "t2chrr";
    
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
