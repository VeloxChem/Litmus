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

#include "t4c_diag_cpu_generator.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "t4c_utils.hpp"

void
T4CDiagCPUGenerator::generate(const std::string& label,
                              const int          angmom) const
{
    if (_is_available(label))
    {
        for (int i = 0; i <= angmom; i++)
        {
            for (int j = i; j <= angmom; j++)
            {
                #pragma omp parallel
                {
                    #pragma omp single nowait
                    {
                        const auto integral = _get_integral(label, i, j);
                        
                        #pragma omp task firstprivate(integral)
                        _write_cpp_header(integral);
                        
                        //#pragma omp task firstprivate(integral)
                        //_write_cpp_file(integral);
                        
                        //_write_cpp_prim_headers(integral);
                        
                        //_write_cpp_prim_files(integral);
                    }
                }
            }
        }
    }
    else
    {
        std::cerr << "*** ERROR *** Unsupported type of diagonal four-center integral: ";
        
        std::cerr << label << " !!!" << std::endl;
        
        std::exit(EXIT_FAILURE);
    }
    
}

bool
T4CDiagCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "electron repulsion") return true;
    
    return false;
}

I4CIntegral
T4CDiagCPUGenerator::_get_integral(const std::string& label,
                                   const int          ang_a,
                                   const int          ang_b) const
{
    // bra and ket sides
    
    const auto bpair = I2CPair("GA", ang_a, "GB", ang_b);
    
    const auto kpair = I2CPair("GC", ang_a, "GD", ang_b);
    
    // electron repulsion integrals
    
    if (fstr::lowercase(label) == "electron repulsion")
    {
        return I4CIntegral(bpair, kpair, Operator("1/|r-r'|"));
    }
    
    return I4CIntegral();
}

std::string
T4CDiagCPUGenerator::_file_name(const I4CIntegral& integral) const
{
    return t4c::integral_label(integral) + "DiagRec" + integral.label();
}

void
T4CDiagCPUGenerator::_write_cpp_header(const I4CIntegral& integral) const
{
    auto fname = _file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_hpp_defines(fstream, integral, true);
        
    //_write_hpp_includes(fstream, integral);
        
    //_write_namespace(fstream, integral, true);
        
    //T2CDocuDriver docs_drv;
        
    //T2CDeclDriver decl_drv;
        
    //if ((integral[0] == integral[1]) && integral.is_simple())
    //{
    //    docs_drv.write_doc_str(fstream, integral, true);
    //
    //    decl_drv.write_func_decl(fstream, integral, true, true);
    //}
        
    //docs_drv.write_doc_str(fstream, integral, false);
        
    //decl_drv.write_func_decl(fstream, integral, false, true);

    //_write_namespace(fstream, integral, false);
        
    _write_hpp_defines(fstream, integral, false);

    fstream.close();
}

void
T4CDiagCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                        const I4CIntegral&   integral,
                                        const bool           start) const
{
    const auto fname = _file_name(integral) + "_hpp";
    
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
