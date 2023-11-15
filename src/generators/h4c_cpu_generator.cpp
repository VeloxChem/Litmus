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

#include "h4c_cpu_generator.hpp"

#include <iostream>

#include "t4c_full_docs.hpp"
#include "t4c_full_decl.hpp"
#include "t4c_full_body.hpp"
#include "t4c_full_prim_body.hpp"
#include "t4c_utils.hpp"
#include "string_formater.hpp"

void
H4CCPUGenerator::generate(const std::string& label,
                          const int          angmom) const
{
    if (_is_available(label))
    {
        for (int i = 1; i <= angmom; i++)
        {
            for (int j = i; j <= angmom; j++)
            {
                for (int k = 1; k <= angmom; k++)
                {
                    for (int l = k; l <= angmom; l++)
                    {
                        #pragma omp parallel
                        {
                            #pragma omp single nowait
                            {
                                if ((i + j + k + l) > 0)
                                {
                                    const auto integral = _get_integral(label, i, j, k, l);
                                    
                                    _write_cpp_prim_headers(integral);
                                    
                                    _write_cpp_prim_files(integral);
                                }
                            }
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
H4CCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "electron repulsion") return true;
    
    return false;
}

I4CIntegral
H4CCPUGenerator::_get_integral(const std::string& label,
                               const int          ang_a,
                               const int          ang_b,
                               const int          ang_c,
                               const int          ang_d) const
{
    // bra and ket sides
    
    const auto bpair = I2CPair("GA", ang_a, "GB", ang_b);
    
    const auto kpair = I2CPair("GC", ang_c, "GD", ang_d);
    
    // electron repulsion integrals
    
    if (fstr::lowercase(label) == "electron repulsion")
    {
        return I4CIntegral(bpair, kpair, Operator("1/|r-r'|"));
    }
    
    return I4CIntegral();
}

std::string
H4CCPUGenerator::_file_name(const I4CIntegral& integral) const
{
    return t4c::integral_label(integral) + "HRRRec" + integral.label();
}

void
H4CCPUGenerator::_write_cpp_prim_headers(const I4CIntegral& integral) const
{
    if ((integral[0] + integral[1] + integral[2] + integral[3]) == 0) return;
    
    for (const auto& tcomp : integral.components<T2CPair, T2CPair>())
    {
        #pragma omp task firstprivate(tcomp, integral)
        {
            auto fname = t4c::full_hrr_file_name(tcomp, integral) + ".hpp";
            
            std::ofstream fstream;
                   
            fstream.open(fname.c_str(), std::ios_base::trunc);
            
            _write_hpp_prim_defines(fstream, t4c::full_hrr_file_name(tcomp, integral), true);
            
            _write_hpp_prim_includes(fstream, integral);
            
            _write_namespace(fstream, integral, true);
            
            T4CFullDocuDriver docs_drv;
            
            docs_drv.write_hrr_doc_str(fstream, tcomp, integral);
            
            T4CFullDeclDriver decl_drv;
            
            decl_drv.write_hrr_func_decl(fstream, tcomp, integral, true);
            
            _write_namespace(fstream, integral, false);
            
            _write_hpp_prim_defines(fstream, t4c::full_hrr_file_name(tcomp, integral), false);

            fstream.close();
        }
    }
}

void
H4CCPUGenerator::_write_cpp_prim_files(const I4CIntegral& integral) const
{
    if ((integral[0] + integral[1] + integral[2] + integral[3]) == 0) return;
    
    for (const auto& tcomp : integral.components<T2CPair, T2CPair>())
    {
        #pragma omp task firstprivate(tcomp, integral)
        {
            auto fname = t4c::full_hrr_file_name(tcomp, integral) + ".cpp";
            
            std::ofstream fstream;
                   
            fstream.open(fname.c_str(), std::ios_base::trunc);
            
            _write_cpp_prim_includes(fstream, tcomp, integral);
            
            _write_namespace(fstream, integral, true);
            
            T4CFullDeclDriver decl_drv;
            
            T4CFullPrimFuncBodyDriver func_drv;
            
            decl_drv.write_hrr_func_decl(fstream, tcomp, integral, false);
            
            func_drv.write_hrr_func_body(fstream, tcomp, integral);
            
            _write_namespace(fstream, integral, false);
            
            fstream.close();
        }
    }
}

void
H4CCPUGenerator::_write_namespace(      std::ofstream& fstream,
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
H4CCPUGenerator::_write_hpp_prim_defines(      std::ofstream& fstream,
                                         const std::string&   fname,
                                         const bool           start) const
{
    const auto flabel = fname + "_hpp";
    
    auto lines = VCodeLines();
 
    if (start)
    {
        lines.push_back({0, 0, 1, "#ifndef " + flabel});
        
        lines.push_back({0, 0, 2, "#define " + flabel});
    }
    else
    {
        lines.push_back({0, 0, 1, "#endif /* " + flabel + " */"});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
H4CCPUGenerator::_write_hpp_prim_includes(      std::ofstream& fstream,
                                          const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include <cstdint>"});
        
    lines.push_back({0, 0, 1, "#include \"Point.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"SimdTypes.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
H4CCPUGenerator::_write_cpp_prim_includes(      std::ofstream& fstream,
                                          const T4CIntegral&   component,
                                          const I4CIntegral&   integral) const
{
    auto fname = t4c::full_hrr_file_name(component, integral) + ".hpp";
    
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + fname + "\""});
    
    lines.push_back({0, 0, 2, "#include <cmath>"});
  
    lines.push_back({0, 0, 1, "#include \"BoysFunc.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"MathConst.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
H4CCPUGenerator::_add_prim_call_includes(      VCodeLines&  lines,
                                         const I4CIntegral& integral) const
{
    for (const auto& tcomp : integral.components<T2CPair, T2CPair>())
    {
        lines.push_back({0, 0, 1, "#include \"" + t4c::full_hrr_file_name(tcomp, integral) + ".hpp\""});
    }
    
    lines.push_back({0, 0, 1, ""});
}