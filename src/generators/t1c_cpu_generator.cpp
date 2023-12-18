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

#include "t1c_cpu_generator.hpp"

#include "string_formater.hpp"
#include "t1c_docs.hpp"
#include "t1c_decl.hpp"
#include "t1c_body.hpp"

void
T1CCPUGenerator::generate(const std::string& label,
                          const int          angmom,
                          const int          gdrv) const
{
    for (int i = 0; i <= angmom; i++)
    {
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
               #pragma omp task firstprivate(i)
               _write_cpp_header(i, gdrv);

               #pragma omp task firstprivate(i)
               _write_cpp_file(i, gdrv);
            }
        }
    }
}

std::string
T1CCPUGenerator::_file_name(const int angmom) const
{
    return "GtoValuesGeomRec" + Tensor(angmom).label();
}

void
T1CCPUGenerator::_write_cpp_header(const int angmom,
                                   const int gdrv) const
{
    auto fname = _file_name(angmom) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, angmom, true);
        
    _write_hpp_includes(fstream);
        
    _write_namespace(fstream, true);
    
    T1CDocuDriver docs_drv;
    
    T1CDeclDriver decl_drv;
    
    for (int i = 1; i <= gdrv; i++)
    {
        docs_drv.write_doc_str(fstream, angmom, i);
       
        decl_drv.write_func_decl(fstream, angmom, i, true);
    }

    _write_namespace(fstream, false);
        
    _write_hpp_defines(fstream, angmom, false);

    fstream.close();
}

void
T1CCPUGenerator::_write_cpp_file(const int angmom,
                                 const int gdrv) const
{
    auto fname = _file_name(angmom) + ".cpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_cpp_includes(fstream, angmom);
        
    _write_namespace(fstream, true);
        
    T1CDeclDriver decl_drv;
        
    T1CFuncBodyDriver func_drv;
    
    for (int i = 1; i <= gdrv; i++)
    {
        decl_drv.write_func_decl(fstream, angmom, i, false);
    
        func_drv.write_func_body(fstream, angmom, i);
    }
    
    _write_namespace(fstream, false);

    fstream.close();
}

void
T1CCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                    const int            angmom,
                                    const bool           start) const
{
    const auto fname = _file_name(angmom) + "_hpp";
    
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
T1CCPUGenerator::_write_hpp_includes(std::ofstream& fstream) const
{
    auto lines = VCodeLines();
  
    lines.push_back({0, 0, 1, "#include <cstdint>"});
        
    lines.push_back({0, 0, 2, "#include <vector>"});
    
    lines.push_back({0, 0, 1, "#include \"GtoBlock.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"Matrix.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T1CCPUGenerator::_write_namespace(      std::ofstream& fstream,
                                  const bool           start) const
{
    const auto label = std::string("gtoval");
    
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
T1CCPUGenerator::_write_cpp_includes(      std::ofstream& fstream,
                                     const int            angmom) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + _file_name(angmom) +  ".hpp\""});
    
    lines.push_back({0, 0, 2, "#include <cmath>"});
    
    lines.push_back({0, 0, 1, "#include \"DftFunc.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"MathFunc.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"MatrixFunc.hpp\""});

    ost::write_code_lines(fstream, lines);
}
