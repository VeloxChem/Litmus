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

#include "t2c_cpu_generator.hpp"

#include <iostream>
#include <iterator>

#include "operator.hpp"
#include "string_formater.hpp"
#include "file_stream.hpp"
#include "spherical_momentum.hpp"

#include "t2c_docs.hpp"
#include "t2c_decl.hpp"
#include "t2c_body.hpp"
#include "t2c_prim_body.hpp"
#include "t2c_utils.hpp"

void
T2CCPUGenerator::generate(const std::string& label,
                          const int          angmom) const
{
    if (_is_available(label))
    {
        for (int i = 0; i <= angmom; i++)
        {
            for (int j = 0; j <= angmom; j++)
            {
                const auto integral = _get_integral(label, i, j);
                
                _write_cpp_header(integral);
                
                _write_cpp_file(integral);
            }
        }
    }
    else
    {
        std::cerr << "*** ERROR *** Unsupported type of two-center integral: ";
        
        std::cerr << label << " !!!" << std::endl;
        
        std::exit(EXIT_FAILURE);
    }
}

bool
T2CCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "overlap") return true;
    
    if (fstr::lowercase(label) == "kinetic_energy") return true;
    
    return false;
}

I2CIntegral
T2CCPUGenerator::_get_integral(const std::string& label,
                               const int          ang_a,
                               const int          ang_b) const
{
    const auto bra = I1CPair("GA", ang_a);
    
    const auto ket = I1CPair("GB", ang_b);
    
    // overlap integrals
    
    if (fstr::lowercase(label) == "overlap")
    {
        return I2CIntegral(bra, ket, Operator("1"));
    }
    
    // kinetic energy integrals
    
    if (fstr::lowercase(label) == "kinetic_energy")
    {
        return I2CIntegral(bra, ket, Operator("T"));
    }
    
    return I2CIntegral();
}

std::string
T2CCPUGenerator::_file_name(const I2CIntegral& integral) const
{
    return t2c::integral_label(integral) + "Rec" + integral.label();
}

void
T2CCPUGenerator::_write_cpp_header(const I2CIntegral& integral) const
{
    auto fname = _file_name(integral) + ".hpp";
    
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, true);
    
    _write_hpp_includes(fstream, integral);
    
    _write_namespace(fstream, integral, true);
    
    T2CDocuDriver docs_drv;
    
    T2CDeclDriver decl_drv;
    
    if (integral[0] == integral[1])
    {
        docs_drv.write_doc_str(fstream, integral, true);
        
        decl_drv.write_func_decl(fstream, integral, true, true);
    }
    
    docs_drv.write_doc_str(fstream, integral, true);
    
    decl_drv.write_func_decl(fstream, integral, false, true);
    
    _write_prim_funcs_to_cpp_header(fstream, integral); 

    _write_namespace(fstream, integral, false);
    
    _write_hpp_defines(fstream, integral, false);

    fstream.close();
}

void
T2CCPUGenerator::_write_cpp_file(const I2CIntegral& integral) const
{
    auto fname = _file_name(integral) + ".cpp";
    
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_cpp_includes(fstream, integral);
    
    _write_namespace(fstream, integral, true);
    
    T2CDeclDriver decl_drv;
    
    T2CFuncBodyDriver func_drv;
    
    if (integral[0] == integral[1])
    {
        decl_drv.write_func_decl(fstream, integral, true, false);
        
        func_drv.write_func_body(fstream, integral, true);
    }
    
    decl_drv.write_func_decl(fstream, integral, false, false);
    
    func_drv.write_func_body(fstream, integral, true);
    
    _write_prim_funcs_to_cpp_file(fstream, integral);

    _write_namespace(fstream, integral, false);

    fstream.close();
}

void
T2CCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                    const I2CIntegral&   integral,
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

void
T2CCPUGenerator::_write_hpp_includes(      std::ofstream& fstream,
                                     const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include <cstdint>"});
    
    lines.push_back({0, 0, 1, "#include \"GtoBlock.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"SubMatrix.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"SimdTypes.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"MatrixType.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_cpp_includes(      std::ofstream& fstream,
                                     const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + _file_name(integral) +  ".hpp\""});
    
    lines.push_back({0, 0, 2, "#include <cmath>"});
    
    lines.push_back({0, 0, 1, "#include \"BatchFunc.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"MathConst.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"T2CDistributor.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_namespace(      std::ofstream& fstream,
                                  const I2CIntegral&   integral,
                                  const bool           start) const
{
    const auto label = t2c::namespace_label(integral);
    
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
T2CCPUGenerator::_write_prim_funcs_to_cpp_header(      std::ofstream& fstream,
                                                 const I2CIntegral&   integral) const
{
    T2CDocuDriver docs_drv;
    
    T2CDeclDriver decl_drv;
    
    if (integral.is_simple_integrand() && integral.is_simple())
    {
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            docs_drv.write_prim_doc_str(fstream, integral);
            
            decl_drv.write_prim_func_decl(fstream, integral, true);
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                for (const auto& bcomp: Tensor(integral[0]).components())
                {
                    docs_drv.write_prim_doc_str(fstream, bcomp, integral, true);
                    
                    decl_drv.write_prim_func_decl(fstream, bcomp, integral, true, true);
                }
            }
            else
            {
                for (const auto& kcomp: Tensor(integral[1]).components())
                {
                    docs_drv.write_prim_doc_str(fstream, kcomp, integral, false);
                    
                    decl_drv.write_prim_func_decl(fstream, kcomp, integral, false, true);
                }
            }
        }
    }
    else
    {
        for (const auto& bcomp: Tensor(integral[0]).components())
        {
            for (const auto& kcomp: Tensor(integral[1]).components())
            {
                docs_drv.write_prim_doc_str(fstream, bcomp, kcomp, integral);
                
                decl_drv.write_prim_func_decl(fstream, bcomp, kcomp, integral, true);
            }
        }
    }
}

void
T2CCPUGenerator::_write_prim_funcs_to_cpp_file(      std::ofstream& fstream,
                                               const I2CIntegral&   integral) const
{
    T2CDeclDriver decl_drv;
    
    T2CPrimFuncBodyDriver func_drv;

    if (integral.is_simple_integrand() && integral.is_simple())
    {
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            decl_drv.write_prim_func_decl(fstream, integral, false);
            
            func_drv.write_prim_func_body(fstream, integral);
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                for (const auto& bcomp: Tensor(integral[0]).components())
                {
                    decl_drv.write_prim_func_decl(fstream, bcomp, integral, true, false);
                    
                    func_drv.write_prim_func_body(fstream, bcomp, integral, true);
                }
            }
            else
            {
                for (const auto& kcomp: Tensor(integral[1]).components())
                {
                    decl_drv.write_prim_func_decl(fstream, kcomp, integral, false, false);
                    
                    func_drv.write_prim_func_body(fstream, kcomp, integral, false);
                }
            }
        }
    }
    else
    {
        for (const auto& bcomp: Tensor(integral[0]).components())
        {
            for (const auto& kcomp: Tensor(integral[1]).components())
            {
                decl_drv.write_prim_func_decl(fstream, bcomp, kcomp, integral, false);
                
                func_drv.write_prim_func_body(fstream, bcomp, kcomp, integral);
            }
        }
    }
}
