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
#include "spherical_momentum.hpp"

#include "t2c_docs.hpp"
#include "t2c_decl.hpp"
#include "t2c_body.hpp"
#include "t2c_prim_body.hpp"
#include "t2c_utils.hpp"

void
T2CCPUGenerator::generate(const std::string& label,
                          const int          angmom,
                          const int          op_gdrv) const
{
    if (_is_available(label))
    {
        for (int i = 0; i <= angmom; i++)
        {
            for (int j = 0; j <= angmom; j++)
            {
                #pragma omp parallel
                {
                    #pragma omp single nowait
                    {
                        const auto integral = _get_integral(label, i, j, op_gdrv);
                        
                        #pragma omp task firstprivate(integral)
                        _write_cpp_header(integral);
                        
                        #pragma omp task firstprivate(integral)
                        _write_cpp_file(integral);
                        
                        _write_cpp_prim_headers(integral);
                        
                        _write_cpp_prim_files(integral);
                    }
                }
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
    
    if (fstr::lowercase(label) == "kinetic energy") return true;
    
    if (fstr::lowercase(label) == "nuclear potential") return true;
    
    if (fstr::lowercase(label) == "nuclear potential geom") return true;
    
    return false;
}

I2CIntegral
T2CCPUGenerator::_get_integral(const std::string& label,
                               const int          ang_a,
                               const int          ang_b,
                               const int          op_gdrv) const
{
    const auto bra = I1CPair("GA", ang_a);
    
    const auto ket = I1CPair("GB", ang_b);
    
    // overlap integrals
    
    if (fstr::lowercase(label) == "overlap")
    {
        return I2CIntegral(bra, ket, Operator("1"));
    }
    
    // kinetic energy integrals
    
    if (fstr::lowercase(label) == "kinetic energy")
    {
        return I2CIntegral(bra, ket, Operator("T"));
    }
    
    if (op_gdrv > 0)
    {
        // nuclear potential geometrical derivative integrals
        
        if (fstr::lowercase(label) == "nuclear potential")
        {
            return I2CIntegral(bra, ket, Operator("AG", Tensor(op_gdrv)));
        }
    }
    else
    {
        // nuclear potential integrals
        
        if (fstr::lowercase(label) == "nuclear potential")
        {
            return I2CIntegral(bra, ket, Operator("A"));
        }
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
    
    docs_drv.write_doc_str(fstream, integral, false);
    
    decl_drv.write_func_decl(fstream, integral, false, true);

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
    
    func_drv.write_func_body(fstream, integral, false);

    _write_namespace(fstream, integral, false);

    fstream.close();
}

void
T2CCPUGenerator::_write_cpp_prim_headers(const I2CIntegral& integral) const
{
    if (integral.is_simple_integrand() && integral.is_simple())
    {
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            #pragma omp task firstprivate(integral)
            _write_cpp_prim_header(integral);
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                for (const auto& bcomp: Tensor(integral[0]).components())
                {
                    #pragma omp task firstprivate(bcomp, integral)
                    _write_cpp_prim_header(bcomp, integral, true);
                }
            }
            else
            {
                for (const auto& kcomp: Tensor(integral[1]).components())
                {
                    #pragma omp task firstprivate(kcomp, integral)
                    _write_cpp_prim_header(kcomp, integral, false);
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
                #pragma omp task firstprivate(bcomp, kcomp, integral)
                _write_cpp_prim_header(bcomp, kcomp, integral);
            }
        }
    }
}

void
T2CCPUGenerator::_write_cpp_prim_files(const I2CIntegral& integral) const
{
    if (integral.is_simple_integrand() && integral.is_simple())
    {
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            #pragma omp task firstprivate(integral)
            _write_cpp_prim_file(integral);
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                for (const auto& bcomp: Tensor(integral[0]).components())
                {
                    #pragma omp task firstprivate(bcomp, integral)
                    _write_cpp_prim_file(bcomp, integral, true);
                }
            }
            else
            {
                for (const auto& kcomp: Tensor(integral[1]).components())
                {
                    #pragma omp task firstprivate(kcomp, integral)
                    _write_cpp_prim_file(kcomp, integral, false);
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
                #pragma omp task firstprivate(bcomp, kcomp, integral)
                _write_cpp_prim_file(bcomp, kcomp, integral);
            }
        }
    }
}

void
T2CCPUGenerator::_write_cpp_prim_header(const I2CIntegral& integral) const
{
    auto fname = t2c::prim_file_name(integral) + ".hpp";
    
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_prim_defines(fstream, t2c::prim_file_name(integral), true);
    
    _write_hpp_prim_includes(fstream, integral);
    
    _write_namespace(fstream, integral, true);
    
    T2CDocuDriver docs_drv;
    
    docs_drv.write_prim_doc_str(fstream, integral);
    
    T2CDeclDriver decl_drv;
    
    decl_drv.write_prim_func_decl(fstream, integral, true);
   
    _write_namespace(fstream, integral, false);
    
    _write_hpp_prim_defines(fstream, t2c::prim_file_name(integral), false);
    
    fstream.close();
}

void
T2CCPUGenerator::_write_cpp_prim_header(const TensorComponent& component,
                                        const I2CIntegral&     integral,
                                        const bool             bra_first) const
{
    auto fname = t2c::prim_file_name(component, integral, bra_first) + ".hpp";
    
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_prim_defines(fstream, t2c::prim_file_name(component, integral, bra_first), true);
    
    _write_hpp_prim_includes(fstream, integral);
    
    _write_namespace(fstream, integral, true);
    
    T2CDocuDriver docs_drv;
    
    docs_drv.write_prim_doc_str(fstream, component, integral, bra_first);
    
    T2CDeclDriver decl_drv;
    
    decl_drv.write_prim_func_decl(fstream, component, integral, bra_first, true);
    
    _write_namespace(fstream, integral, false);
    
    _write_hpp_prim_defines(fstream, t2c::prim_file_name(component, integral, bra_first), false);

    fstream.close();
}

void
T2CCPUGenerator::_write_cpp_prim_header(const TensorComponent& bra_component,
                                        const TensorComponent& ket_component,
                                        const I2CIntegral&     integral) const
{
    auto fname = t2c::prim_file_name(bra_component, ket_component, integral) + ".hpp";
    
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_prim_defines(fstream, t2c::prim_file_name(bra_component, ket_component, integral), true);
    
    _write_hpp_prim_includes(fstream, integral);
    
    _write_namespace(fstream, integral, true);
    
    T2CDocuDriver docs_drv;
    
    docs_drv.write_prim_doc_str(fstream, bra_component, ket_component, integral);
    
    T2CDeclDriver decl_drv;

    decl_drv.write_prim_func_decl(fstream, bra_component, ket_component, integral, true);
    
    _write_namespace(fstream, integral, false);
    
    _write_hpp_prim_defines(fstream, t2c::prim_file_name(bra_component, ket_component, integral), false);
    
    fstream.close();
}

void
T2CCPUGenerator::_write_cpp_prim_file(const I2CIntegral& integral) const
{
    auto fname = t2c::prim_file_name(integral) + ".cpp";
    
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_cpp_prim_includes(fstream, integral);
    
    _write_namespace(fstream, integral, true);
    
    T2CDeclDriver decl_drv;
    
    decl_drv.write_prim_func_decl(fstream, integral, false);
    
    T2CPrimFuncBodyDriver func_drv;
    
    func_drv.write_prim_func_body(fstream, integral);
   
    _write_namespace(fstream, integral, false);
    
    fstream.close();
}

void
T2CCPUGenerator::_write_cpp_prim_file(const TensorComponent& component,
                                      const I2CIntegral&     integral,
                                      const bool             bra_first) const
{
    auto fname = t2c::prim_file_name(component, integral, bra_first) + ".cpp";
    
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_cpp_prim_includes(fstream, component, integral, bra_first);
    
    _write_namespace(fstream, integral, true);
    
    T2CDeclDriver decl_drv;
    
    decl_drv.write_prim_func_decl(fstream, component, integral, bra_first, false);
    
    T2CPrimFuncBodyDriver func_drv;
    
    func_drv.write_prim_func_body(fstream, component, integral, bra_first);
    
    _write_namespace(fstream, integral, false);

    fstream.close();
}

void
T2CCPUGenerator::_write_cpp_prim_file(const TensorComponent& bra_component,
                                      const TensorComponent& ket_component,
                                      const I2CIntegral&     integral) const
{
    auto fname = t2c::prim_file_name(bra_component, ket_component, integral) + ".cpp";
    
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_cpp_prim_includes(fstream, bra_component, ket_component, integral);
    
    _write_namespace(fstream, integral, true);
    
    T2CDeclDriver decl_drv;

    decl_drv.write_prim_func_decl(fstream, bra_component, ket_component, integral, false);
    
    T2CPrimFuncBodyDriver func_drv;
    
    func_drv.write_prim_func_body(fstream, bra_component, ket_component, integral);
    
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
T2CCPUGenerator::_write_hpp_prim_defines(      std::ofstream& fstream,
                                         const std::string&   fname,
                                         const bool           start) const
{
    const auto flabel = fname + "_hpp";
    
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
    
    if (integral[0] == integral[1])
    {
        lines.push_back({0, 0, 1, "#include \"MatrixType.hpp\""});
    }
    
    lines.push_back({0, 0, 2, "#include \"SubMatrix.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_hpp_prim_includes(      std::ofstream& fstream,
                                          const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include <cstdint>"});
    
    lines.push_back({0, 0, 1, "#include \"SimdTypes.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"Point.hpp\""});
        
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_cpp_includes(      std::ofstream& fstream,
                                     const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + _file_name(integral) +  ".hpp\""});
    
    if ((integral[0] > 1) || (integral[1] > 1))
    {
        lines.push_back({0, 0, 2, "#include <cmath>"});
    }
    
    lines.push_back({0, 0, 1, "#include \"BatchFunc.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"T2CDistributor.hpp\""});
    
    _add_prim_call_includes(lines, integral);

    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_add_prim_call_includes(      VCodeLines&  lines,
                                         const I2CIntegral& integral) const
{
    if (integral.is_simple_integrand() && integral.is_simple())
    {
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            lines.push_back({0, 0, 2, "#include \"" + t2c::prim_file_name(integral) + ".hpp\""});
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                for (const auto& bcomp: Tensor(integral[0]).components())
                {
                    lines.push_back({0, 0, 1, "#include \"" + t2c::prim_file_name(bcomp, integral, true) + ".hpp\""});
                }
            }
            else
            {
                for (const auto& kcomp: Tensor(integral[1]).components())
                {
                    lines.push_back({0, 0, 1, "#include \"" + t2c::prim_file_name(kcomp, integral, false) + ".hpp\""});
                }
            }
            
            lines.push_back({0, 0, 1, ""});
        }
    }
    else
    {
        for (const auto& bcomp: Tensor(integral[0]).components())
        {
            for (const auto& kcomp: Tensor(integral[1]).components())
            {
                lines.push_back({0, 0, 1, "#include \"" + t2c::prim_file_name(bcomp, kcomp, integral) + ".hpp\""});
            }
        }
        
        lines.push_back({0, 0, 1, ""});
    }
}

void
T2CCPUGenerator::_write_cpp_prim_includes(      std::ofstream& fstream,
                                          const I2CIntegral&   integral) const
{
    auto fname = t2c::prim_file_name(integral) + ".hpp";
    
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + fname + "\""});
    
    lines.push_back({0, 0, 2, "#include <cmath>"});
    
    if (t2c::need_boys(integral))
    {
        lines.push_back({0, 0, 1, "#include \"BoysFunc.hpp\""});
    }
    
    lines.push_back({0, 0, 2, "#include \"MathConst.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_cpp_prim_includes(       std::ofstream&  fstream,
                                          const TensorComponent& component,
                                          const I2CIntegral&     integral,
                                          const bool             bra_first) const
{
    auto fname = t2c::prim_file_name(component, integral, bra_first) + ".hpp";
    
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + fname + "\""});
    
    lines.push_back({0, 0, 2, "#include <cmath>"});
    
    if (t2c::need_boys(integral))
    {
        lines.push_back({0, 0, 1, "#include \"BoysFunc.hpp\""});
    }
    
    lines.push_back({0, 0, 2, "#include \"MathConst.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_cpp_prim_includes(      std::ofstream&   fstream,
                                          const TensorComponent& bra_component,
                                          const TensorComponent& ket_component,
                                          const I2CIntegral&     integral) const
{
    auto fname = t2c::prim_file_name(bra_component, ket_component, integral) + ".hpp";
    
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + fname + "\""});
    
    lines.push_back({0, 0, 2, "#include <cmath>"});
    
    if (t2c::need_boys(integral))
    {
        lines.push_back({0, 0, 1, "#include \"BoysFunc.hpp\""});
    }
    
    lines.push_back({0, 0, 2, "#include \"MathConst.hpp\""});
    
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
