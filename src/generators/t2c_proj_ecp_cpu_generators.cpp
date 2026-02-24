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

#include "t2c_proj_ecp_cpu_generators.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "file_stream.hpp"
#include "t2c_utils.hpp"
#include "v2i_proj_ecp_driver.hpp"
#include "t2c_prim_docs.hpp"
#include "t2c_prim_decl.hpp"
#include "t2c_proj_prim_body.hpp"
#include "t2c_docs.hpp"
#include "t2c_decl.hpp"
#include "t2c_proj_ecp_body.hpp"

void
T2CProjECPCPUGenerator::generate(const std::string& label,
                                 const int          max_ang_mom,
                                 const int          proj_ang_mom) const
{
    if (_is_available(label))
    {
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                for (int l = 0; l <= proj_ang_mom; l++)
                {
                    for (int i = 0; i <= max_ang_mom; i++)
                    {
                        for (int j = 0; j <= max_ang_mom; j++)
                        {
                            #pragma omp task firstprivate(i,j,l)
                            {
                                const auto integral = _get_integral(label, {i, j}, l);
                                
                                const auto integrals = _generate_integral_group(integral);
                                
                                std::cout << " *** " << integral.second.label() << "_" << integral.second.order() << " *** " << std::endl;
                                
                                for (const auto& [order, tint] : integrals)
                                {
                                    std::cout << "> " << tint.label() << "_" << tint.order() << " : " ;
                                    
                                    std::cout << "(" << order[0] << ",";
                                    
                                    std::cout << order[1] << ","  << order[2] << ")" << std::endl;
                                }
                                
                                if ((i + j) > 0)
                                {
                                    _write_cpp_header(integrals, integral);
                                    
                                    _write_prim_cpp_header(integral);
                                    
                                    _write_prim_cpp_file(integral);
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
        std::cerr << "*** ERROR *** Unsupported type of two-center ECP integral: ";
        
        std::cerr << label << " !!!" << std::endl;
        
        std::exit(EXIT_FAILURE);
    }
}

bool
T2CProjECPCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "projected") return true;
    
    return false;
}

M2Integral
T2CProjECPCPUGenerator::_get_integral(const std::string&        label,
                                      const std::array<int, 2>& ang_moms,
                                      const int                 proj_ang_mom) const
{
    // bra and ket sides
    
    const auto bra = I1CPair("GA", ang_moms[0]);
    
    const auto ket = I1CPair("GB", ang_moms[1]);
    
    // projected core potential
    
    if (fstr::lowercase(label) == "projected")
    {
        M2Integral tint = {{0, 0, 0}, I2CIntegral(bra, ket, Operator("U_l"), proj_ang_mom, {})};
        
        return tint;
    }
    
    M2Integral tint = {{0, 0, 0}, I2CIntegral()};
    
    return tint;
}

SM2Integrals
T2CProjECPCPUGenerator::_generate_integral_group(const M2Integral& integral) const
{
    SM2Integrals tints;

    // Projected potential integrals
    
    if (integral.second.integrand() == Operator("U_l"))
    {
        V2IProjectedECPDriver ecp_drv;
        
        if (integral.second.is_simple())
        {
            tints = ecp_drv.create_recursion({integral,});
        }
        else
        {
            tints = ecp_drv.create_recursion(tints);
        }
    }
    
    return tints;
}

void
T2CProjECPCPUGenerator::_write_cpp_header(const SM2Integrals& integrals,
                                          const M2Integral&   integral) const
{
    auto fname = _file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, false, true);
    
    _write_hpp_includes(fstream, integrals, integral);
    
    _write_namespace(fstream, integral, true);
   
    T2CDocuDriver docs_drv;

    T2CDeclDriver decl_drv;
    
    T2CProjECPFuncBodyDriver func_drv;

    docs_drv.write_proj_ecp_doc_str(fstream, integral);

    decl_drv.write_proj_ecp_func_decl(fstream, integral, false);

    func_drv.write_func_body(fstream, integrals, integral);
    
    fstream << std::endl;

    _write_namespace(fstream, integral, false);
        
    _write_hpp_defines(fstream, integral, false, false);
    
    fstream.close();
}

void
T2CProjECPCPUGenerator::_write_hpp_includes(      std::ofstream& fstream,
                                            const SM2Integrals&  integrals,
                                            const M2Integral&    integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "#include <cstddef>"});
    
    lines.push_back({0, 0, 1, "#include <array>"});
    
    lines.push_back({0, 0, 1, "#include <vector>"});

    lines.push_back({0, 0, 2, "#include <utility>"});
    
    lines.push_back({0, 0, 1, "#include \"GtoBlock.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"BaseCorePotential.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"SimdArray.hpp\""});
    
    std::set<std::string> plabels;
    
    plabels.insert("ProjectedCorePotentialPrimRecSS");
    
    for (const auto& tint : integrals)
    {
        if ((tint.second[0] + tint.second[1]) > 0)
        {
            plabels.insert(t2c::prim_file_name(tint));
        }
    }
    
    for (const auto& plabel : plabels)
    {
        lines.push_back({0, 0, 1, "#include \"" + plabel + ".hpp\""});
    }
    
    lines.push_back({0, 0, 1, "#include \"T2CUtils.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T2CTransform.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"BatchFunc.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CProjECPCPUGenerator::_write_prim_cpp_header(const M2Integral& integral) const
{
    auto fname = t2c::prim_file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, true, true);
    
    _write_prim_hpp_includes(fstream, integral);
    
    _write_namespace(fstream, integral, true);
    
    T2CPrimDocuDriver docs_drv;
    
    docs_drv.write_doc_str(fstream, integral);
    
    T2CPrimDeclDriver decl_drv;
    
    decl_drv.write_func_decl(fstream, integral, true);
    
    _write_namespace(fstream, integral, false);
    
    _write_hpp_defines(fstream, integral, true, false);
    
    fstream.close();
}

void
T2CProjECPCPUGenerator::_write_prim_hpp_includes(      std::ofstream& fstream,
                                                 const M2Integral&    integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "#include \"SimdArray.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"Point.hpp\""});
        
    ost::write_code_lines(fstream, lines);
}

void
T2CProjECPCPUGenerator::_write_prim_cpp_file(const M2Integral& integral) const
{
    auto fname = t2c::prim_file_name(integral) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_prim_cpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T2CPrimDeclDriver decl_drv;
    
    decl_drv.write_func_decl(fstream, integral, false);

    T2CProjECPPrimFuncBodyDriver func_drv;

    func_drv.write_func_body(fstream, integral);
    
    fstream << std::endl;
    
    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
T2CProjECPCPUGenerator::_write_prim_cpp_includes(      std::ofstream& fstream,
                                                 const M2Integral&    integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + t2c::prim_file_name(integral) +  ".hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CProjECPCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                           const M2Integral&    integral,
                                           const bool           is_prim_rec,
                                           const bool           start) const
{
    auto fname = (is_prim_rec) ? t2c::prim_file_name(integral) : _file_name(integral) + "_hpp";
    
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

std::string
T2CProjECPCPUGenerator::_file_name(const M2Integral& integral) const
{
    auto label = integral.second.label() + "For" + Tensor(integral.second.order()).label();
    
    return t2c::integral_label(integral.second) + label;
}

void
T2CProjECPCPUGenerator::_write_namespace(      std::ofstream& fstream,
                                         const M2Integral&    integral,
                                         const bool           start) const
{
    const auto label = t2c::namespace_label(integral.second);
    
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
