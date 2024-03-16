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

#include "cold_cpu_generator.hpp"

#include <iostream>
#include <iterator>

#include "operator.hpp"
#include "string_formater.hpp"
#include "spherical_momentum.hpp"

#include "t2c_defs.hpp"
#include "t2c_utils.hpp"
#include "t2c_docs.hpp"
#include "t2c_decl.hpp"
#include "c2c_body.hpp"
#include "c2c_auxilary_body.hpp"

#include "cold_ovl_driver.hpp"
#include "cold_kin_driver.hpp"
#include "cold_npot_driver.hpp"

void
ColdCPUGenerator::generate(const std::string& label,
                           const int          angmom,
                           const int          bra_gdrv,
                           const int          ket_gdrv,
                           const int          op_gdrv,
                           const bool         sum_form) const
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
                        const auto integral = _get_integral(label, i, j, bra_gdrv, ket_gdrv, op_gdrv);
                        
                        #pragma omp task firstprivate(integral)
                        {
                            const auto rgroup = _generate_integral_group(integral);
                            
                            _write_cpp_header(integral, sum_form);
                            
                            _write_cpp_file(rgroup, integral, sum_form);
                            
                            _write_auxilary_header(rgroup, integral);
                            
                            _write_auxilary_file(rgroup, integral, sum_form);
                        }
                    }
                }
            }
        }
        
        //_write_func_header(label, angmom, bra_gdrv, ket_gdrv, op_gdrv, sum_form);
    }
    else
    {
        std::cerr << "*** ERROR *** Unsupported type of two-center integral: ";
        
        std::cerr << label << " !!!" << std::endl;
        
        std::exit(EXIT_FAILURE);
    }
}

bool
ColdCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "overlap") return true;
    
    if (fstr::lowercase(label) == "kinetic energy") return true;
    
    if (fstr::lowercase(label) == "nuclear potential") return true;
    
    return false;
}

I2CIntegral
ColdCPUGenerator::_get_integral(const std::string& label,
                                const int          ang_a,
                                const int          ang_b,
                                const int          bra_gdrv,
                                const int          ket_gdrv,
                                const int          op_gdrv) const
{
    // bra and ket sides
    
    const auto bra = I1CPair("GA", ang_a);
    
    const auto ket = I1CPair("GB", ang_b);
    
    // prefixes of integral bra, ket order
    
    VOperators prefixes;
    
    if (bra_gdrv > 0) prefixes.push_back(Operator("d/dR", Tensor(bra_gdrv)));
    
    if (ket_gdrv > 0) prefixes.push_back(Operator("d/dR", Tensor(ket_gdrv)));
    
    // overlap integrals
    
    if (fstr::lowercase(label) == "overlap")
    {
        return I2CIntegral(bra, ket, Operator("1"), 0, prefixes);
    }
    
    // kinetic energy integrals
    
    if (fstr::lowercase(label) == "kinetic energy")
    {
        return I2CIntegral(bra, ket, Operator("T"), 0, prefixes);
    }
    
    // nuclear potential integrals
    
    if (fstr::lowercase(label) == "nuclear potential")
    {
        return I2CIntegral(bra, ket, Operator("A"), 0, prefixes);
    }
    
    return I2CIntegral();
}

std::string
ColdCPUGenerator::_file_name(const I2CIntegral& integral,
                            const bool         sum_form) const
{
    if (sum_form)
    {
        return t2c::integral_label(integral) + "SumColdRec" + integral.label();
    }
    else
    {
        return t2c::integral_label(integral) + "ColdRec" + integral.label();
    }
}

std::string
ColdCPUGenerator::_func_file_name(const std::string& label,
                                  const int          bra_gdrv,
                                  const int          ket_gdrv,
                                  const int          op_gdrv,
                                  const bool         sum_form) const
{
    std::string fname = (sum_form) ? "SumFunc" : "Func";
    
    if ((bra_gdrv + ket_gdrv + op_gdrv) > 0)
    {
        if (_need_all_gdrv(label))
        {
            fname = "Geom" + std::to_string(bra_gdrv) + std::to_string(op_gdrv) + std::to_string(ket_gdrv) + fname;
        }
        else
        {
            fname = "Geom" + std::to_string(bra_gdrv) + std::to_string(ket_gdrv) + fname;
        }
    }
    
    if (fstr::lowercase(label) == "overlap")
    {
        fname = "Overlap" + fname;
    }
    
    return fname;
}

void
ColdCPUGenerator::_write_cpp_header(const I2CIntegral& integral,
                                    const bool         sum_form) const
{
    auto fname = _file_name(integral, sum_form) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_hpp_defines(fstream, integral, false, sum_form, true);

    _write_hpp_includes(fstream, integral, sum_form);

    _write_namespace(fstream, integral, true);
        
    T2CDocuDriver docs_drv;

    T2CDeclDriver decl_drv;

    if ((integral[0] == integral[1]) && integral.is_simple())
    {
        docs_drv.write_doc_str(fstream, integral, sum_form, true);

        decl_drv.write_func_decl(fstream, integral, sum_form, true, true);
    }

    docs_drv.write_doc_str(fstream, integral, sum_form, false);

    decl_drv.write_func_decl(fstream, integral, sum_form, false, true);

    _write_namespace(fstream, integral, false);

    _write_hpp_defines(fstream, integral, false, sum_form, false);

    fstream.close();
}

void
ColdCPUGenerator::_write_cpp_file(const R2Group&     rgroup,
                                  const I2CIntegral& integral,
                                  const bool         sum_form) const
{
    auto fname = _file_name(integral, sum_form) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_cpp_includes(fstream, integral, sum_form);
        
    _write_namespace(fstream, integral, true);
        
    T2CDeclDriver decl_drv;
        
    C2CFuncBodyDriver func_drv;
        
    if ((integral[0] == integral[1]) && (integral.is_simple()))
    {
        decl_drv.write_func_decl(fstream, integral, sum_form, true, false);
            
        func_drv.write_func_body(fstream, rgroup, integral, sum_form, true);
    }
        
    decl_drv.write_func_decl(fstream, integral, sum_form, false, false);
        
    func_drv.write_func_body(fstream, rgroup, integral, sum_form, false);
        
    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
ColdCPUGenerator::_write_auxilary_header(const R2Group&     rgroup,
                                         const I2CIntegral& integral) const
{
    auto fname = t2c::auxilary_file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_hpp_defines(fstream, integral, true, false, true);

    _write_auxilary_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T2CDocuDriver docs_drv;

    T2CDeclDriver decl_drv;

    if ((integral[0] == integral[1]) && integral.is_simple())
    {
        docs_drv.write_auxilary_doc_str(fstream, integral, true);

        decl_drv.write_auxilary_func_decl(fstream, rgroup, integral, true, true);
    }

    docs_drv.write_auxilary_doc_str(fstream, integral, false);

    decl_drv.write_auxilary_func_decl(fstream, rgroup, integral, false, true);

    _write_namespace(fstream, integral, false);

    _write_hpp_defines(fstream, integral, true, false, false);

    fstream.close();
}

void
ColdCPUGenerator::_write_auxilary_file(const R2Group&     rgroup,
                                       const I2CIntegral& integral,
                                       const bool         sum_form) const
{
    auto fname = t2c::auxilary_file_name(integral) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_auxilary_file_includes(fstream, integral);
        
    _write_namespace(fstream, integral, true);
        
    T2CDeclDriver decl_drv;
        
    C2CAuxilaryBodyDriver func_drv;
        
    if ((integral[0] == integral[1]) && (integral.is_simple()))
    {
        decl_drv.write_auxilary_func_decl(fstream, rgroup, integral, true, false);
            
        func_drv.write_aux_body(fstream, rgroup, integral, sum_form, true);
    }
        
    decl_drv.write_auxilary_func_decl(fstream, rgroup, integral, false, false);
        
    func_drv.write_aux_body(fstream, rgroup, integral, sum_form, false);
        
    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
ColdCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                     const I2CIntegral&   integral,
                                     const bool           is_auxilary,
                                     const bool           sum_form,
                                     const bool           start) const
{
    auto fname = _file_name(integral, sum_form) + "_hpp";
    
    if (is_auxilary)
    {
        fname = t2c::auxilary_file_name(integral) + "_hpp";
    }
    
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
ColdCPUGenerator::_write_namespace(      std::ofstream& fstream,
                                   const I2CIntegral&   integral,
                                   const bool           start) const
{
    const auto label = "cold_" + t2c::namespace_label(integral);
    
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
ColdCPUGenerator::_write_hpp_includes(      std::ofstream& fstream,
                                      const I2CIntegral&   integral,
                                      const bool           sum_form) const
{
    auto lines = VCodeLines();
    
    if (sum_form)
    {
        lines.push_back({0, 0, 2, "#include <vector>"});
    }
    
    lines.push_back({0, 0, 1, "#include \"GtoBlock.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"Matrix.hpp\""});

    if (integral.integrand().name() == "A")
    {
        lines.push_back({0, 0, 1, "#include \"TPoints.hpp\""});
    }
    
    if (integral.integrand().name() == "AG")
    {
        lines.push_back({0, 0, 1, "#include \"TPoints.hpp\""});
        
        if (integral.integrand().shape().order() > 1)
        {
            lines.push_back({0, 0, 1, "#include \"TensorTypes.hpp\""});
        }
    }
    
    lines.push_back({0, 0, 1, "#include \"TIndexes.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"T1EDistributor.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
ColdCPUGenerator::_write_auxilary_includes(      std::ofstream& fstream,
                                           const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
 
    lines.push_back({0, 0, 1, "#include \"GtoBlock.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"TIndexes.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"TSimd.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
ColdCPUGenerator::_write_cpp_includes(      std::ofstream& fstream,
                                      const I2CIntegral&   integral,
                                      const bool           sum_form) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + _file_name(integral, sum_form) +  ".hpp\""});
    
    if ((integral[0] > 1) || (integral[1] > 1))
    {
        lines.push_back({0, 0, 2, "#include <cmath>"});
    }
    
    lines.push_back({0, 0, 1, "#include \"BatchFunc.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"TSimd.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"TransformationFunc.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"" + t2c::auxilary_file_name(integral) + ".hpp\""});

    ost::write_code_lines(fstream, lines);
}

void
ColdCPUGenerator::_write_auxilary_file_includes(      std::ofstream& fstream,
                                                const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + t2c::auxilary_file_name(integral) +  ".hpp\""});
    
    lines.push_back({0, 0, 2, "#include <cmath>"});
    
    lines.push_back({0, 0, 2, "#include \"MathConst.hpp\""});

    ost::write_code_lines(fstream, lines);
}

R2Group
ColdCPUGenerator::_generate_integral_group(const I2CIntegral& integral) const
{
    R2Group rgroup;
    
    // Overlap integrals
    
    if (integral.integrand() == Operator("1"))
    {
        ColdOverlapDriver cold_ovl_drv;
        
        if (integral.is_simple())
        {
            rgroup = cold_ovl_drv.create_recursion(integral.components<T1CPair, T1CPair>());
        }
        else
        {
            cold_ovl_drv.apply_recursion(rgroup);
        }
    }
    
    // Kinetic energy integrals
    
    if (integral.integrand() == Operator("T"))
    {
        ColdKineticEnergyDriver cold_kin_drv;
        
        if (integral.is_simple())
        {
            rgroup = cold_kin_drv.create_recursion(integral.components<T1CPair, T1CPair>());
        }
        else
        {
            cold_kin_drv.apply_recursion(rgroup);
        }
    }
    
    // Nuclear potential integrals
    
    if (const auto integrand = integral.integrand();
        (integrand.name() == "A") && (integrand.shape() == Tensor(0)))
    {
        ColdNuclearPotentialDriver cold_npot_drv;
        
        if (integral.is_simple())
        {
            rgroup = cold_npot_drv.create_recursion(integral.components<T1CPair, T1CPair>());
        }
        else
        {
            cold_npot_drv.apply_recursion(rgroup);
        }
    }
    
    // ... other integrals
    
    rgroup.simplify();
    
    return rgroup;
}

void
ColdCPUGenerator::_write_func_header(const std::string& label,
                                     const int          angmom,
                                     const int          bra_gdrv,
                                     const int          ket_gdrv,
                                     const int          op_gdrv,
                                     const bool         sum_form) const
{
    auto fname = _func_file_name(label, bra_gdrv, ket_gdrv, op_gdrv, sum_form) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
//    _write_hpp_defines(fstream, integral, false, sum_form, true);
//
//    _write_hpp_includes(fstream, integral, sum_form);
//
//    _write_namespace(fstream, integral, true);
//
//    T2CDocuDriver docs_drv;
//
//    T2CDeclDriver decl_drv;
//
//    if ((integral[0] == integral[1]) && integral.is_simple())
//    {
//        docs_drv.write_doc_str(fstream, integral, sum_form, true);
//
//        decl_drv.write_func_decl(fstream, integral, sum_form, true, true);
//    }
//
//    docs_drv.write_doc_str(fstream, integral, sum_form, false);
//
//    decl_drv.write_func_decl(fstream, integral, sum_form, false, true);
//
//    _write_namespace(fstream, integral, false);
//
//    _write_hpp_defines(fstream, integral, false, sum_form, false);

    fstream.close();
}

bool
ColdCPUGenerator::_need_all_gdrv(const std::string& label) const
{
//    if (fstr::lowercase(label) == "overlap")
//    {
//        fname = "Overlap" + fname;
//    }
    
    return false;
}
