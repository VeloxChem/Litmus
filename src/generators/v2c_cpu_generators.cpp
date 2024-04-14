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

#include "v2c_cpu_generators.hpp"

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
V2CCPUGenerator::generate(const std::string& label,
                           const int         angmom,
                           const int         bra_gdrv,
                           const int         ket_gdrv,
                           const int         op_gdrv,
                           const bool        sum_form,
                           const bool        diag_form) const
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
//                            const auto rgroup = _generate_integral_group(integral);

                            _write_cpp_header(integral, sum_form, diag_form);

                            _write_cpp_file(integral, sum_form, diag_form);
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
V2CCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "overlap") return true;
    
    if (fstr::lowercase(label) == "kinetic energy") return true;
    
    if (fstr::lowercase(label) == "nuclear potential") return true;
    
    return false;
}

I2CIntegral
V2CCPUGenerator::_get_integral(const std::string& label,
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
V2CCPUGenerator::_file_name(const I2CIntegral& integral,
                            const bool         sum_form,
                            const bool         diag_form) const
{
    std::string label = "Rec" + integral.label();
    
    if (sum_form) label = "Sum" + label;
    
    if (diag_form) label = "Diag" + label;
    
    return t2c::integral_label(integral) + label;
}

void
V2CCPUGenerator::_write_cpp_header(const I2CIntegral& integral,
                                   const bool         sum_form,
                                   const bool         diag_form) const
{
    auto fname = _file_name(integral, sum_form, diag_form) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_hpp_defines(fstream, integral, sum_form, diag_form, true);
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
    _write_hpp_defines(fstream, integral, sum_form, diag_form, false);

    fstream.close();
}

void
V2CCPUGenerator::_write_cpp_file(const I2CIntegral& integral,
                                 const bool         sum_form,
                                 const bool         diag_form) const
{
    auto fname = _file_name(integral, sum_form, diag_form) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
//    _write_cpp_includes(fstream, integral, sum_form);
//
//    _write_namespace(fstream, integral, true);
//
//    T2CDeclDriver decl_drv;
//
//    C2CFuncBodyDriver func_drv;
//
//    if ((integral[0] == integral[1]) && (integral.is_simple()))
//    {
//        decl_drv.write_func_decl(fstream, integral, sum_form, true, false);
//
//        func_drv.write_func_body(fstream, rgroup, integral, sum_form, true);
//    }
//
//    decl_drv.write_func_decl(fstream, integral, sum_form, false, false);
//
//    func_drv.write_func_body(fstream, rgroup, integral, sum_form, false);
//
//    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
V2CCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                    const I2CIntegral&   integral,
                                    const bool           sum_form,
                                    const bool           diag_form,
                                    const bool           start) const
{
    auto fname = _file_name(integral, sum_form, diag_form) + "_hpp";
    
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
