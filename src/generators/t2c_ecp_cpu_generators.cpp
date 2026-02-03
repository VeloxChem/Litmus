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

#include "t2c_ecp_cpu_generators.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "file_stream.hpp"
#include "t2c_utils.hpp"
#include "t2c_docs.hpp"
#include "t2c_decl.hpp"
#include "t2c_ecp_body.hpp"

void
T2CECPCPUGenerator::generate(const std::string& label,
                             const int          max_ang_mom) const
{
    if (_is_available(label))
    {
        #pragma omp parallel
        {
            #pragma omp single nowait
            {
                for (int i = 0; i <= max_ang_mom; i++)
                {
                    for (int j = 0; j <= max_ang_mom; j++)
                    {
                        #pragma omp task firstprivate(i,j)
                        {
                            const auto integral = _get_integral(label, {i, j});
                            
                            const auto integrals = _generate_integral_group(integral);
                            
                            const auto hrr_integrals = _filter_hrr_integrals(integrals);
                            
                            const auto vrr_integrals = _filter_vrr_integrals(integrals);
                            
                            _write_cpp_header(hrr_integrals, vrr_integrals, integral);
                            
                            std::cout << "*** BASE INTEGRAL : " << integral.label() << " *** " << std::endl;
                            
                            std::cout << "HRR Integrals : " << std::endl;
                           
                            for (const auto& tint : hrr_integrals)
                            {
                                std::cout << " < > " << tint.label() << std::endl;
                            }
                            
                            std::cout << "VRR Integrals : " << std::endl;
                            
                            for (const auto& tint : vrr_integrals)
                            {
                                std::cout << " < > " << tint.label() << std::endl;
                            }
                            
//
//                            _write_cpp_header(integrals, integral, rec_form, use_rs);
//                            
//                            if (((i + j) >= 0) && (!use_rs))
//                            {
//                                _write_prim_cpp_header(integral, rec_form);
//                                    
//                                _write_prim_cpp_file(integral);
//                            }
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
T2CECPCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "local") return true;

    if (fstr::lowercase(label) == "projected") return true;
    
    return false;
}

I2CIntegral
T2CECPCPUGenerator::_get_integral(const std::string&        label,
                                  const std::array<int, 2>& ang_moms) const
{
    // bra and ket sides
    
    const auto bra = I1CPair("GA", ang_moms[0]);
    
    const auto ket = I1CPair("GB", ang_moms[1]);
    
    // local core potential
    
    if (fstr::lowercase(label) == "local")
    {
        return I2CIntegral(bra, ket, Operator("U_L"), 0, {});
    }
    
    // projected core potential
    
    if (fstr::lowercase(label) == "projected")
    {
        return I2CIntegral(bra, ket, Operator("U_l"), 0, {});
    }
    
    return I2CIntegral();
}

SI2CIntegrals
T2CECPCPUGenerator::_generate_integral_group(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;

    // Local core potential integrals
    
    if (integral.integrand() == Operator("U_L"))
    {
        // set up VRR recursions
        
        if (integral[0] > integral[1])
        {
            for (int i = 0; i <= integral[0] + integral[1]; i++)
            {
                tints.insert(_get_integral("local", {i, 0}));
            }
        } else {
            for (int i = 0; i <= integral[0] + integral[1]; i++)
            {
                tints.insert(_get_integral("local", {0, i}));
            }
        }
        
        // set up HRR recursions
        
        if ((integral[0] > 0) && (integral[1] > 0))
        {
            if (integral[0] > integral[1])
            {
                for (int i = 1; i <= integral[1]; i++)
                {
                    for (int j = integral[0]; j <= integral[0] + integral[1] - i; j++)
                    {
                        tints.insert(_get_integral("local", {j, i}));
                    }
                }
            } else {
                for (int i = 1; i <= integral[0]; i++)
                {
                    for (int j = integral[1]; j <= integral[0] + integral[1] - i; j++)
                    {
                        tints.insert(_get_integral("local", {i, j}));
                    }
                }
            }
        }
    }
    
    // Projected core potential 

    return tints;
}

SI2CIntegrals
T2CECPCPUGenerator::_filter_hrr_integrals(const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints;
    
    for (const auto& tint : integrals)
    {
        if ((tint[0] > 0) && (tint[1] > 0))
        {
            tints.insert(tint);
        }
    }
    return tints;
}

SI2CIntegrals
T2CECPCPUGenerator::_filter_vrr_integrals(const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints;
    
    for (const auto& tint : integrals)
    {
        if ((tint[0] == 0) || (tint[1] == 0))
        {
            tints.insert(tint);
        }
    }
    
    return tints;
}

void
T2CECPCPUGenerator::_write_cpp_header(const SI2CIntegrals& hrr_integrals,
                                      const SI2CIntegrals& vrr_integrals,
                                      const I2CIntegral&   integral) const
{
    auto fname = _file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, false, true);
    
    _write_hpp_includes(fstream, hrr_integrals, vrr_integrals, integral);
    
    _write_namespace(fstream, integral, true);
   
    T2CDocuDriver docs_drv;
    
    T2CDeclDriver decl_drv;
    
    T2CECPFuncBodyDriver func_drv;

    docs_drv.write_ecp_doc_str(fstream, integral);
    
    decl_drv.write_ecp_func_decl(fstream, integral, false);
    
    func_drv.write_func_body(fstream, hrr_integrals, vrr_integrals, integral);
    
    fstream << std::endl;

    _write_namespace(fstream, integral, false);
        
    _write_hpp_defines(fstream, integral, false, false);
    
    fstream.close();
}

std::string
T2CECPCPUGenerator::_file_name(const I2CIntegral& integral) const
{
    auto label = integral.label();
    
    return t2c::integral_label(integral) + label;
}

void
T2CECPCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                       const I2CIntegral&   integral,
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

void
T2CECPCPUGenerator::_write_hpp_includes(      std::ofstream& fstream,
                                        const SI2CIntegrals& hrr_integrals,
                                        const SI2CIntegrals& vrr_integrals,
                                        const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "#include <cstddef>"});
    
    lines.push_back({0, 0, 1, "#include <array>"});
    
    lines.push_back({0, 0, 1, "#include <vector>"});

    lines.push_back({0, 0, 2, "#include <utility>"});
    
    lines.push_back({0, 0, 1, "#include \"GtoBlock.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"BaseCorePotential.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"SimdArray.hpp\""});
    
    SI2CIntegrals rints;
    
    for (const auto& tint : vrr_integrals)
    {
        lines.push_back({0, 0, 1, "#include \"" + t2c::prim_file_name(tint) + ".hpp\""});
    }
    
    // need hrr recursion fixing 
    
    for (const auto& tint : hrr_integrals)
    {
       lines.push_back({0, 0, 1, "#include \"" + t2c::hrr_file_name(tint) + ".hpp\""});
    }
    
    lines.push_back({0, 0, 1, "#include \"T2CUtils.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T2CTransform.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"BatchFunc.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CECPCPUGenerator::_write_namespace(      std::ofstream& fstream,
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
