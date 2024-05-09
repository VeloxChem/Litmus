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

#include "t2c_cpu_generators.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "file_stream.hpp"

#include "t2c_defs.hpp"
#include "t2c_utils.hpp"
#include "t2c_docs.hpp"
#include "t2c_decl.hpp"
#include "t2c_body.hpp"
#include "t2c_prim_docs.hpp"
#include "t2c_prim_decl.hpp"
#include "t2c_prim_body.hpp"

#include "v2i_ovl_driver.hpp"
#include "v2i_kin_driver.hpp"
#include "v2i_dip_driver.hpp"
#include "v2i_npot_driver.hpp"

void
T2CCPUGenerator::generate(const std::string&           label,
                          const int                    max_ang_mom,
                          const std::array<int, 3>&    geom_drvs,
                          const std::pair<bool, bool>& rec_form) const
{
    if (_is_available(label))
    {
        SI2CIntegrals all_integrals;
        
        for (int i = 0; i <= max_ang_mom; i++)
        {
            for (int j = 0; j <= max_ang_mom; j++)
            {
                const auto integral = _get_integral(label, {i, j}, geom_drvs);
                
                const auto integrals = _generate_integral_group(integral);

                _write_cpp_header(integral, rec_form);
                
                _write_cpp_file(integrals, integral, rec_form);
                
                _write_prim_cpp_header(integral, rec_form);
                
                _write_prim_cpp_file(integral);
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

// MR: Need changes here
bool
T2CCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "overlap") return true;

    if (fstr::lowercase(label) == "dipole moment") return true;
    
    if (fstr::lowercase(label) == "kinetic energy") return true;
    
    if (fstr::lowercase(label) == "nuclear potential") return true;
        
    return false;
}

// MR: Needs changes for new cases
I2CIntegral
T2CCPUGenerator::_get_integral(const std::string&        label,
                               const std::array<int, 2>& ang_moms,
                               const std::array<int, 3>& geom_drvs) const
{
    // bra and ket sides
    
    const auto bra = I1CPair("GA", ang_moms[0]);
    
    const auto ket = I1CPair("GB", ang_moms[1]);
    
    // prefixes of integral bra, ket order
    
    VOperators prefixes;
    
    if (geom_drvs[0] > 0)
    {
        prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[0])));
    }
    
    if (geom_drvs[2] > 0)
    {
        prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[2])));
    }
    
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

    // dipole moment integrals

    if (fstr::lowercase(label) == "dipole moment")
    {
    /// Operator takes one or two arguments; first name, then optionally rank of tensor associated with operator
        return I2CIntegral(bra, ket, Operator("r", Tensor(1)), 0, prefixes);
    }
    
    // nuclear potential integrals
    
    if (fstr::lowercase(label) == "nuclear potential")
    {
        return I2CIntegral(bra, ket, Operator("A"), 0, prefixes);
    }
    
    return I2CIntegral();
}

// MR: Changes here for new integral cases
SI2CIntegrals
T2CCPUGenerator::_generate_integral_group(const I2CIntegral& integral) const
{
    SI2CIntegrals tints;
    
    // Overlap integrals
    
    if (integral.integrand() == Operator("1"))
    {
        V2IOverlapDriver ovl_drv;
        
        if (integral.is_simple())
        {
            tints = ovl_drv.create_recursion({integral,});
        }
        else
        {
            /// TODO: ...
        }
    }

    // Dipole moment integrals

    if (integral.integrand() == Operator("r", Tensor(1)))
    {
        V2IDipoleDriver dip_drv;

        tints = dip_drv.create_recursion({integral,});

        V2IOverlapDriver ovl_drv;

        tints = ovl_drv.create_recursion(tints);
    }

    // Kinetic energy integrals
    
    if (integral.integrand() == Operator("T"))
    {
        V2IKineticEnergyDriver kin_drv;
        
        if (integral.is_simple())
        {
            tints = kin_drv.create_recursion({integral,});
            
            V2IOverlapDriver ovl_drv;
            
            tints = ovl_drv.create_recursion(tints);
        }
        else
        {
            /// TODO: ...
        }
    }
    
    // Nuclear potential integrals
    
    if (integral.integrand() == Operator("A"))
    {
        V2INuclearPotentialDriver npot_drv;
        
        if (integral.is_simple())
        {
            tints = npot_drv.create_recursion({integral,});
        }
        else
        {
            /// TODO: ...
        }
    }
    
    return tints;
}

std::string
T2CCPUGenerator::_file_name(const I2CIntegral&           integral,
                            const std::pair<bool, bool>& rec_form) const
{
    std::string label = "Rec" + integral.label();
    
    if (rec_form.first) label = "Sum" + label;
    
    if (rec_form.second) label = "Conv" + label;
    
    return t2c::integral_label(integral) + label;
}

void
T2CCPUGenerator::_write_cpp_header(const I2CIntegral&           integral,
                                   const std::pair<bool, bool>& rec_form) const
{
    auto fname = _file_name(integral, rec_form) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, rec_form, false, true);
    
    _write_hpp_includes(fstream, integral, rec_form);
    
    _write_namespace(fstream, integral, true);
    
    T2CDocuDriver docs_drv;
    
    T2CDeclDriver decl_drv;

    if ((integral[0] == integral[1]) && integral.is_simple())
    {
        docs_drv.write_doc_str(fstream, integral, rec_form, true);
        
        decl_drv.write_func_decl(fstream, integral, rec_form, true, true);
        
        fstream << std::endl;
    }

    docs_drv.write_doc_str(fstream, integral, rec_form, false);
    
    decl_drv.write_func_decl(fstream, integral, rec_form, false, true);

    _write_namespace(fstream, integral, false);
        
    _write_hpp_defines(fstream, integral, rec_form, false, false);
    
    fstream.close();
}

void
T2CCPUGenerator::_write_hpp_defines(      std::ofstream&         fstream,
                                    const I2CIntegral&           integral,
                                    const std::pair<bool, bool>& rec_form,
                                    const bool                   is_prim_rec, 
                                    const bool                   start) const
{
    auto fname = (is_prim_rec) ? t2c::prim_file_name(integral) : _file_name(integral, rec_form) + "_hpp";
    
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
T2CCPUGenerator::_write_hpp_includes(      std::ofstream&         fstream,
                                     const I2CIntegral&           integral,
                                     const std::pair<bool, bool>& rec_form) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include <array>"});
    
    lines.push_back({0, 0, 1, "#include \"GtoBlock.hpp\""});

    if (integral[0] == integral[1])
    {
       lines.push_back({0, 0, 1, "#include \"Matrix.hpp\""});
    }
    
    lines.push_back({0, 0, 2, "#include \"SubMatrix.hpp\""});
        
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
T2CCPUGenerator::_write_cpp_file(const SI2CIntegrals&         integrals,
                                 const I2CIntegral&           integral,
                                 const std::pair<bool, bool>& rec_form) const
{
    auto fname = _file_name(integral, rec_form) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_cpp_includes(fstream, integrals, integral, rec_form);

    _write_namespace(fstream, integral, true);

    T2CDeclDriver decl_drv;

    T2CFuncBodyDriver func_drv;

    if ((integral[0] == integral[1]) && (integral.is_simple()))
    {
        decl_drv.write_func_decl(fstream, integral, rec_form, true, false);

        func_drv.write_func_body(fstream, integrals, integral, rec_form, true);
        
        fstream << std::endl; 
    }

    decl_drv.write_func_decl(fstream, integral, rec_form, false, false);

    func_drv.write_func_body(fstream, integrals, integral, rec_form, false);

    fstream << std::endl;
    
    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
T2CCPUGenerator::_write_cpp_includes(      std::ofstream&         fstream,
                                     const SI2CIntegrals&         integrals,
                                     const I2CIntegral&           integral,
                                     const std::pair<bool, bool>& rec_form) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + _file_name(integral, rec_form) +  ".hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"SimdArray.hpp\""});
    
    for (const auto& tint : integrals)
    {
        lines.push_back({0, 0, 1, "#include \"" + t2c::prim_file_name(tint) + ".hpp\""});
    }
    
    if (integral.integrand().name() == "A")
    {
        lines.push_back({0, 0, 1, "#include \"BoysFunc.hpp\""});
    }
    
    lines.push_back({0, 0, 1, "#include \"T2CDistributor.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"T2CUtils.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_prim_cpp_header(const I2CIntegral&           integral,
                                        const std::pair<bool, bool>& rec_form) const
{
    auto fname = t2c::prim_file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, rec_form, true, true);
    
    _write_prim_hpp_includes(fstream, integral);
    
    _write_namespace(fstream, integral, true);
    
    T2CPrimDocuDriver docs_drv;
    
    docs_drv.write_doc_str(fstream, integral);
    
    T2CPrimDeclDriver decl_drv;
    
    decl_drv.write_func_decl(fstream, integral, true);
    
    _write_namespace(fstream, integral, false);
    
    _write_hpp_defines(fstream, integral, rec_form, true, false);
    
    fstream.close();
}

void
T2CCPUGenerator::_write_prim_hpp_includes(      std::ofstream& fstream,
                                          const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"SimdArray.hpp\""});
        
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_prim_cpp_file(const I2CIntegral& integral) const
{
    auto fname = t2c::prim_file_name(integral) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_prim_cpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T2CPrimDeclDriver decl_drv;
    
    decl_drv.write_func_decl(fstream, integral, false);

    T2CPrimFuncBodyDriver func_drv;

    func_drv.write_func_body(fstream, integral);
    
    fstream << std::endl; 
    
    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
T2CCPUGenerator::_write_prim_cpp_includes(      std::ofstream& fstream,
                                          const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + t2c::prim_file_name(integral) +  ".hpp\""});
    
    ost::write_code_lines(fstream, lines);
}
