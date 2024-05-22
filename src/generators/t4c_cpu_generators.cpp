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

#include "t4c_cpu_generators.hpp"

#include <iostream>

#include "string_formater.hpp"
#include "file_stream.hpp"

#include "t4c_utils.hpp"
#include "t4c_docs.hpp"
#include "t4c_decl.hpp"
#include "t4c_body.hpp"
#include "t4c_prim_docs.hpp"
#include "t4c_prim_decl.hpp"
#include "t4c_prim_body.hpp"
#include "t4c_hrr_docs.hpp"
#include "t4c_hrr_decl.hpp"
#include "t4c_hrr_body.hpp"

#include "v4i_eri_driver.hpp"

void
T4CCPUGenerator::generate(const std::string& label,
                          const int          max_ang_mom) const
{
    if (_is_available(label))
    {
        for (int i = 0; i <= max_ang_mom; i++)
        {
            for (int j = i; j <= max_ang_mom; j++)
            {
                for (int k = 0; k <= max_ang_mom; k++)
                {
                    for (int l = k; l <= max_ang_mom; l++)
                    {
                        const auto integral = _get_integral(label, {i, j, k, l});
                        
                        const auto bra_integrals = _generate_bra_hrr_integral_group(integral);
                        
                        const auto ket_integrals = _generate_ket_hrr_integral_group(integral, bra_integrals);
                        
                        auto hrr_integrals = bra_integrals;
                        
                        hrr_integrals.insert(ket_integrals.begin(), ket_integrals.end());
                        
                        const auto vrr_integrals = _generate_vrr_integral_group(integral, hrr_integrals);
                        
                        _write_cpp_header(bra_integrals, ket_integrals, vrr_integrals, integral);
                    }
                }
            }
        }
        
//        for (int i = 0; i <= 2 * max_ang_mom; i++)
//        {
//            for (int j = 0; j <= 2 * max_ang_mom; j++)
//            {
//                const auto integral = _get_integral(label, {0, i, 0, j});
//                
//                _write_prim_cpp_header(integral);
//                
//                _write_prim_cpp_file(integral);
//            }
//        }
//        
//        for (int i = 1; i <= max_ang_mom; i++)
//        {
//            for (int j = i; j <= (2 * max_ang_mom - i) ; j++)
//            {
//                const auto integral = _get_integral(label, {0, 0, i, j});
//                
//                _write_ket_hrr_cpp_header(integral);
//                
//                _write_ket_hrr_cpp_file(integral);
//            }
//        }
//        
//        for (int i = 1; i <= max_ang_mom; i++)
//        {
//            for (int j = i; j <= (2 * max_ang_mom - i) ; j++)
//            {
//                const auto integral = _get_integral(label, {i, j, 0, 0});
//                
//                _write_bra_hrr_cpp_header(integral);
//                
//                _write_bra_hrr_cpp_file(integral);
//            }
//        }
    }
    else
    {
        std::cerr << "*** ERROR *** Unsupported type of four-center integral: ";
        
        std::cerr << label << " !!!" << std::endl;
        
        std::exit(EXIT_FAILURE);
    }
}

bool
T4CCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "electron repulsion") return true;
    
    return false;
}

I4CIntegral
T4CCPUGenerator::_get_integral(const std::string&        label,
                               const std::array<int, 4>& ang_moms) const
{
    // bra and ket sides
    
    const auto bpair = I2CPair("GA", ang_moms[0], "GB", ang_moms[1]);
    
    const auto kpair = I2CPair("GC", ang_moms[2], "GD", ang_moms[3]);
    
    // electron repulsion integrals
    
    if (fstr::lowercase(label) == "electron repulsion")
    {
        return I4CIntegral(bpair, kpair, Operator("1/|r-r'|"));
    }
    
    return I4CIntegral();
}

SI4CIntegrals
T4CCPUGenerator::_generate_bra_hrr_integral_group(const I4CIntegral& integral) const
{
    SI4CIntegrals tints;
    
    // Electron repulsion integrals
    
    if (integral.integrand() == Operator("1/|r-r'|"))
    {
        V4IElectronRepulsionDriver eri_drv;
        
        if (integral.is_simple())
        {
            tints = eri_drv.create_bra_hrr_recursion({integral,});
        }
        else
        {
            /// TODO: ...
        }
    }
    
    return tints;
}

SI4CIntegrals
T4CCPUGenerator::_generate_ket_hrr_integral_group(const I4CIntegral&   integral,
                                                  const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    // Electron repulsion integrals
    
    if (integral.integrand() == Operator("1/|r-r'|"))
    {
        V4IElectronRepulsionDriver eri_drv;
        
        for (const auto& tint : integrals)
        {
            if ((tint[0] == 0) && (tint[2] > 0))
            {
                const auto ctints = eri_drv.create_ket_hrr_recursion({tint, });
                
                tints.insert(ctints.cbegin(), ctints.cend());
            }
        }
    }
    
    return tints;
}

SI4CIntegrals
T4CCPUGenerator::_generate_vrr_integral_group(const I4CIntegral&   integral,
                                              const SI4CIntegrals& integrals) const
{
    SI4CIntegrals tints;
    
    // Electron repulsion integrals
    
    if (integral.integrand() == Operator("1/|r-r'|"))
    {
        V4IElectronRepulsionDriver eri_drv;
        
        for (const auto& tint : integrals)
        {
            if ((tint[0] == 0) && (tint[2] == 0))
            {
                const auto ctints = eri_drv.create_vrr_recursion({tint, });
                
                tints.insert(ctints.cbegin(), ctints.cend());
            }
        }
    }
    
    return tints;
}

std::string
T4CCPUGenerator::_file_name(const I4CIntegral& integral) const
{
    std::string label = "Rec" + integral.label();
    
    return t4c::integral_label(integral) + label;
}

void
T4CCPUGenerator::_write_cpp_header(const SI4CIntegrals& bra_integrals,
                                   const SI4CIntegrals& ket_integrals,
                                   const SI4CIntegrals& vrr_integrals,
                                   const I4CIntegral&   integral) const
{
    auto fname = _file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, true);
    
    _write_hpp_includes(fstream, bra_integrals, ket_integrals, vrr_integrals, integral);
    
    _write_namespace(fstream, integral, true);
    
    T4CDocuDriver docs_drv;
    
    T4CDeclDriver decl_drv;
    
    T4CFuncBodyDriver func_drv;

    if ((integral[0] == integral[2]) && (integral[1] == integral[3]))
    {
        docs_drv.write_doc_str(fstream, integral, true);
        
        decl_drv.write_func_decl(fstream, integral, true, false);
        
        func_drv.write_func_body(fstream, bra_integrals, ket_integrals, vrr_integrals, integral, true);

        fstream << std::endl;
    }

    docs_drv.write_doc_str(fstream, integral, false);
    
    decl_drv.write_func_decl(fstream, integral, false, false);
    
    func_drv.write_func_body(fstream, bra_integrals, ket_integrals, vrr_integrals, integral, false);
    
    fstream << std::endl;

    _write_namespace(fstream, integral, false);
        
    _write_hpp_defines(fstream, integral, false);
    
    fstream.close();
}

void
T4CCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                    const I4CIntegral&   integral,
                                    const bool           start) const
{
    auto fname = _file_name(integral) + "_hpp";
    
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
T4CCPUGenerator::_write_hpp_includes(      std::ofstream& fstream,
                                     const SI4CIntegrals& bra_integrals,
                                     const SI4CIntegrals& ket_integrals,
                                     const SI4CIntegrals& vrr_integrals,
                                     const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include <array>"});
    
    std::set<std::string> labels;
    
    for (const auto& tint : vrr_integrals)
    {
        if ((tint[0] + tint[2]) == 0)
        {
            labels.insert(t4c::prim_file_name(tint));
        }
    }
    
    for (const auto& tint : ket_integrals)
    {
        if ((tint[0] == 0) && (tint[2] > 0))
        {
            labels.insert(t4c::ket_hrr_file_name(tint));
        }
    }
    
    for (const auto& tint : bra_integrals)
    {
        if ((tint[0] > 0) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
        {
            labels.insert(t4c::bra_hrr_file_name(tint));
        }
    }
    
    for (const auto& label : labels)
    {
        lines.push_back({0, 0, 1, "#include \"" + label + ".hpp\""});
    }
    
    lines.push_back({0, 0, 1, "#include \"SimdArray.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"BoysFunc.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T4CUtils.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T2CUtils.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"GtoPairBlock.hpp\""});
   
    
        
    ost::write_code_lines(fstream, lines);
}

void
T4CCPUGenerator::_write_namespace(      std::ofstream& fstream,
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
T4CCPUGenerator::_write_cpp_file(const SI4CIntegrals& bra_integrals,
                                 const SI4CIntegrals& ket_integrals,
                                 const SI4CIntegrals& vrr_integrals,
                                 const I4CIntegral&   integral) const
{
    auto fname = _file_name(integral) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_cpp_includes(fstream, bra_integrals, ket_integrals, vrr_integrals, integral);

    _write_namespace(fstream, integral, true);

    T4CDeclDriver decl_drv;

    T4CFuncBodyDriver func_drv;

    if ((integral[0] == integral[2]) && (integral[1] == integral[3]))
    {
        decl_drv.write_func_decl(fstream, integral, true, false);

        func_drv.write_func_body(fstream, bra_integrals, ket_integrals, vrr_integrals, integral, true);
        
        fstream << std::endl;
    }

    decl_drv.write_func_decl(fstream, integral, false, false);

    func_drv.write_func_body(fstream, bra_integrals, ket_integrals, vrr_integrals, integral, false);

    fstream << std::endl;
    
    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
T4CCPUGenerator::_write_cpp_includes(      std::ofstream& fstream,
                                     const SI4CIntegrals& bra_integrals,
                                     const SI4CIntegrals& ket_integrals,
                                     const SI4CIntegrals& vrr_integrals,
                                     const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + _file_name(integral) +  ".hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"SimdArray.hpp\""});
    
    std::set<std::string> labels;
    
    for (const auto& tint : vrr_integrals)
    {
        if ((tint[0] + tint[2]) == 0)
        {
            labels.insert(t4c::prim_file_name(tint));
        }
    }
    
    for (const auto& tint : ket_integrals)
    {
        if ((tint[0] == 0) && (tint[2] > 0))
        {
            labels.insert(t4c::ket_hrr_file_name(tint));
        }
    }
    
    for (const auto& tint : bra_integrals)
    {
        if ((tint[0] > 0) && (tint[2] == integral[2]) && (tint[3] == integral[3]))
        {
            labels.insert(t4c::bra_hrr_file_name(tint));
        }
    }
    
    for (const auto& label : labels)
    {
        lines.push_back({0, 0, 1, "#include \"" + label + ".hpp\""});
    }
    
    lines.push_back({0, 0, 1, "#include \"BoysFunc.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T4CUtils.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"T2CUtils.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T4CCPUGenerator::_write_prim_cpp_header(const I4CIntegral& integral) const
{
    auto fname = t4c::prim_file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_prim_hpp_defines(fstream, integral, true);
    
    _write_prim_hpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T4CPrimDocuDriver docs_drv;

    docs_drv.write_doc_str(fstream, integral);

    T4CPrimDeclDriver decl_drv;

    decl_drv.write_func_decl(fstream, integral, true);

    _write_namespace(fstream, integral, false);
    
    _write_prim_hpp_defines(fstream, integral, false);
    
    fstream.close();
}

void
T4CCPUGenerator::_write_prim_hpp_defines(      std::ofstream& fstream,
                                         const I4CIntegral&   integral,
                                         const bool           start) const
{
    auto fname = t4c::prim_file_name(integral) + "_hpp";
    
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
T4CCPUGenerator::_write_prim_hpp_includes(      std::ofstream& fstream,
                                          const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"SimdArray.hpp\""});
        
    ost::write_code_lines(fstream, lines);
}

void
T4CCPUGenerator::_write_prim_cpp_file(const I4CIntegral& integral) const
{
    auto fname = t4c::prim_file_name(integral) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_prim_cpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T4CPrimDeclDriver decl_drv;
    
    decl_drv.write_func_decl(fstream, integral, false);

    T4CPrimFuncBodyDriver func_drv;

    func_drv.write_func_body(fstream, integral);
    
    fstream << std::endl;
    
    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
T4CCPUGenerator::_write_prim_cpp_includes(      std::ofstream& fstream,
                                          const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + t4c::prim_file_name(integral) +  ".hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T4CCPUGenerator::_write_ket_hrr_cpp_header(const I4CIntegral& integral) const
{
    auto fname = t4c::ket_hrr_file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_ket_hrr_hpp_defines(fstream, integral, true);
    
    _write_ket_hrr_hpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T4CHrrDocuDriver docs_drv;

    docs_drv.write_ket_doc_str(fstream, integral);

    T4CHrrDeclDriver decl_drv;

    decl_drv.write_ket_func_decl(fstream, integral, true);

    _write_namespace(fstream, integral, false);
    
    _write_ket_hrr_hpp_defines(fstream, integral, false);
    
    fstream.close();
}

void
T4CCPUGenerator::_write_ket_hrr_hpp_defines(      std::ofstream& fstream,
                                            const I4CIntegral&   integral,
                                            const bool           start) const
{
    auto fname = t4c::ket_hrr_file_name(integral) + "_hpp";
    
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
T4CCPUGenerator::_write_ket_hrr_hpp_includes(      std::ofstream& fstream,
                                             const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"SimdArray.hpp\""});
        
    ost::write_code_lines(fstream, lines);
}

void
T4CCPUGenerator::_write_ket_hrr_cpp_file(const I4CIntegral& integral) const
{
    auto fname = t4c::ket_hrr_file_name(integral) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_ket_hrr_cpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T4CHrrDeclDriver decl_drv;
    
    decl_drv.write_ket_func_decl(fstream, integral, false);

    T4CHrrFuncBodyDriver func_drv;

    func_drv.write_ket_func_body(fstream, integral);
    
    fstream << std::endl;
    
    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
T4CCPUGenerator::_write_ket_hrr_cpp_includes(      std::ofstream& fstream,
                                          const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + t4c::ket_hrr_file_name(integral) +  ".hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"TensorComponents.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T4CCPUGenerator::_write_bra_hrr_cpp_header(const I4CIntegral& integral) const
{
    auto fname = t4c::bra_hrr_file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_bra_hrr_hpp_defines(fstream, integral, true);
    
    _write_bra_hrr_hpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T4CHrrDocuDriver docs_drv;

    docs_drv.write_bra_doc_str(fstream, integral);

    T4CHrrDeclDriver decl_drv;

    decl_drv.write_bra_func_decl(fstream, integral, true);

    _write_namespace(fstream, integral, false);
    
    _write_bra_hrr_hpp_defines(fstream, integral, false);
    
    fstream.close();
}

void
T4CCPUGenerator::_write_bra_hrr_hpp_defines(      std::ofstream& fstream,
                                            const I4CIntegral&   integral,
                                            const bool           start) const
{
    auto fname = t4c::bra_hrr_file_name(integral) + "_hpp";
    
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
T4CCPUGenerator::_write_bra_hrr_hpp_includes(      std::ofstream& fstream,
                                             const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"SimdArray.hpp\""});
        
    ost::write_code_lines(fstream, lines);
}

void
T4CCPUGenerator::_write_bra_hrr_cpp_file(const I4CIntegral& integral) const
{
    auto fname = t4c::bra_hrr_file_name(integral) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_bra_hrr_cpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T4CHrrDeclDriver decl_drv;
    
    decl_drv.write_bra_func_decl(fstream, integral, false);

    T4CHrrFuncBodyDriver func_drv;

    func_drv.write_bra_func_body(fstream, integral);
    
    fstream << std::endl;
    
    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
T4CCPUGenerator::_write_bra_hrr_cpp_includes(      std::ofstream& fstream,
                                          const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + t4c::bra_hrr_file_name(integral) +  ".hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"TensorComponents.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}
