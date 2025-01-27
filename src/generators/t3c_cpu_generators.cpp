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

#include "t3c_cpu_generators.hpp"

#include "string_formater.hpp"
#include "file_stream.hpp"

#include "v3i_eri_driver.hpp"
#include "t3c_utils.hpp"
#include "t3c_docs.hpp"
#include "t3c_decl.hpp"
#include "t3c_body.hpp"
#include "t3c_prim_docs.hpp"
#include "t3c_prim_decl.hpp"
#include "t3c_prim_body.hpp"
#include "t3c_hrr_docs.hpp"
#include "t3c_hrr_decl.hpp"
#include "t3c_hrr_body.hpp"

void
T3CCPUGenerator::generate(const std::string& label,
                          const int          max_ang_mom,
                          const int          max_aux_ang_mom) const
{
    if (_is_available(label))
    {
        for (int i = 0; i <= max_aux_ang_mom; i++)
        {
            for (int j = 0; j <= max_ang_mom; j++)
            {
                for (int k = j; k <= max_ang_mom; k++)
                {
                    const auto integral = _get_integral(label, {i, j, k});
                    
                    const auto hrr_integrals = _generate_ket_hrr_integral_group(integral);
                    
                    const auto vrr_integrals = _generate_vrr_integral_group(integral, hrr_integrals);
                    
                    _write_cpp_header(hrr_integrals, vrr_integrals, integral);
                }
            }
        }
        
        for (int i = 0; i <= max_aux_ang_mom; i++)
        {
            for (int j = 0; j <= 2 * max_ang_mom; j++)
            {
                if ((i + j) == 0) continue;
                
                const auto integral = _get_integral(label, {i, 0, j});

                _write_prim_cpp_header(integral);

                _write_prim_cpp_file(integral);
            }
        }
        
        for (int i = 1; i <= max_ang_mom; i++)
        {
            for (int j = 0; j <= (2 * max_ang_mom - i) ; j++)
            {
                const auto integral = _get_integral(label, {0, i, j});
                
                _write_hrr_cpp_header(integral);
                
                _write_hrr_cpp_file(integral);
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
T3CCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "electron repulsion") return true;
    
    return false;
}

I3CIntegral
T3CCPUGenerator::_get_integral(const std::string&        label,
                               const std::array<int, 3>& ang_moms) const
{
    // bra and ket sides
    
    const auto bpair = I1CPair("GA", ang_moms[0]);
    
    const auto kpair = I2CPair("GC", ang_moms[1], "GD", ang_moms[2]);
    
    // electron repulsion integrals
    
    if (fstr::lowercase(label) == "electron repulsion")
    {
        return I3CIntegral(bpair, kpair, Operator("1/|r-r'|"));
    }
    
    return I3CIntegral();
}

SI3CIntegrals
T3CCPUGenerator::_generate_ket_hrr_integral_group(const I3CIntegral& integral) const
{
    SI3CIntegrals tints;
    
    // Electron repulsion integrals
    
    if (integral.integrand() == Operator("1/|r-r'|"))
    {
        V3IElectronRepulsionDriver eri_drv;
        
        if (integral.is_simple())
        {
            tints = eri_drv.create_ket_hrr_recursion({integral,});
        }
        else
        {
            /// TODO: ...
        }
    }
    
    return tints;
}

SI3CIntegrals
T3CCPUGenerator::_generate_vrr_integral_group(const I3CIntegral&   integral,
                                              const SI3CIntegrals& integrals) const
{
    SI3CIntegrals tints;
    
    // Electron repulsion integrals
    
    if (integral.integrand() == Operator("1/|r-r'|"))
    {
        V3IElectronRepulsionDriver eri_drv;
        
        for (const auto& tint : integrals)
        {
            if (tint[1] == 0)
            {
                const auto ctints = eri_drv.create_vrr_recursion({tint, });
                
                tints.insert(ctints.cbegin(), ctints.cend());
            }
        }
    }
    
    return tints;
}

void
T3CCPUGenerator::_write_cpp_header(const SI3CIntegrals& hrr_integrals,
                                   const SI3CIntegrals& vrr_integrals,
                                   const I3CIntegral&   integral) const
{
    auto fname = _file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, true);
    
    _write_hpp_includes(fstream, hrr_integrals, vrr_integrals, integral);

    _write_namespace(fstream, integral, true);
    
    T3CDocuDriver docs_drv;
    
    T3CDeclDriver decl_drv;
    
    T3CFuncBodyDriver func_drv;

    docs_drv.write_doc_str(fstream, integral);
    
    decl_drv.write_func_decl(fstream, integral, false);
    
    func_drv.write_func_body(fstream, hrr_integrals, vrr_integrals, integral);
    
    fstream << std::endl;

    _write_namespace(fstream, integral, false);
        
    _write_hpp_defines(fstream, integral, false);
    
    fstream.close();
}

std::string
T3CCPUGenerator::_file_name(const I3CIntegral& integral) const
{
    std::string label = "Rec" + integral.label();
    
    return t3c::integral_label(integral) + label;
}

void
T3CCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                    const I3CIntegral&   integral,
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
T3CCPUGenerator::_write_hpp_includes(      std::ofstream& fstream,
                                     const SI3CIntegrals& hrr_integrals,
                                     const SI3CIntegrals& vrr_integrals,
                                     const I3CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "#include <array>"});
    
    lines.push_back({0, 0, 1, "#include <cstddef>"});
    
    lines.push_back({0, 0, 2, "#include <utility>"});
    
    std::set<std::string> labels;
    
    for (const auto& tint : vrr_integrals)
    {
        if (tint[1] == 0)
        {
            labels.insert(t3c::prim_file_name(tint));
        }
    }
    
    for (const auto& tint : hrr_integrals)
    {
        if (tint[1] != 0)
        {
            labels.insert(t3c::hrr_file_name(tint));
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
    
    lines.push_back({0, 0, 1, "#include \"GtoPairBlock.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"GtoBlock.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"BatchFunc.hpp\""});
        
    ost::write_code_lines(fstream, lines);
}

void
T3CCPUGenerator::_write_namespace(      std::ofstream& fstream,
                                  const I3CIntegral&   integral,
                                  const bool           start) const
{
    const auto label = t3c::namespace_label(integral);
    
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
T3CCPUGenerator::_write_prim_cpp_header(const I3CIntegral& integral) const
{
    auto fname = t3c::prim_file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_prim_hpp_defines(fstream, integral, true);
    
    _write_prim_hpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T3CPrimDocuDriver docs_drv;

    docs_drv.write_doc_str(fstream, integral);

    T3CPrimDeclDriver decl_drv;

    decl_drv.write_func_decl(fstream, integral, true);

    _write_namespace(fstream, integral, false);
    
    _write_prim_hpp_defines(fstream, integral, false);
    
    fstream.close();
}

void
T3CCPUGenerator::_write_prim_cpp_file(const I3CIntegral& integral) const
{
    auto fname = t3c::prim_file_name(integral) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_prim_cpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T3CPrimDeclDriver decl_drv;
    
    decl_drv.write_func_decl(fstream, integral, false);

    T3CPrimFuncBodyDriver func_drv;

    func_drv.write_func_body(fstream, integral);
    
    fstream << std::endl;
    
    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
T3CCPUGenerator::_write_prim_hpp_defines(      std::ofstream& fstream,
                                         const I3CIntegral&   integral,
                                         const bool           start) const
{
    auto fname = t3c::prim_file_name(integral) + "_hpp";
    
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
T3CCPUGenerator::_write_prim_hpp_includes(      std::ofstream& fstream,
                                          const I3CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include <cstddef>"});
    
    lines.push_back({0, 0, 1, "#include \"Point.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"SimdArray.hpp\""});
        
    ost::write_code_lines(fstream, lines);
}

void
T3CCPUGenerator::_write_prim_cpp_includes(      std::ofstream& fstream,
                                          const I3CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + t3c::prim_file_name(integral) +  ".hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T3CCPUGenerator::_write_hrr_cpp_header(const I3CIntegral& integral) const
{
    auto fname = t3c::hrr_file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hrr_hpp_defines(fstream, integral, true);
    
    _write_hrr_hpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T3CHrrDocuDriver docs_drv;

    docs_drv.write_doc_str(fstream, integral);

    T3CHrrDeclDriver decl_drv;

    decl_drv.write_func_decl(fstream, integral, true);

    _write_namespace(fstream, integral, false);
    
    _write_hrr_hpp_defines(fstream, integral, false);
    
    fstream.close();
}

void
T3CCPUGenerator::_write_hrr_cpp_file(const I3CIntegral& integral) const
{
    auto fname = t3c::hrr_file_name(integral) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_hrr_cpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T3CHrrDeclDriver decl_drv;
    
    decl_drv.write_func_decl(fstream, integral, false);

    T3CHrrFuncBodyDriver func_drv;

    func_drv.write_func_body(fstream, integral);
    
    fstream << std::endl;
    
    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
T3CCPUGenerator::_write_hrr_hpp_defines(      std::ofstream& fstream,
                                        const I3CIntegral&   integral,
                                        const bool           start) const
{
    auto fname = t3c::hrr_file_name(integral) + "_hpp";
    
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
T3CCPUGenerator::_write_hrr_hpp_includes(      std::ofstream& fstream,
                                         const I3CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include <cstddef>"});
    
    lines.push_back({0, 0, 2, "#include \"SimdArray.hpp\""});
        
    ost::write_code_lines(fstream, lines);
}

void
T3CCPUGenerator::_write_hrr_cpp_includes(      std::ofstream& fstream,
                                          const I3CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + t3c::hrr_file_name(integral) +  ".hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"TensorComponents.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}
