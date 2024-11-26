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

#include "t4c_geom_hrr_cpu_generators.hpp"

#include "string_formater.hpp"
#include "file_stream.hpp"

#include "t4c_utils.hpp"
#include "t4c_docs.hpp"
#include "t4c_decl.hpp"
#include "t4c_hrr_docs.hpp"
#include "t4c_hrr_decl.hpp"
#include "t4c_geom_hrr_body.hpp"

void
T4CGeomHrrCPUGenerator::generate(const std::string&        label,
                                 const int                 max_ang_mom,
                                 const std::array<int, 4>& geom_drvs) const
{
    if (_is_available(label))
    {
        if (geom_drvs == std::array<int, 4>({0, 1, 0, 0}))
        {
            for (int i = 0; i <= 3 * max_ang_mom; i++)
            {
                for (int j = 0; j <= 3 * max_ang_mom; j++)
                {
                    const auto integral = _get_integral(label, {i, j, 0, 0}, geom_drvs);
                    
                    _write_bra_hrr_cpp_header(integral);
                    
                    _write_bra_hrr_cpp_file(integral);
                }
            }
        }
        else if (geom_drvs == std::array<int, 4>({1, 1, 0, 0}))
        {
            for (int i = 0; i <= 2 * max_ang_mom; i++)
            {
                for (int j = 0; j <= 2 * max_ang_mom; j++)
                {
                    const auto integral = _get_integral(label, {i, j, 0, 0}, geom_drvs);
                    
                    _write_bra_hrr_cpp_header(integral);
                    
                    _write_bra_hrr_cpp_file(integral);
                }
            }
        }
        else
        {
            for (int i = 1; i <= 2 * max_ang_mom; i++)
            {
                for (int j = 0; j <= 2 * max_ang_mom; j++)
                {
                    const auto integral = _get_integral(label, {i, j, 0, 0}, geom_drvs);
                    
                    _write_bra_hrr_cpp_header(integral);
                    
                    _write_bra_hrr_cpp_file(integral);
                }
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
T4CGeomHrrCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "electron repulsion") return true;
    
    return false;
}

I4CIntegral
T4CGeomHrrCPUGenerator::_get_integral(const std::string&        label,
                                      const std::array<int, 4>& ang_moms,
                                      const std::array<int, 4>& geom_drvs) const
{
    // bra and ket sides

    const auto bpair = I2CPair("GA", ang_moms[0], "GB", ang_moms[1]);

    const auto kpair = I2CPair("GC", ang_moms[2], "GD", ang_moms[3]);

    VOperators prefixes;

    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[0])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[1])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[2])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[3])));

    // electron repulsion integrals

    if (fstr::lowercase(label) == "electron repulsion")
    {
        return I4CIntegral(bpair, kpair, Operator("1/|r-r'|"), 0, prefixes);
    }
    
    return I4CIntegral();
}

void
T4CGeomHrrCPUGenerator::_write_bra_hrr_cpp_header(const I4CIntegral& integral) const
{
    auto fname = t4c::bra_geom_hrr_file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_bra_hrr_hpp_defines(fstream, integral, true);
    
    _write_bra_hrr_hpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T4CHrrDocuDriver docs_drv;

    docs_drv.write_bra_geom_doc_str(fstream, integral);

    T4CHrrDeclDriver decl_drv;

    decl_drv.write_bra_geom_func_decl(fstream, integral, true);

    _write_namespace(fstream, integral, false);
   
    _write_bra_hrr_hpp_defines(fstream, integral, false);
    
    fstream.close();
}

void
T4CGeomHrrCPUGenerator::_write_bra_hrr_cpp_file(const I4CIntegral& integral) const
{
    auto fname = t4c::bra_geom_hrr_file_name(integral) + ".cpp";
        
    std::ofstream fstream;
        
    fstream.open(fname.c_str(), std::ios_base::trunc);
        
    _write_bra_hrr_cpp_includes(fstream, integral);

    _write_namespace(fstream, integral, true);

    T4CHrrDeclDriver decl_drv;
    
    decl_drv.write_bra_geom_func_decl(fstream, integral, false);

    T4CGeomHrrFuncBodyDriver func_drv;

    func_drv.write_bra_func_body(fstream, integral);
    
    fstream << std::endl;
    
    _write_namespace(fstream, integral, false);
        
    fstream.close();
}

void
T4CGeomHrrCPUGenerator::_write_bra_hrr_hpp_defines(      std::ofstream& fstream,
                                                   const I4CIntegral&   integral,
                                                   const bool           start) const
{
    auto fname = t4c::bra_geom_hrr_file_name(integral) + "_hpp";
    
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
T4CGeomHrrCPUGenerator::_write_bra_hrr_hpp_includes(      std::ofstream& fstream,
                                                    const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include <cstddef>"});
    
    lines.push_back({0, 0, 1, "#include \"Point.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"SimdArray.hpp\""});
        
    ost::write_code_lines(fstream, lines);
}

void
T4CGeomHrrCPUGenerator::_write_namespace(      std::ofstream& fstream,
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
T4CGeomHrrCPUGenerator::_write_bra_hrr_cpp_includes(      std::ofstream& fstream,
                                                    const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + t4c::bra_geom_hrr_file_name(integral) +  ".hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"TensorComponents.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}
