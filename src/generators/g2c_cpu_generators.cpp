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

#include "g2c_cpu_generators.hpp"

#include "string_formater.hpp"
#include "file_stream.hpp"

#include "t2c_defs.hpp"
#include "t2c_utils.hpp"
#include "g2c_docs.hpp"
#include "g2c_decl.hpp"
#include "g2c_body.hpp"

#include "v2i_npot_driver.hpp"

void
G2CCPUGenerator::generate(const std::string&           label,
                          const int                    max_ang_mom,
                          const std::array<int, 3>&    geom_drvs,
                          const bool                   use_rs) const
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
                            const auto integral = _get_integral(label, {i, j}, geom_drvs);

                            const auto integrals = _generate_integral_group(integral, geom_drvs);

                            _write_cpp_header(integrals, integral, use_rs);
                            
//                            if (((i + j) > 0) && (!use_rs))
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
        std::cerr << "*** ERROR *** Unsupported type of two-center integral: ";
        
        std::cerr << label << " !!!" << std::endl;
        
        std::exit(EXIT_FAILURE);
    }
}

// MR: Need changes here
bool
G2CCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "nuclear potential") return true;

    return false;
}

I2CIntegral
G2CCPUGenerator::_get_integral(const std::string&        label,
                               const std::array<int, 2>& ang_moms,
                               const std::array<int, 3>& geom_drvs) const
{
    // bra and ket sides
    
    const auto bra = I1CPair("GA", ang_moms[0]);
    
    const auto ket = I1CPair("GB", ang_moms[1]);
    
    // nuclear potential integrals and it's operator derivatives
    
    if (fstr::lowercase(label) == "nuclear potential")
    {
        if (geom_drvs[1] == 0)
        {
            return I2CIntegral(bra, ket, Operator("A"), 0, {});
        }
        else
        {
            return I2CIntegral(bra, ket, Operator("AG", Tensor(geom_drvs[1])), 0, {});
        }
    }
    
    return I2CIntegral();
}

SI2CIntegrals
G2CCPUGenerator::_generate_integral_group(const I2CIntegral&        integral,
                                          const std::array<int, 3>& geom_drvs) const
{
    SI2CIntegrals tints;

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
            tints = npot_drv.create_recursion(tints);
        }
    }
    
    SI2CIntegrals rints;
    
    for (const auto& tint : tints)
    {
        if (tint.integrand().name() != "1")
        {
            rints.insert(tint);
        }
    }

    return rints;
}

std::string
G2CCPUGenerator::_file_name(const I2CIntegral& integral,
                            const bool         use_rs) const
{
    std::string label = (use_rs) ? "GridErfRec" : "GridRec";
    
    label += integral.label();
    
    return t2c::integral_label(integral) + label;
}

void
G2CCPUGenerator::_write_cpp_header(const SI2CIntegrals&         integrals,
                                   const I2CIntegral&           integral,
                                   const bool                   use_rs) const
{
    auto fname = _file_name(integral, use_rs) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral,  use_rs, false, true);
    
    _write_hpp_includes(fstream, integrals, integral, use_rs);
    
    _write_namespace(fstream, integral, true);
    
    G2CDocuDriver docs_drv;

    G2CDeclDriver decl_drv;

    G2CFuncBodyDriver func_drv;

    docs_drv.write_doc_str(fstream, integral, use_rs);
    
    decl_drv.write_func_decl(fstream, integral, use_rs, false);
    
    auto geom_drvs = std::array<int, 3>{0, 0, 0};
    
    func_drv.write_func_body(fstream, {}, integrals, integral, geom_drvs, use_rs);
    
    fstream << std::endl;

    _write_namespace(fstream, integral, false);
        
    _write_hpp_defines(fstream, integral, use_rs, false, false);
    
    fstream.close();
}

void
G2CCPUGenerator::_write_hpp_defines(      std::ofstream&         fstream,
                                    const I2CIntegral&           integral,
                                    const bool                   use_rs,
                                    const bool                   is_prim_rec,
                                    const bool                   start) const
{
    auto fname = (is_prim_rec) ? t2c::grid_prim_file_name(integral) : _file_name(integral, use_rs) + "_hpp";
    
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
G2CCPUGenerator::_write_hpp_includes(      std::ofstream&         fstream,
                                     const SI2CIntegrals&         integrals,
                                     const I2CIntegral&           integral,
                                     const bool                   use_rs) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "#include <cstddef>"});
    
    lines.push_back({0, 0, 1, "#include <array>"});
    
    lines.push_back({0, 0, 1, "#include <utility>"});
    
    lines.push_back({0, 0, 1, "#include <cmath>"});
        
    lines.push_back({0, 0, 1, "#include \"GtoBlock.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T2CUtils.hpp\""});
    
    if ((integral.integrand().name() == "A")  ||
        (integral.integrand().name() == "AG") ||
        (integral.integrand().name() == "1/|r-r'|"))
    {
        lines.push_back({0, 0, 1, "#include \"BoysFunc.hpp\""});
    }
    
    SI2CIntegrals rints;
    
    for (const auto& tint : integrals)
    {
        auto rint = tint;
        
        rint.set_order(0);
        
        rints.insert(rint);
    }
    
    for (const auto& rint : rints)
    {
        lines.push_back({0, 0, 1, "#include \"" + t2c::grid_prim_file_name(rint) + ".hpp\""});
    }
    
    lines.push_back({0, 0, 2, "#include \"MathConst.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
G2CCPUGenerator::_write_namespace(      std::ofstream& fstream,
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
