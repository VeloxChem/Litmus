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


#include "t2c_geom_ecp_generators.hpp"

#include "string_formater.hpp"
#include "v2i_center_driver.hpp"
#include "v2i_translation_driver.hpp"
#include "v2i_loc_ecp_driver.hpp"
#include "t2c_utils.hpp"
#include "file_stream.hpp"
#include "t2c_docs.hpp"
#include "t2c_decl.hpp"
#include "t2c_ecp_body.hpp"

void
T2CECPGeomCPUGenerator::generate(const std::string&        label,
                                 const int                 max_ang_mom,
                                 const std::array<int, 3>& geom_drvs) const
{
    if (_is_available(label))
    {
        for (int i = 0; i <= max_ang_mom; i++)
        {
            for (int j = 0; j <= max_ang_mom; j++)
            {
                const auto integral = _get_integral(label, {i, j}, geom_drvs);
                        
                const auto geom_integrals = _generate_geom_integral_group(integral);
             
                const auto vrr_integrals = _generate_vrr_integral_group(integral, geom_integrals);
                       
                _write_cpp_header(geom_integrals, vrr_integrals, integral, geom_drvs);
                        
                std::cout << " *** REFERENCE: " << integral.prefix_label() << " | " << integral.label() << std::endl;
                        
                for (const auto& tint : geom_integrals)
                {
                    std::cout << " <>" << tint.prefix_label() << " | " << tint.label()  << " OP : " << tint.integrand().name() << std::endl;
                }
                       
                std::cout << " --- VRR --- " << std::endl;
                        
                for (const auto& tint : vrr_integrals)
                {
                    std::cout << " <>" << tint.prefix_label() << " | " << tint.label() << "_"  << tint.order() << " OP : " << tint.integrand().name() << std::endl;
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
T2CECPGeomCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "local") return true;
    
    return false;
}

I2CIntegral
T2CECPGeomCPUGenerator::_get_integral(const std::string&        label,
                                      const std::array<int, 2>& ang_moms,
                                      const std::array<int, 3>& geom_drvs) const
{
    // bra and ket sides
    
    const auto bra = I1CPair("GA", ang_moms[0]);
    
    const auto ket = I1CPair("GB", ang_moms[1]);
    
    // prefixes of integral bra, ket order
    
    VOperators prefixes;
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[0])));

    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[2])));
    
    // local ECP potential
    
    if (fstr::lowercase(label) == "local")
    {
        if (geom_drvs[1] > 0)
        {
            return I2CIntegral(bra, ket, Operator("U_L", Tensor(geom_drvs[1])), 0, prefixes);
        }
        else
        {
            return I2CIntegral(bra, ket, Operator("U_L"), 0, prefixes);
        }
    }
    
    return I2CIntegral();
}

SI2CIntegrals
T2CECPGeomCPUGenerator::_generate_geom_integral_group(const I2CIntegral& integral) const
{
    V2ICenterDriver geom_drv;
    
    SI2CIntegrals tints, rints;
    
    const auto gorder = integral.prefixes_sum_order();
    
    if (integral.integrand().shape().order() > 0)
    {
        V2ITranslationDriver trans_drv;
        
        rints = trans_drv.operator_vrr(integral);
        
        if ((integral.integrand().shape().order() == 1) && (gorder == 0))
        {
            return rints;
        }
        
        rints = geom_drv.apply_recursion(rints);
    }
    
    if (rints.empty()) rints = geom_drv.apply_recursion({integral,});
    
    for (auto rint : rints)
    {
        if (rint.prefixes().empty()) tints.insert(rint);
    }
    
    return tints;
}

SI2CIntegrals
T2CECPGeomCPUGenerator::_generate_vrr_integral_group(const I2CIntegral&   integral,
                                                     const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints(integrals);

    // Local ECP integrals
    
    if (integral.integrand().name() == "U_L")
    {
        V2ILocalECPDriver ecp_drv;
        
        tints = ecp_drv.create_full_recursion(tints);
    }
    
    return tints;
}

void
T2CECPGeomCPUGenerator::_write_cpp_header(const SI2CIntegrals&      geom_integrals,
                                          const SI2CIntegrals&      vrr_integrals,
                                          const I2CIntegral&        integral,
                                          const std::array<int, 3>& geom_drvs) const
{
    auto fname = _file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, true);

    _write_hpp_includes(fstream, vrr_integrals, integral, geom_drvs);

    _write_namespace(fstream, integral, true);

    T2CDocuDriver docs_drv;

    T2CDeclDriver decl_drv;

    T2CECPFuncBodyDriver func_drv;

    docs_drv.write_ecp_doc_str(fstream, integral);
    
    decl_drv.write_ecp_func_decl(fstream, integral, false);

    func_drv.write_func_body(fstream, geom_integrals, vrr_integrals, integral, geom_drvs);

    fstream << std::endl;

    _write_namespace(fstream, integral, false);

    _write_hpp_defines(fstream, integral, false);
    
    fstream.close();
}

std::string
T2CECPGeomCPUGenerator::_file_name(const I2CIntegral& integral) const
{
    std::string label = integral.label();
        
    return t2c::integral_label(integral) + label;
}

void
T2CECPGeomCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                           const I2CIntegral&   integral,
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
T2CECPGeomCPUGenerator::_write_hpp_includes(      std::ofstream&      fstream,
                                            const SI2CIntegrals&      integrals,
                                            const I2CIntegral&        integral,
                                            const std::array<int, 3>& geom_drvs) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "#include <cstddef>"});
    
    lines.push_back({0, 0, 1, "#include <array>"});
        
    lines.push_back({0, 0, 2, "#include <utility>"});
    
    lines.push_back({0, 0, 1, "#include \"GtoBlock.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"SimdArray.hpp\""});
    
    SI2CIntegrals rints;
    
    for (const auto& tint : integrals)
    {
        auto rint = tint;
        
        rint.set_order(0);
        
        rints.insert(rint);
    }
    
    for (const auto& rint : rints)
    {
        lines.push_back({0, 0, 1, "#include \"" + t2c::prim_file_name(rint) + ".hpp\""});
    }
    
    lines.push_back({0, 0, 1, "#include \"" + t2c::geom_file_name(integral, geom_drvs) +  ".hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T2CUtils.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T2CTransform.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"BatchFunc.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CECPGeomCPUGenerator::_write_namespace(      std::ofstream& fstream,
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
