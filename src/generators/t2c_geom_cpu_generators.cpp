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

#include "t2c_geom_cpu_generators.hpp"

#include <iostream>

#include "v2i_center_driver.hpp"
#include "v2i_ovl_driver.hpp"
#include "v2i_kin_driver.hpp"
#include "v2i_dip_driver.hpp"
#include "v2i_npot_driver.hpp"
#include "v2i_linmom_driver.hpp"
#include "v2i_el_field_driver.hpp"
#include "v2i_eri_driver.hpp"
#include "v3i_ovl_driver.hpp"
#include "v3i_ovl_grad_driver.hpp"

#include "string_formater.hpp"
#include "file_stream.hpp"
#include "t2c_utils.hpp"
#include "t2c_docs.hpp"
#include "t2c_decl.hpp"
#include "t2c_body.hpp"

void
T2CGeomCPUGenerator::generate(const std::string&           label,
                              const int                    max_ang_mom,
                              const std::array<int, 3>&    geom_drvs,
                              const std::pair<bool, bool>& rec_form,
                              const bool                   use_rs) const
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
                       
                _write_cpp_header(geom_integrals, vrr_integrals, integral, geom_drvs, rec_form, use_rs);
                        
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
T2CGeomCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "overlap") return true;
    
    if (fstr::lowercase(label) == "kinetic energy") return true;
    
    if (fstr::lowercase(label) == "dipole momentum") return true;
    
    if (fstr::lowercase(label) == "nuclear potential") return true;
    
    if (fstr::lowercase(label) == "electron repulsion") return true;
    
    if (fstr::lowercase(label) == "three center overlap") return true;
    
    return false;
}

I2CIntegral
T2CGeomCPUGenerator::_get_integral(const std::string&        label,
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

    if (fstr::lowercase(label) == "dipole momentum")
    {
        return I2CIntegral(bra, ket, Operator("r", Tensor(1)), 0, prefixes);
    }

    // nuclear potential integrals and it's operator derivatives
    
    if (fstr::lowercase(label) == "nuclear potential")
    {
        if (geom_drvs[1] == 0)
        {
            return I2CIntegral(bra, ket, Operator("A"), 0, prefixes);
        }
        else
        {
            return I2CIntegral(bra, ket, Operator("AG", Tensor(geom_drvs[1])), 0, prefixes);
        }
    }
    
    // electron repulsion integrals
    
    if (fstr::lowercase(label) == "electron repulsion")
    {
        return I2CIntegral(bra, ket, Operator("1/|r-r'|"), 0, prefixes);
    }
    
    // three center overlap integrals
    
    if (fstr::lowercase(label) == "three center overlap")
    {
        if (geom_drvs[1] == 0)
        {
            return I2CIntegral(bra, ket, Operator("G(r)"), 0, prefixes);
        }
        
        if (geom_drvs[1] == 1)
        {
            return I2CIntegral(bra, ket, Operator("GX(r)", Tensor(1)), 0, prefixes);
        }
    }
    
    return I2CIntegral();
}

SI2CIntegrals
T2CGeomCPUGenerator::_generate_geom_integral_group(const I2CIntegral& integral) const
{
    V2ICenterDriver geom_drv;

    SI2CIntegrals tints;
    
    for (auto tint : geom_drv.apply_recursion({integral,}))
    {
        if (tint.prefixes().empty()) tints.insert(tint);
    }
    
    return tints;
}

SI2CIntegrals
T2CGeomCPUGenerator::_generate_vrr_integral_group(const I2CIntegral&   integral,
                                                  const SI2CIntegrals& integrals) const
{
    SI2CIntegrals tints(integrals);

    // Overlap integrals
    
    if (integral.integrand() == Operator("1"))
    {
        V2IOverlapDriver ovl_drv;
        
        tints = ovl_drv.create_recursion(tints);
    }

    // Dipole momentum integrals

    if (integral.integrand() == Operator("r", Tensor(1)))
    {
        V2IDipoleDriver dip_drv;

        tints = dip_drv.create_recursion(tints);

        V2IOverlapDriver ovl_drv;

        tints = ovl_drv.create_recursion(tints);
    }

    // Kinetic energy integrals
    
    if (integral.integrand() == Operator("T"))
    {
        V2IKineticEnergyDriver kin_drv;
        
        tints = kin_drv.create_recursion(tints);

        V2IOverlapDriver ovl_drv;

        tints = ovl_drv.create_recursion(tints);
    }
    
    // Nuclear potential integrals
    
    if (integral.integrand() == Operator("A"))
    {
        V2INuclearPotentialDriver npot_drv;
        
        tints = npot_drv.create_recursion(tints);
    }

    if (integral.integrand() == Operator("AG", integral.integrand().shape()))
    {
        V2IElectricFieldDriver el_field_drv;

        tints = el_field_drv.create_recursion(tints);
    
        V2INuclearPotentialDriver npot_drv;
        
        tints = npot_drv.create_recursion(tints);
    }
    
    if (integral.integrand() == Operator("1/|r-r'|"))
    {
        V2IElectronRepulsionDriver eri_drv;
        
        tints = eri_drv.create_recursion(tints);
    }
    
    if (integral.integrand() == Operator("G(r)"))
    {
        V3IOverlapDriver ovl_drv;
        
        tints = ovl_drv.create_recursion(tints);
    }
    
    if (integral.integrand() == Operator("G(r)"))
    {
        V3IOverlapDriver ovl_drv;
        
        tints = ovl_drv.create_recursion(tints);
    }
    
    if (integral.integrand() == Operator("GX(r)", integral.integrand().shape()))
    {
        V3IOverlapGradientDriver ovl_grad_drv;
        
        SI2CIntegrals cints;
        
        for (const auto& tint : tints)
        {
            auto rints = ovl_grad_drv.aux_vrr(tint);
            
            cints.insert(tint); 
            
            cints.insert(rints.begin(), rints.end());
        }
    
        V3IOverlapDriver ovl_drv;
        
        tints = ovl_drv.create_recursion(cints);
    }
    
    return tints;
}

void
T2CGeomCPUGenerator::_write_cpp_header(const SI2CIntegrals& geom_integrals,
                                       const SI2CIntegrals& vrr_integrals,
                                       const I2CIntegral&   integral,
                                       const std::array<int, 3>&    geom_drvs,
                                       const std::pair<bool, bool>& rec_form,
                                       const bool                   use_rs) const
{
    auto fname = _file_name(integral, rec_form, use_rs) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, rec_form, use_rs, false, true);

    _write_hpp_includes(fstream, vrr_integrals, integral, geom_drvs, rec_form, use_rs);

    _write_namespace(fstream, integral, true);

    T2CDocuDriver docs_drv;

    T2CDeclDriver decl_drv;

    T2CFuncBodyDriver func_drv;

    docs_drv.write_doc_str(fstream, integral, rec_form, use_rs);

    decl_drv.write_func_decl(fstream, integral, rec_form, use_rs, false);

    func_drv.write_func_body(fstream, geom_integrals, vrr_integrals, integral, geom_drvs, rec_form, use_rs);

    fstream << std::endl;

    _write_namespace(fstream, integral, false);

    _write_hpp_defines(fstream, integral, rec_form, use_rs, false, false);
    
    fstream.close();
}

std::string
T2CGeomCPUGenerator::_file_name(const I2CIntegral&           integral,
                                const std::pair<bool, bool>& rec_form,
                                const bool                   use_rs) const
{
    std::string label = (use_rs) ? "ErfRec" : "Rec";
    
    label += integral.label();
    
    if (rec_form.first) label = "Sum" + label;
    
    if (rec_form.second) label = "Conv" + label;
    
    return t2c::integral_label(integral) + label;
}


void
T2CGeomCPUGenerator::_write_hpp_defines(      std::ofstream&         fstream,
                                        const I2CIntegral&           integral,
                                        const std::pair<bool, bool>& rec_form,
                                        const bool                   use_rs,
                                        const bool                   is_prim_rec,
                                        const bool                   start) const
{
    auto fname = (is_prim_rec) ? t2c::prim_file_name(integral) : _file_name(integral, rec_form, use_rs) + "_hpp";
    
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
T2CGeomCPUGenerator::_write_hpp_includes(      std::ofstream&         fstream,
                                         const SI2CIntegrals&         integrals,
                                         const I2CIntegral&           integral,
                                         const std::array<int, 3>&    geom_drvs,
                                         const std::pair<bool, bool>& rec_form,
                                         const bool                   use_rs) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "#include <cstddef>"});
    
    lines.push_back({0, 0, 1, "#include <array>"});
    
    if (use_rs)
    {
        lines.push_back({0, 0, 1, "#include <vector>"});
    }
    
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
    
    lines.push_back({0, 0, 2, "#include \"" + t2c::geom_file_name(integral, geom_drvs) +  ".hpp\""});
    
    if ((integral.integrand().name() == "A") || (integral.integrand().name() == "AG"))
    {
        lines.push_back({0, 0, 1, "#include \"BoysFunc.hpp\""});
    }
    
    lines.push_back({0, 0, 1, "#include \"T2CUtils.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"T2CTransform.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"BatchFunc.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CGeomCPUGenerator::_write_namespace(      std::ofstream& fstream,
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
