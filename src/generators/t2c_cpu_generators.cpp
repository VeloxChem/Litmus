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

#include <iostream>

#include "t2c_defs.hpp"
#include "t2c_utils.hpp"
#include "t2c_cpu_generators.hpp"
#include "string_formater.hpp"
#include "file_stream.hpp"

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
                
                _write_cpp_header(integral, rec_form);
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
T2CCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "overlap") return true;
        
    return false;
}

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
    
    // nuclear potential integrals
    
    if (fstr::lowercase(label) == "nuclear potential")
    {
        return I2CIntegral(bra, ket, Operator("A"), 0, prefixes);
    }
    
    return I2CIntegral();
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
    
    _write_hpp_defines(fstream, integral, rec_form, true);
    
    _write_hpp_includes(fstream, integral, rec_form);
    
    _write_namespace(fstream, integral, true);
    
    T2CDocuDriver docs_drv;

    if ((integral[0] == integral[1]) && integral.is_simple())
    {
        docs_drv.write_doc_str(fstream, integral, rec_form,  true);
    }

    docs_drv.write_doc_str(fstream, integral, rec_form, false);

    _write_namespace(fstream, integral, false);
        
    _write_hpp_defines(fstream, integral, rec_form, false);
    
    fstream.close();
}

void
T2CCPUGenerator::_write_hpp_defines(      std::ofstream&         fstream,
                                    const I2CIntegral&           integral,
                                    const std::pair<bool, bool>& rec_form,
                                    const bool                   start) const
{
    auto fname = _file_name(integral, rec_form) + "_hpp";
    
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