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

#include "t2c_cpu_generator.hpp"

#include <iostream>
#include <iterator>

#include "operator.hpp"
#include "string_formater.hpp"
#include "file_stream.hpp"
#include "spherical_momentum.hpp"

#include "t2c_ovl_driver.hpp"
#include "t2c_docs.hpp"
#include "t2c_utils.hpp"

void
T2CCPUGenerator::generate(const std::string& label,
                          const int          angmom) const
{
    if (_is_available(label))
    {
        for (int i = 0; i <= angmom; i++)
        {
            for (int j = 0; j <= angmom; j++)
            {
                const auto integral = _get_integral(label, i, j);
                
                _write_cpp_header(integral);
                
                _write_cpp_file(integral);
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

int
T2CCPUGenerator::_get_tensor_component_index(const TensorComponent& component) const
{
    const auto tcomps = Tensor(component).components();
    
    for (int i = 0; i < tcomps.size(); i++)
    {
        if (tcomps[i] == component) return  i;
    }
    
    return -1;
}

std::string
T2CCPUGenerator::_combine_factors(const std::string& bra_factor,
                                  const std::string& ket_factor) const
{
    auto bra_label = bra_factor;
    
    auto ket_label = ket_factor;
    
    // get signs
    
    int bra_sign = (bra_factor[0] == '-') ? -1 : 1;
    
    int ket_sign = (ket_factor[0] == '-') ? -1 : 1;
    
    //  combinne symbolic factor
    
    if (bra_sign < 0) bra_label.erase(0, 1);
    
    if (ket_sign < 0) ket_label.erase(0, 1);
    
    std::string label;
    
    if (bra_label != "1.0") label = bra_label;
    
    if (ket_label != "1.0") label = (label.empty()) ? ket_label : label + " * " + ket_label;
    
    if (label.empty()) label = "1.0";
    
    if (bra_sign * ket_sign < 0) label.insert(0, "-");
    
    return label;
}

bool
T2CCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "overlap") return true;
    
    return false;
}

std::string
T2CCPUGenerator::_get_label(const I2CIntegral& integral) const
{
    if (integral.integrand() == Operator("1"))
    {
        return "Overlap";
    }
    
    return std::string();
}

std::string
T2CCPUGenerator::_get_namespace_label(const Operator& integrand) const
{
    auto labels = std::map<Operator, std::string>({ // list of operators
                                                   {Operator("1"), "ovlrec"},
                                                   });
    
    return labels[integrand];
}

std::string
T2CCPUGenerator::_get_matrix_symmetry(const Operator& integrand) const
{
    auto labels = std::map<Operator, std::string>({ // list of operators
                                                   {Operator("1"), "mat_t::symm"},
                                                   });
    
    return labels[integrand];
}

VT2CIntegrals
T2CCPUGenerator::_select_integral_components(const TensorComponent& component,
                                             const I2CIntegral&     integral,
                                             const bool             bra_first) const
{
    VT2CIntegrals tcomps;
    
    for (const auto& tcomp : integral.components<T1CPair, T1CPair>())
    {
        if (bra_first)
        {
            if (tcomp.bra().shape() == component) tcomps.push_back(tcomp);
        }
        else
        {
            if (tcomp.ket().shape() == component) tcomps.push_back(tcomp);
        }
    }
        
    return tcomps;
}

VT2CIntegrals
T2CCPUGenerator::_select_integral_components(const TensorComponent& bra_component,
                                             const TensorComponent& ket_component,
                                             const I2CIntegral&     integral) const
{
    VT2CIntegrals tcomps;
    
    for (const auto& tcomp : integral.components<T1CPair, T1CPair>())
    {
        if ((tcomp.bra().shape() == bra_component) &&
            (tcomp.ket().shape() == ket_component))
        {
            tcomps.push_back(tcomp);
        }
    }
        
    return tcomps;
}

I2CIntegral
T2CCPUGenerator::_get_integral(const std::string& label,
                               const int          ang_a,
                               const int          ang_b) const
{
    const auto bra = I1CPair("GA", ang_a);
    
    const auto ket = I1CPair("GB", ang_b);
    
    // overlap integrals
    
    if (fstr::lowercase(label) == "overlap")
    {
        return I2CIntegral(bra, ket, Operator("1"));
    }
    
    return I2CIntegral();
}

std::string
T2CCPUGenerator::_get_factor_label(const R2CTerm& rterm,
                                   const bool     first) const
{
    const auto pre_fact = rterm.prefactor();
    
    auto plabel = pre_fact.label();
    
    if (plabel == "1.0")  plabel = "";
    
    if (plabel == "-1.0") plabel = "-";
    
    if (pre_fact.denominator() != 1)
    {
        if (pre_fact.numerator() < 0) plabel.erase(0, 1);
        
        plabel = "(" + plabel + ")";
        
        if (pre_fact.numerator() < 0) plabel = "-" + plabel;
    }
    
    const auto facts = rterm.factors();
    
    std::string flabel;
    
    for (const auto& fact : facts)
    {
        const auto norder = rterm.factor_order(fact);
        
        for (size_t n = 0; n < norder; n++)
        {
            flabel += " * " + fact.label();
        }
    }
    
    // remove multiplication for special cases
    
    if ((pre_fact == Fraction(1)) || (pre_fact == Fraction(-1)))
    {
        flabel.erase(0, 3);
    }
    
    // merge labels
    
    flabel = plabel + flabel;
    
    if (!first)
    {
        if (flabel[0] == '-')
        {
            flabel.insert(1, " ");
        }
        else
        {
            flabel = "+ " + flabel;
        }
        
        flabel = " " + flabel;
    }
    
    return flabel;
}

bool
T2CCPUGenerator::_find_factor(const R2Group&     rgroup,
                              const std::string& label) const
{
    for (const auto& fact : rgroup.factors())
    {
        if (fact.label() == label) return true;
    }
    
    return false;
}

std::string
T2CCPUGenerator::_file_name(const I2CIntegral& integral) const
{
    return _get_label(integral) + "Rec" + integral.label();
}

void
T2CCPUGenerator::_write_cpp_header(const I2CIntegral& integral) const
{
    auto fname = _file_name(integral) + ".hpp";
    
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_hpp_defines(fstream, integral, true);
    
    _write_hpp_includes(fstream, integral);
    
    _write_namespace(fstream, integral, true);
    
    T2CDocuDriver docs_drv;
    
    if (integral[0] == integral[1])
    {
        docs_drv.write_doc_str(fstream, integral, true);
        
        _write_func_decl(fstream, integral, true, true);
    }
    
    docs_drv.write_doc_str(fstream, integral, true);
    
    _write_func_decl(fstream, integral, false, true);
    
    _write_prim_funcs_to_cpp_header(fstream, integral); 

    _write_namespace(fstream, integral, false);
    
    _write_hpp_defines(fstream, integral, false);

    fstream.close();
}

void
T2CCPUGenerator::_write_cpp_file(const I2CIntegral& integral) const
{
    auto fname = _file_name(integral) + ".cpp";
    
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    _write_cpp_includes(fstream, integral);
    
    _write_namespace(fstream, integral, true);
    
    if (integral[0] == integral[1])
    {
        _write_func_decl(fstream, integral, true, false);
        
        fstream << "{" << std::endl;
        
        _write_angmom_decl(fstream, integral);
        
        _write_gtos_decl(fstream, true);
        
        _write_ket_data_decl(fstream);
        
        _write_buffers_decl(fstream, integral);
        
        _write_batches_loop_start_decl(fstream, true);
        
        _write_main_call_tree_decl(fstream, integral, true);
        
        _write_batches_loop_end_decl(fstream);
        
        fstream << "}" << std::endl << std::endl;
    }
    
    _write_func_decl(fstream, integral, false, false);
    
    fstream << "{" << std::endl;
    
    _write_angmom_decl(fstream, integral);
    
    _write_gtos_decl(fstream, false);
    
    _write_ket_data_decl(fstream);
    
    _write_buffers_decl(fstream, integral);
    
    _write_batches_loop_start_decl(fstream, false);
    
    _write_main_call_tree_decl(fstream, integral, false);
    
    _write_batches_loop_end_decl(fstream);
    
    fstream << "}" << std::endl << std::endl;
    
    _write_prim_funcs_to_cpp_file(fstream, integral);

    _write_namespace(fstream, integral, false);

    fstream.close();
}

void
T2CCPUGenerator::_write_hpp_defines(      std::ofstream& fstream,
                                    const I2CIntegral&   integral,
                                    const bool           start) const
{
    const auto fname = _file_name(integral) + "_hpp";
    
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
T2CCPUGenerator::_write_hpp_includes(      std::ofstream& fstream,
                                     const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include <cstdint>"});
    
    lines.push_back({0, 0, 1, "#include \"GtoBlock.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"SubMatrix.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"SimdTypes.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"MatrixType.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_cpp_includes(      std::ofstream& fstream,
                                     const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 2, "#include \"" + _file_name(integral) +  ".hpp\""});
    
    lines.push_back({0, 0, 2, "#include <cmath>"});
    
    lines.push_back({0, 0, 1, "#include \"BatchFunc.hpp\""});
    
    lines.push_back({0, 0, 1, "#include \"MathConst.hpp\""});
    
    lines.push_back({0, 0, 2, "#include \"T2CDistributor.hpp\""});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_namespace(      std::ofstream& fstream,
                                  const I2CIntegral&   integral,
                                  const bool           start) const
{
    const auto label = _get_namespace_label(integral.integrand());
    
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
T2CCPUGenerator::_write_func_decl(      std::ofstream& fstream,
                                  const I2CIntegral&   integral,
                                  const bool           diagonal,
                                  const bool           terminus) const
{
    auto fname = "comp" + _get_label(integral) + integral.label();
    
    const auto fsize = fname.size() + 1;
    
    const auto padding = std::string(6, ' ');
    
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    if (const auto labels = t2c::integrand_components(integral.integrand(), "matrix"); labels.size() == 1)
    {
        lines.push_back({0, 0, 1, fname + "(" + padding + "CSubMatrix* matrix,"});
    }
    else
    {
        for (size_t i = 0; i < labels.size(); i++)
        {
            if (i == 0)
            {
                lines.push_back({0, 0, 1, fname + "(" + padding + "CSubMatrix* " + labels[i]});
            }
            else
            {
                lines.push_back({0, fsize, 1, padding + "CSubMatrix* " + labels[i]});
            }
        }
    }
    
    if (diagonal)
    {
        lines.push_back({0, fsize, 1, "const CGtoBlock&  gto_block,"});
    }
    else
    {
        lines.push_back({0, fsize, 1, "const CGtoBlock&  bra_gto_block,"});
        
        lines.push_back({0, fsize, 1, "const CGtoBlock&  ket_gto_block,"});
    }
    
    if (integral[0] != integral[1])
    {
        lines.push_back({0, fsize, 1, "const bool        ang_order,"});
    }
    
    lines.push_back({0, fsize, 1, "const int64_t     bra_first,"});
    
    fname = (terminus) ? ";" : "";
    
    if ((!diagonal) && (integral[0] == integral[1]))
    {
        lines.push_back({0, fsize, 1, "const int64_t     bra_last,"});
        
        lines.push_back({0, fsize, 2, "const mat_t       mat_type) -> void" + fname});
    }
    else
    {
        lines.push_back({0, fsize, 2, "const int64_t     bra_last) -> void" + fname});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_prim_funcs_to_cpp_header(      std::ofstream& fstream,
                                                 const I2CIntegral&   integral) const
{
    T2CDocuDriver docs_drv;
    
    if (integral.is_simple_integrand() && integral.is_simple())
    {
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            docs_drv.write_prim_doc_str(fstream, integral);
            
            _write_prim_func_decl(fstream, integral, true);
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                const auto bra = Tensor(integral[0]);
                
                for (const auto& bcomp: bra.components())
                {
                    _write_prim_func_docstr(fstream, bcomp, integral, true);
                    
                    _write_prim_func_decl(fstream, bcomp, integral, true, true);
                }
            }
            else
            {
                const auto ket = Tensor(integral[1]);
                
                for (const auto& kcomp: ket.components())
                {
                    _write_prim_func_docstr(fstream, kcomp, integral, false);
                    
                    _write_prim_func_decl(fstream, kcomp, integral, false, true);
                }
            }
        }
    }
    else
    {
        const auto bra = Tensor(integral[0]);
        
        const auto ket = Tensor(integral[1]);
        
        for (const auto& bcomp: bra.components())
        {
            for (const auto& kcomp: ket.components())
            {
                _write_prim_func_docstr(fstream, bcomp, kcomp, integral);
                
                _write_prim_func_decl(fstream, bcomp, kcomp, integral, true);
            }
        }
    }
}

void
T2CCPUGenerator::_write_prim_funcs_to_cpp_file(      std::ofstream& fstream,
                                                 const I2CIntegral&   integral) const
{
    if (const auto labels = t2c::integrand_components(integral.integrand(), "buffer"); labels.size() == 1)
    {
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            _write_prim_func_decl(fstream, integral, false);
            
            _write_prim_func_body(fstream, integral); 
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                const auto bra = Tensor(integral[0]);
                
                for (const auto& bcomp: bra.components())
                {
                    _write_prim_func_decl(fstream, bcomp, integral, true, false);
                    
                    _write_prim_func_body(fstream, bcomp, integral, true);
                }
            }
            else
            {
                const auto ket = Tensor(integral[1]);
                
                for (const auto& kcomp: ket.components())
                {
                    _write_prim_func_decl(fstream, kcomp, integral, false, false);
                    
                    _write_prim_func_body(fstream, kcomp, integral, false);
                }
            }
        }
    }
    else
    {
        const auto bra = Tensor(integral[0]);
        
        const auto ket = Tensor(integral[1]);
        
        for (const auto& bcomp: bra.components())
        {
            for (const auto& kcomp: ket.components())
            {
                _write_prim_func_decl(fstream, bcomp, kcomp, integral, false);
                
                _write_prim_func_body(fstream, bcomp, kcomp, integral);
            }
        }
    }
}



void
T2CCPUGenerator::_write_prim_func_docstr(      std::ofstream&   fstream,
                                         const TensorComponent& bra_component,
                                         const TensorComponent& ket_component,
                                         const I2CIntegral&     integral) const
{
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    const auto integrand = integral.integrand();
    
    // generate function name
    
    std::string fname = "<" + bra.label();
    
    fname += "_" + fstr::upcase(bra_component.label());
    
    fname += "|" + t2c::integrand_label(integrand) + "|";
    
    fname += ket.label() + "_" + fstr::upcase(ket_component.label());
    
    fname += ">";
    
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "/**"});
    
    lines.push_back({0, 1, 2, "Evaluates block of primitive " + fname + " integrals."});
    
    for (const auto& label : t2c::integrand_components(integrand, "buffer"))
    {
        lines.push_back({0, 1, 1, "@param " + label + " the partial integrals buffer."});
    }
    
    ost::write_code_lines(fstream, lines);
            
    _write_prim_data_docstr(fstream);
}

void
T2CCPUGenerator::_write_prim_func_decl(      std::ofstream& fstream,
                                       const I2CIntegral&   integral,
                                       const bool           terminus) const
{
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    const auto fname = "compPrimitive" + _get_label(integral) + integral.label();
    
    const auto fsize = fname.size() + 1;
    
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    std::vector<std::string> labels({"buffer", });
    
    if (integral[0] > 0) labels = t2c::tensor_components(bra, "buffer");
    
    if (integral[1] > 0) labels = t2c::tensor_components(ket, "buffer");
    
    lines.push_back({0, 0, 1, fname + "(      TDoubleArray& " + labels[0] + ","});
   
    for (size_t i = 1; i < labels.size(); i++)
    {
        lines.push_back({0, fsize + 6, 1, "TDoubleArray& " + labels[i] + ","});
    }
    
    ost::write_code_lines(fstream, lines);
    
    _write_prim_data_decl(fstream, fsize, terminus);
}

void
T2CCPUGenerator::_write_prim_func_decl(      std::ofstream&   fstream,
                                       const TensorComponent& component,
                                       const I2CIntegral&     integral,
                                       const bool             bra_first,
                                       const bool             terminus) const
{
    auto fname = "compPrimitive" + _get_label(integral) + integral.label();
    
    if (bra_first)
    {
        fname += "_" + fstr::upcase(component.label()) + "_T";
    }
    else
    {
        fname += "_T_" + fstr::upcase(component.label());
    }
    
    const auto fsize = fname.size() + 1;
    
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    const auto labels = (bra_first) ? t2c::tensor_components(integral[1], "buffer")
                                    : t2c::tensor_components(integral[0], "buffer");
    
    lines.push_back({0, 0, 1, fname + "(      TDoubleArray& " + labels[0] + ","});
   
    for (size_t i = 1; i < labels.size(); i++)
    {
        lines.push_back({0, fsize + 6, 1, "TDoubleArray& " + labels[i] + ","});
    }
    
    ost::write_code_lines(fstream, lines);
    
    _write_prim_data_decl(fstream, fsize, terminus);
}

void
T2CCPUGenerator::_write_prim_func_decl(      std::ofstream&   fstream,
                                       const TensorComponent& bra_component,
                                       const TensorComponent& ket_component,
                                       const I2CIntegral&     integral,
                                       const bool             terminus) const
{
    auto fname = "compPrimitive" + _get_label(integral) + integral.label();
    
    fname += "_" + fstr::upcase(bra_component.label());
   
    fname += "_" + fstr::upcase(ket_component.label());
    
    const auto fsize = fname.size() + 1;
    
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    const auto labels = t2c::integrand_components(integral.integrand(), "buffer");
    
    lines.push_back({0, 0, 1, fname + "(      TDoubleArray& " + labels[0] + ","});
   
    for (size_t i = 1; i < labels.size(); i++)
    {
        lines.push_back({0, fsize + 6, 1, "TDoubleArray& " + labels[i] + ","});
    }
    
    ost::write_code_lines(fstream, lines);
    
    _write_prim_data_decl(fstream, fsize, terminus);
}

void
T2CCPUGenerator::_write_prim_data_decl(      std::ofstream& fstream,
                                       const size_t         spacer,
                                       const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, spacer, 1, "const double        bra_exp,"});
    
    lines.push_back({0, spacer, 1, "const double        bra_norm,"});
    
    lines.push_back({0, spacer, 1, "const TPoint3D&     bra_coord,"});
    
    lines.push_back({0, spacer, 1, "const TDoubleArray& ket_exps,"});
    
    lines.push_back({0, spacer, 1, "const TDoubleArray& ket_norms,"});
    
    lines.push_back({0, spacer, 1, "const TDoubleArray& ket_coords_x,"});
    
    lines.push_back({0, spacer, 1, "const TDoubleArray& ket_coords_y,"});
    
    lines.push_back({0, spacer, 1, "const TDoubleArray& ket_coords_z,"});
    
    if (terminus)
    {
        lines.push_back({0, spacer, 2, "const int64_t       ket_dim) -> void;"});
    }
    else
    {
        lines.push_back({0, spacer, 1, "const int64_t       ket_dim) -> void"});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_gtos_decl(      std::ofstream& fstream,
                                  const bool           diagonal) const
{
    auto lines = VCodeLines();
    
    if (diagonal)
    {
        lines.push_back({1, 0, 2, "// intialize GTOs data"});
        
        lines.push_back({1, 0, 2, "const auto gto_coords = gto_block.getCoordinates();"});
       
        lines.push_back({1, 0, 2, "const auto gto_exps = gto_block.getExponents();"});
       
        lines.push_back({1, 0, 2, "const auto gto_norms = gto_block.getNormalizationFactors();"});
        
        lines.push_back({1, 0, 2, "const auto gto_indexes = gto_block.getOrbitalIndexes();"});
       
        lines.push_back({1, 0, 2, "const auto ncgtos = gto_block.getNumberOfBasisFunctions();"});

        lines.push_back({1, 0, 2, "const auto npgtos = gto_block.getNumberOfPrimitives();"});
    }
    else
    {
        lines.push_back({1, 0, 2, "// intialize GTOs data on bra side"});
        
        lines.push_back({1, 0, 2, "const auto bra_gto_coords = bra_gto_block.getCoordinates();"});
       
        lines.push_back({1, 0, 2, "const auto bra_gto_exps = bra_gto_block.getExponents();"});
       
        lines.push_back({1, 0, 2, "const auto bra_gto_norms = bra_gto_block.getNormalizationFactors();"});
        
        lines.push_back({1, 0, 2, "const auto bra_gto_indexes = bra_gto_block.getOrbitalIndexes();"});
       
        lines.push_back({1, 0, 2, "const auto bra_ncgtos = bra_gto_block.getNumberOfBasisFunctions();"});

        lines.push_back({1, 0, 2, "const auto bra_npgtos = bra_gto_block.getNumberOfPrimitives();"});
        
        lines.push_back({1, 0, 2, "// intialize GTOs data on ket side"});
        
        lines.push_back({1, 0, 2, "const auto ket_gto_coords = ket_gto_block.getCoordinates();"});
       
        lines.push_back({1, 0, 2, "const auto ket_gto_exps = ket_gto_block.getExponents();"});
       
        lines.push_back({1, 0, 2, "const auto ket_gto_norms = ket_gto_block.getNormalizationFactors();"});
        
        lines.push_back({1, 0, 2, "const auto ket_gto_indexes = ket_gto_block.getOrbitalIndexes();"});
       
        lines.push_back({1, 0, 2, "const auto ket_ncgtos = ket_gto_block.getNumberOfBasisFunctions();"});

        lines.push_back({1, 0, 2, "const auto ket_npgtos = ket_gto_block.getNumberOfPrimitives();"});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_angmom_decl(      std::ofstream& fstream,
                                    const I2CIntegral&   integral) const
{
    if ((integral[0] > 1) || (integral[1] > 1))
    {
        const auto angmom = SphericalMomentum(0);
        
        auto lines = VCodeLines();
        
        lines.push_back({1, 0, 2, "// spherical transformation factors"});
        
        if (integral[0] > 1)
        {
            for (const auto& fact : angmom.get_factors(integral[0]))
            {
                lines.push_back({1, 0, 2, "const double " + fact + ";"});
            }
        }
        
        if ((integral[1] > 1) && (integral[0] != integral[1]))
        {
            for (const auto& fact : angmom.get_factors(integral[1]))
            {
                lines.push_back({1, 0, 2, "const double " + fact + ";"});
            }
        }
        
        ost::write_code_lines(fstream, lines);
    }
}

void
T2CCPUGenerator::_write_ket_data_decl(std::ofstream& fstream) const
{
    auto lines = VCodeLines();
    
    lines.push_back({1, 0, 2, "// initialize aligned arrays for ket side"});
    
    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray ket_coords_x;"});
    
    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray ket_coords_y;"});
    
    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray ket_coords_z;"});
    
    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray ket_exps;"});
    
    lines.push_back({1, 0, 2, "alignas(64) TDoubleArray ket_norms;"});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_buffers_decl(      std::ofstream&   fstream,
                                     const I2CIntegral&     integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({1, 0, 2, "// initialize contracted integrals buffer"});
    
    auto labels = t2c::integrand_components(integral.integrand(), "buffer");

    if (labels.size() == 1)
    {
        const auto bra = Tensor(integral[0]);
        
        const auto ket = Tensor(integral[1]);
        
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            if (integral[0] > 0) labels = t2c::tensor_components(bra, "buffer");
            
            if (integral[1] > 0) labels = t2c::tensor_components(ket, "buffer");
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                labels = t2c::tensor_components(ket, "buffer");
            }
            else
            {
                labels = t2c::tensor_components(bra, "buffer");
            }
        }
    }
    
    for (const auto& label : labels)
    {
        lines.push_back({1, 0, 2, "alignas(64) TDoubleArray " + label + ";"});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_batches_loop_start_decl(      std::ofstream& fstream,
                                                const bool           diagonal) const
{
    auto lines = VCodeLines();
    
    lines.push_back({1, 0, 2, "// loop over integral batches"});
    
    if (diagonal)
    {
        lines.push_back({1, 0, 2, "const auto nbatches = batch::getNumberOfBatches(ncgtos, simd_width);"});
    }
    else
    {
        lines.push_back({1, 0, 2, "const auto nbatches = batch::getNumberOfBatches(ket_ncgtos, simd_width);"});
    }
        
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < nbatches; i++)"});
        
    lines.push_back({1, 0, 1, "{"});
    
    if (diagonal)
    {
        lines.push_back({2, 0, 2, "const auto [ket_first, ket_last] = batch::getBatchRange(i, ncgtos, simd_width);"});
    }
    else
    {
        lines.push_back({2, 0, 2, "const auto [ket_first, ket_last] = batch::getBatchRange(i, ket_ncgtos, simd_width);"});
    }
            
    lines.push_back({2, 0, 2, "const auto ket_dim = ket_last - ket_first;"});
            
    lines.push_back({2, 0, 1, "simd::loadCoordinates(ket_coords_x,"});
        
    lines.push_back({2, 22, 1, "ket_coords_y,"});
        
    lines.push_back({2, 22, 1, "ket_coords_z,"});
    
    if (diagonal)
    {
        lines.push_back({2, 22, 1, "gto_coords,"});
    }
    else
    {
        lines.push_back({2, 22, 1, "ket_gto_coords,"});
    }
        
    lines.push_back({2, 22, 1, "ket_first,"});
        
    lines.push_back({2, 22, 2, "ket_last);"});
        
    lines.push_back({2, 0, 1, "for (int64_t j = bra_first; j < bra_last; j++) "});
        
    lines.push_back({2, 0, 1, "{"});
        
    if (diagonal)
    {
        lines.push_back({3, 0, 2, "const auto bra_coord = gto_coords[j];"});
    }
    else
    {
        lines.push_back({3, 0, 2, "const auto bra_coord = bra_gto_coords[j];"});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_batches_loop_end_decl(std::ofstream& fstream) const
{
    auto lines = VCodeLines();
    
    lines.push_back({2, 0, 1, "}"});
    
    lines.push_back({1, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_main_call_tree_decl(      std::ofstream& fstream,
                                            const I2CIntegral&   integral,
                                            const bool           diagonal) const
{
    if (const auto labels = t2c::integrand_components(integral.integrand(), "buffer"); labels.size() == 1)
    {
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            _write_prim_call_tree_block_decl(fstream, integral, diagonal);
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                const auto bra = Tensor(integral[0]);
                
                for (const auto& bcomp: bra.components())
                {
                    _write_prim_call_tree_block_decl(fstream, bcomp, integral, true, diagonal);
                    
                    fstream << std::endl;
                }
            }
            else
            {
                const auto ket = Tensor(integral[1]);
                
                for (const auto& kcomp: ket.components())
                {
                    _write_prim_call_tree_block_decl(fstream, kcomp, integral, false, diagonal);
                    
                    fstream << std::endl;
                }
            }
        }
    }
    else
    {
        const auto bra = Tensor(integral[0]);
        
        const auto ket = Tensor(integral[1]);
        
        for (const auto& bcomp: bra.components())
        {
            for (const auto& kcomp: ket.components())
            {
                _write_prim_call_tree_block_decl(fstream, bcomp, kcomp, integral, diagonal);
                
                fstream << std::endl; 
            }
        }
    }
}

void
T2CCPUGenerator::_write_prim_call_tree_block_decl(      std::ofstream& fstream,
                                                  const I2CIntegral&   integral,
                                                  const bool           diagonal) const
{
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    std::vector<std::string> labels({"buffer", });
    
    if (integral[0] > 0) labels = t2c::tensor_components(bra, "buffer");
    
    if (integral[1] > 0) labels = t2c::tensor_components(ket, "buffer");
    
    auto lines = VCodeLines();
    
    lines.push_back({3, 0, 2, "// compute primitive integrals block"});
    
    for (const auto& label : labels)
    {
        lines.push_back({3, 0, 2, "simd::zero(" + label + ");"});
    }
    
    ost::write_code_lines(fstream, lines);
    
    lines.clear();
    
    _write_primitives_loop_start_decl(fstream, diagonal);
    
    auto fname = _get_namespace_label(integral.integrand());
    
    fname += "::compPrimitive" + _get_label(integral) + integral.label();
    
    const auto fsize = fname.size() + 1;
    
    lines.push_back({5, 0, 1, fname + "(" + labels[0] + ","});
   
    for (size_t i = 1; i < labels.size(); i++)
    {
        lines.push_back({5, fsize, 1, labels[i] + ","});
    }
    
    ost::write_code_lines(fstream, lines);
    
    _write_primitives_call_data_decl(fstream, fsize);
    
    // write primitives loop end
    
    _write_primitives_loop_end_decl(fstream);
    
    _write_block_distributor_decl(fstream, integral, diagonal);
}

void
T2CCPUGenerator::_write_prim_call_tree_block_decl(      std::ofstream&   fstream,
                                                  const TensorComponent& component,
                                                  const I2CIntegral&     integral,
                                                  const bool             bra_first,
                                                  const bool             diagonal) const
{
    const auto labels = (bra_first) ? t2c::tensor_components(integral[1], "buffer")
                                    : t2c::tensor_components(integral[0], "buffer");
    
    auto lines = VCodeLines();
    
    lines.push_back({3, 0, 2, "// compute primitive integrals block (" + fstr::upcase(component.label()) + ")"});
    
    for (const auto& label : labels)
    {
        lines.push_back({3, 0, 2, "simd::zero(" + label + ");"});
    }
    
    ost::write_code_lines(fstream, lines);
    
    lines.clear();
    
    _write_primitives_loop_start_decl(fstream, diagonal);
    
    auto fname = _get_namespace_label(integral.integrand());
    
    fname += "::compPrimitive" + _get_label(integral) + integral.label();
    
    if (bra_first)
    {
        fname += "_" + fstr::upcase(component.label()) + "_T";
    }
    else
    {
        fname += "_T_" + fstr::upcase(component.label());
    }
    
    const auto fsize = fname.size() + 1;
    
    lines.push_back({5, 0, 1, fname + "(" + labels[0] + ","});
   
    for (size_t i = 1; i < labels.size(); i++)
    {
        lines.push_back({5, fsize, 1, labels[i] + ","});
    }
    
    ost::write_code_lines(fstream, lines);
    
    _write_primitives_call_data_decl(fstream, fsize);
    
    // write primitives loop end
    
    _write_primitives_loop_end_decl(fstream);
    
    _write_block_distributor_decl(fstream, component, integral, bra_first, diagonal);
}

void
T2CCPUGenerator::_write_prim_call_tree_block_decl(      std::ofstream&   fstream,
                                                  const TensorComponent& bra_component,
                                                  const TensorComponent& ket_component,
                                                  const I2CIntegral&     integral,
                                                  const bool             diagonal) const
{
   
    
    const auto labels = t2c::integrand_components(integral.integrand(), "buffer");
    
    auto lines = VCodeLines();
    
    lines.push_back({3, 0, 2, "// compute primitive integrals block (" +
                   
                    fstr::upcase(bra_component.label()) + "_" +
                   
                    fstr::upcase(ket_component.label()) + ")"});
    
    for (const auto& label : labels)
    {
        lines.push_back({3, 0, 2, "simd::zero(" + label + ");"});
    }
    
    ost::write_code_lines(fstream, lines);
    
    lines.clear();
    
    _write_primitives_loop_start_decl(fstream, diagonal);
    
    auto fname = _get_namespace_label(integral.integrand());
    
    fname += "compPrimitive" + _get_label(integral) + integral.label();
    
    fname += "_" + fstr::upcase(bra_component.label());
   
    fname += "_" + fstr::upcase(ket_component.label());
    
    const auto fsize = fname.size() + 1;
    
    lines.push_back({5, 0, 1, fname + "(" + labels[0] + ","});
   
    for (size_t i = 1; i < labels.size(); i++)
    {
        lines.push_back({5, fsize, 1, labels[i] + ","});
    }
    
    ost::write_code_lines(fstream, lines);
    
    _write_primitives_call_data_decl(fstream, fsize);
    
    // write primitives loop end
    
    _write_primitives_loop_end_decl(fstream);
    
    _write_block_distributor_decl(fstream, bra_component, ket_component, integral, diagonal);
}

void
T2CCPUGenerator::_write_primitives_loop_start_decl(      std::ofstream& fstream,
                                                   const bool           diagonal) const
{
    auto lines = VCodeLines();
    
    if (diagonal)
    {
        lines.push_back({3, 0, 1, "for (int64_t k = 0; k < npgtos; k++)"});
        
        lines.push_back({3, 0, 1, "{"});
        
        lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);"});
            
        lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);"});
            
        lines.push_back({4, 0, 1, "for (int64_t l = 0; l < npgtos; l++)"});
        
        lines.push_back({4, 0, 1, "{"});
        
        lines.push_back({5, 0, 2, "const auto bra_index = l * ncgtos + j;"});
                
        lines.push_back({5, 0, 2, "const auto bra_exp = gto_exps[bra_index];"});
                
        lines.push_back({5, 0, 2, "const auto bra_norm = gto_norms[bra_index];"});
    }
    else
    {
        lines.push_back({3, 0, 1, "for (int64_t k = 0; k < ket_npgtos; k++)"});
        
        lines.push_back({3, 0, 1, "{"});
        
        lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);"});
            
        lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);"});
            
        lines.push_back({4, 0, 1, "for (int64_t l = 0; l < bra_npgtos; l++)"});
        
        lines.push_back({4, 0, 1, "{"});
        
        lines.push_back({5, 0, 2, "const auto bra_index = l * bra_ncgtos + j;"});
                
        lines.push_back({5, 0, 2, "const auto bra_exp = bra_gto_exps[bra_index];"});
                
        lines.push_back({5, 0, 2, "const auto bra_norm = bra_gto_norms[bra_index];"});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_primitives_loop_end_decl(std::ofstream& fstream) const
{
    auto lines = VCodeLines();
    
    lines.push_back({4, 0, 1, "}"});
    
    lines.push_back({3, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_primitives_call_data_decl(      std::ofstream& fstream,
                                                  const size_t         spacer) const
{
    auto lines = VCodeLines();
   
    lines.push_back({5, spacer, 1, "bra_exp,"});
    
    lines.push_back({5, spacer, 1, "bra_norm,"});
    
    lines.push_back({5, spacer, 1, "bra_coord,"});
    
    lines.push_back({5, spacer, 1, "ket_exps,"});
    
    lines.push_back({5, spacer, 1, "ket_norms,"});
    
    lines.push_back({5, spacer, 1, "ket_coords_x,"});
    
    lines.push_back({5, spacer, 1, "ket_coords_y,"});
    
    lines.push_back({5, spacer, 1, "ket_coords_z,"});
    
    lines.push_back({5, spacer, 1, "ket_dim);"});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_block_distributor_decl(      std::ofstream& fstream,
                                               const I2CIntegral&   integral,
                                               const bool           diagonal) const
{
    auto lines = VCodeLines();
    
    if ((integral[0] + integral[1]) == 0)
    {
        if (diagonal)
        {
            lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, buffer, gto_indexes,"});
                            
            lines.push_back({3, 20, 1, "0, 0, j, ket_first, ket_last);"});
        }
        else
        {
            lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, buffer, bra_gto_indexes, ket_gto_indexes,"});
            
            lines.push_back({3, 20, 1, "0, 0, j, ket_first, ket_last, mat_type);"});
        }
    }
    
    if (integral[0] > 0)
    {
        const auto bra_mom = SphericalMomentum(integral[0]);
        
        const auto bra = Tensor(integral[0]);
        
        const auto labels = t2c::tensor_components(bra, "buffer");
        
        for (int i = 0; i < labels.size(); i++)
        {
            for (const auto& pair : bra_mom.select_pairs(i))
            {
                auto flabel = labels[i];
                
                flabel += (pair.second == "1.0")  ? "" : ", "  + pair.second;
                
                const auto findex = std::to_string(pair.first);
                  
                lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, " + flabel +
                                          ", bra_gto_indexes, ket_gto_indexes,"});
                    
                lines.push_back({3, 20, 2, findex + ", 0, j, ket_first, ket_last, ang_order);"});
            }
        }
    }
    
    if (integral[1] > 0)
    {
        const auto ket_mom = SphericalMomentum(integral[1]);
        
        const auto ket = Tensor(integral[1]);
        
        const auto labels = t2c::tensor_components(ket, "buffer");
        
        for (int i = 0; i < labels.size(); i++)
        {
            for (const auto& pair : ket_mom.select_pairs(i))
            {
                auto flabel = labels[i];
                
                flabel += (pair.second == "1.0")  ? "" : ", "  + pair.second;
                
                const auto findex = std::to_string(pair.first);
              
                lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, " + flabel +
                                          ", bra_gto_indexes, ket_gto_indexes,"});
                    
                lines.push_back({3, 20, 2, "0, " + findex + ", j, ket_first, ket_last, ang_order);"});
            }
        }
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_block_distributor_decl(      std::ofstream&   fstream,
                                               const TensorComponent& component,
                                               const I2CIntegral&     integral,
                                               const bool             bra_first,
                                               const bool             diagonal) const
{
    const auto labels = (bra_first) ? t2c::tensor_components(integral[1], "buffer")
                                    : t2c::tensor_components(integral[0], "buffer");
    
    const auto bra_mom = (bra_first) ? SphericalMomentum(integral[0])
                                     : SphericalMomentum(integral[1]);
    
    const auto ket_mom = (bra_first) ? SphericalMomentum(integral[1])
                                     : SphericalMomentum(integral[0]);
    
    const auto index = _get_tensor_component_index(component);
    
    const auto bra_pairs = bra_mom.select_pairs(index);
    
    auto lines = VCodeLines();
    
    for (int i = 0; i < labels.size(); i++)
    {
        for (const auto& ket_pair : ket_mom.select_pairs(i))
        {
            for (const auto& bra_pair : bra_pairs)
            {
                const auto lfactor = _combine_factors(bra_pair.second,
                                                      ket_pair.second);
                
                auto flabel = labels[i];
                
                flabel += (lfactor == "1.0")  ? "" : ", "  + lfactor;
                
                auto ijlabel = (bra_first) ? std::to_string(bra_pair.first) + ", " + std::to_string(ket_pair.first)
                                           : std::to_string(ket_pair.first) + ", " + std::to_string(bra_pair.first);
                
                if  (integral[0] ==  integral[1])
                {
                    if (diagonal)
                    {
                        lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, " + flabel + ", gto_indexes,"});
                                        
                        lines.push_back({3, 20, 2, ijlabel + ", j, ket_first, ket_last);"});
                    }
                    else
                    {
                        lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, " + flabel + ", bra_gto_indexes, ket_gto_indexes,"});
                        
                        lines.push_back({3, 20, 2, ijlabel + ", j, ket_first, ket_last, mat_type);"});
                    }
                }
                else
                {
                    lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, " + flabel +
                                               ", bra_gto_indexes, ket_gto_indexes,"});
                        
                    lines.push_back({3, 20, 2, ijlabel + ", j, ket_first, ket_last, ang_order);"});
                    
                }
            }
        }
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_block_distributor_decl(      std::ofstream&   fstream,
                                               const TensorComponent& bra_component,
                                               const TensorComponent& ket_component,
                                               const I2CIntegral&     integral,
                                               const bool             diagonal) const
{
    const auto labels = t2c::integrand_components(integral.integrand(), "buffer");
    
    const auto matrices = t2c::integrand_components(integral.integrand(), "matrix");
    
    const auto bra_mom = SphericalMomentum(integral[0]);
    
    const auto bra_index = _get_tensor_component_index(bra_component);
    
    const auto bra_pairs = bra_mom.select_pairs(bra_index);
    
    const auto ket_mom = SphericalMomentum(integral[1]);
    
    const auto ket_index = _get_tensor_component_index(ket_component);
    
    const auto ket_pairs = ket_mom.select_pairs(ket_index);
    
    auto lines = VCodeLines();
    
    for (int i = 0; i < labels.size(); i++)
    {
        for (const auto& bra_pair : bra_pairs)
        {
            for (const auto& ket_pair : ket_pairs)
            {
                const auto lfactor = _combine_factors(bra_pair.second,
                                                      ket_pair.second);
                
                auto flabel = labels[i];
                
                flabel += (lfactor == "1.0")  ? "" : ", "  + lfactor;
                
                auto ijlabel = std::to_string(bra_pair.first) + ", " + std::to_string(ket_pair.first);
                
                if  (integral[0] ==  integral[1])
                {
                    if (diagonal)
                    {
                        lines.push_back({3, 0, 1, "t2cfunc::distribute(" + matrices[i] + ", " + flabel +
                                                  ", gto_indexes,"});
                        
                        lines.push_back({3, 20, 2, ijlabel + ", j, ket_first, ket_last);"});
                    }
                    else
                    {
                        lines.push_back({3, 0, 1, "t2cfunc::distribute(" + matrices[i] + ", " + flabel +
                                                  ", bra_gto_indexes, ket_gto_indexes,"});
                        
                        lines.push_back({3, 20, 2, ijlabel + ", j, ket_first, ket_last, mat_type);"});
                    }
                }
                else
                {
                    lines.push_back({3, 0, 1, "t2cfunc::distribute(" + matrices[i] + ", " + flabel +
                                              ", bra_gto_indexes, ket_gto_indexes,"});
                    
                    lines.push_back({3, 20, 2, ijlabel + ", j, ket_first, ket_last, ang_order);"});
                }
            }
        }
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_prim_func_body(      std::ofstream& fstream,
                                       const I2CIntegral&   integral) const
{
    fstream << "{" << std::endl;
    
    _write_prim_func_common_data(fstream);
    
    _write_prim_func_buffers(fstream, integral);
    
    _write_prim_func_pragma(fstream, integral);
    
    _write_prim_func_loop_start(fstream, integral);
    
    const auto tcomps = integral.components<T1CPair, T1CPair>();
    
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    std::vector<std::string> labels({"fints", });
    
    if (integral[0] > 0) labels = t2c::tensor_components(bra, "fints");
    
    if (integral[1] > 0) labels = t2c::tensor_components(ket, "fints");
    
    _write_simd_code(fstream, labels, tcomps, integral);
    
    _write_prim_func_loop_end(fstream);
    
    fstream << "}" << std::endl << std::endl;
}

void
T2CCPUGenerator::_write_prim_func_body(      std::ofstream&   fstream,
                                       const TensorComponent& component,
                                       const I2CIntegral&     integral,
                                       const bool             bra_first) const
{
    fstream << "{" << std::endl;
    
    _write_prim_func_common_data(fstream);
    
    _write_prim_func_buffers(fstream, component, integral, bra_first);
    
    _write_prim_func_pragma(fstream, component, integral, bra_first);
    
    _write_prim_func_loop_start(fstream, integral);
    
    const auto tcomps = _select_integral_components(component, integral, bra_first);
    
    const auto labels = (bra_first) ? t2c::tensor_components(integral[1], "fints")
                                    : t2c::tensor_components(integral[0], "fints");
    
    _write_simd_code(fstream, labels, tcomps, integral);
    
    _write_prim_func_loop_end(fstream);
    
    fstream << "}" << std::endl << std::endl;
}

void
T2CCPUGenerator::_write_prim_func_body(      std::ofstream&   fstream,
                                       const TensorComponent& bra_component,
                                       const TensorComponent& ket_component,
                                       const I2CIntegral&     integral) const
{
    fstream << "{" << std::endl;
    
    _write_prim_func_common_data(fstream);
    
    _write_prim_func_buffers(fstream, bra_component, ket_component, integral);
    
    _write_prim_func_pragma(fstream, bra_component, ket_component, integral);
    
    _write_prim_func_loop_start(fstream, integral);
    
    const auto tcomps = _select_integral_components(bra_component, ket_component, integral);
    
    const auto labels = t2c::integrand_components(integral.integrand(), "fints");
    
    _write_simd_code(fstream, labels, tcomps, integral);
    
    _write_prim_func_loop_end(fstream);
    
    fstream << "}" << std::endl << std::endl;
}

void
T2CCPUGenerator::_write_prim_func_common_data(std::ofstream& fstream) const
{
    auto lines = VCodeLines();
    
    lines.push_back({1, 0, 2, "// set up math constants"});
    
    lines.push_back({1, 0, 2, "const auto fpi = mathconst::getPiValue();"});
    
    lines.push_back({1, 0, 2, "// set up coordinates for bra side"});
    
    lines.push_back({1, 0, 2, "const auto bra_rx = bra_coord[0];"});
    
    lines.push_back({1, 0, 2, "const auto bra_ry = bra_coord[1];"});
    
    lines.push_back({1, 0, 2, "const auto bra_rz = bra_coord[2];"});
    
    lines.push_back({1, 0, 2, "// set up coordinates for ket side"});
    
    lines.push_back({1, 0, 2, "auto ket_rx = ket_coords_x.data();"});
    
    lines.push_back({1, 0, 2, "auto ket_ry = ket_coords_y.data();"});
    
    lines.push_back({1, 0, 2, "auto ket_rz = ket_coords_z.data();"});
    
    lines.push_back({1, 0, 2, "// set exponents and normalization factors on ket side"});
    
    lines.push_back({1, 0, 2, "auto ket_fe = ket_exps.data();"});
    
    lines.push_back({1, 0, 2, "auto ket_fn = ket_norms.data();"});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_prim_func_buffers(      std::ofstream& fstream,
                                          const I2CIntegral&   integral) const
{
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    std::vector<std::string> labels({"buffer", });
    
    if (integral[0] > 0) labels = t2c::tensor_components(bra, "buffer");
    
    if (integral[1] > 0) labels = t2c::tensor_components(ket, "buffer");
    
    std::vector<std::string> flabels({"fints", });
    
    if (integral[0] > 0) flabels = t2c::tensor_components(bra, "fints");
    
    if (integral[1] > 0) flabels = t2c::tensor_components(ket, "fints");
    
    auto lines = VCodeLines();
    
    lines.push_back({1, 0, 2, "// set up pointer to integrals buffer(s)"});
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        lines.push_back({1, 0, 2, "auto " + flabels[i] + " = " + labels[i] + ".data();"});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_prim_func_buffers(      std::ofstream&   fstream,
                                          const TensorComponent& component,
                                          const I2CIntegral&     integral,
                                          const bool             bra_first) const
{
    const auto labels  = (bra_first) ? t2c::tensor_components(integral[1], "buffer")
                                     : t2c::tensor_components(integral[0], "buffer");
    
    const auto flabels = (bra_first) ? t2c::tensor_components(integral[1], "fints")
                                     : t2c::tensor_components(integral[0], "fints");
    
    auto lines = VCodeLines();
    
    lines.push_back({1, 0, 2, "// set up pointer to integrals buffer(s)"});
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        lines.push_back({1, 0, 2, "auto " + flabels[i] + " = " + labels[i] + ".data();"});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_prim_func_buffers(      std::ofstream&   fstream,
                                          const TensorComponent& bra_component,
                                          const TensorComponent& ket_component,
                                          const I2CIntegral&     integral) const
{
    const auto labels  = t2c::integrand_components(integral.integrand(), "buffer");
    
    const auto flabels = t2c::integrand_components(integral.integrand(), "fints");
    
    auto lines = VCodeLines();
    
    lines.push_back({1, 0, 2, "// set up pointer to integrals buffer(s)"});
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        lines.push_back({1, 0, 2, "auto " + flabels[i] + " = " + labels[i] + ".data();"});
    }
    
    ost::write_code_lines(fstream, lines);
}


void
T2CCPUGenerator::_write_prim_func_pragma(      std::ofstream& fstream,
                                         const I2CIntegral&   integral) const
{
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    std::vector<std::string> labels({"fints", });
    
    if (integral[0] > 0) labels = t2c::tensor_components(bra, "fints");
    
    if (integral[1] > 0) labels = t2c::tensor_components(ket, "fints");
    
    auto lines = VCodeLines();
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        if (i == 0)
        {
            lines.push_back({1, 0, 1, "#pragma omp simd aligned(" + labels[i] + ",\\"});
        }
        else
        {
            lines.push_back({1, 25, 1, labels[i] + ",\\"});
        }
    }
    
    ost::write_code_lines(fstream, lines);
    
    _write_prim_func_common_pragma(fstream);
}

void
T2CCPUGenerator::_write_prim_func_pragma(      std::ofstream&   fstream,
                                         const TensorComponent& component,
                                         const I2CIntegral&     integral,
                                         const bool             bra_first) const
{
    const auto labels = (bra_first) ? t2c::tensor_components(integral[1], "fints")
                                    : t2c::tensor_components(integral[0], "fints");
    
    auto lines = VCodeLines();
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        if (i == 0)
        {
            lines.push_back({1, 0, 1, "#pragma omp simd aligned(" + labels[i] + ",\\"});
        }
        else
        {
            lines.push_back({1, 25, 1, labels[i] + ",\\"});
        }
    }
    
    ost::write_code_lines(fstream, lines);
    
    _write_prim_func_common_pragma(fstream);
}

void
T2CCPUGenerator::_write_prim_func_pragma(      std::ofstream&   fstream,
                                         const TensorComponent& bra_component,
                                         const TensorComponent& ket_component,
                                         const I2CIntegral&     integral) const
{
    const auto labels = t2c::integrand_components(integral.integrand(), "fints");
    
    auto lines = VCodeLines();
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        if (i == 0)
        {
            lines.push_back({1, 0, 1, "#pragma omp simd aligned(" + labels[i] + ",\\"});
        }
        else
        {
            lines.push_back({1, 25, 1, labels[i] + ",\\"});
        }
    }
    
    ost::write_code_lines(fstream, lines);
 
    _write_prim_func_common_pragma(fstream);
}

void
T2CCPUGenerator::_write_prim_func_common_pragma(std::ofstream& fstream) const
{
    auto lines = VCodeLines();
    
    lines.push_back({1, 25, 1, "ket_fe,\\"});
    
    lines.push_back({1, 25, 1, "ket_fn,\\"});
    
    lines.push_back({1, 25, 1, "ket_rx,\\"});
    
    lines.push_back({1, 25, 1, "ket_ry,\\"});
    
    lines.push_back({1, 25, 1, "ket_rz : 64)"});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_prim_func_loop_start(      std::ofstream& fstream,
                                             const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();

    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < ket_dim; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "const auto ab_x = bra_rx - ket_rx[i];"});
    
    lines.push_back({2, 0, 2, "const auto ab_y = bra_ry - ket_ry[i];"});
    
    lines.push_back({2, 0, 2, "const auto ab_z = bra_rz - ket_rz[i];"});
    
    lines.push_back({2, 0, 2, "const auto fe_0 = 1.0 / (bra_exp + ket_fe[i]);"});
    
    lines.push_back({2, 0, 2, "auto fz_0 = bra_exp * ket_fe[i] * fe_0;"});
    
    lines.push_back({2, 0, 2, "fz_0 *= (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);"});
    
    if ((integral.integrand() == Operator("1")) && ((integral[0] + integral[1]) == 0))
    {
        lines.push_back({2, 0, 1, "fints[i] += bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);"});
    }
    else
    {
        lines.push_back({2, 0, 2, "const auto fss = bra_norm * ket_fn[i] * std::pow(fe_0 * fpi, 1.50) * std::exp(-fz_0);"});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_prim_func_loop_end(std::ofstream& fstream) const
{
    auto lines = VCodeLines();
    
    lines.push_back({1, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_simd_code(      std::ofstream&            fstream,
                                  const std::vector<std::string>& labels,
                                  const VT2CIntegrals&            components,
                                  const I2CIntegral&              integral) const
{
    R2Group rgroup;
    
    // Overlap inntegrals
    
    if (integral.integrand() == Operator("1"))
    {
        T2COverlapDriver t2c_ovl_drv;
        
        rgroup = t2c_ovl_drv.create_recursion(components);
    }
    
    // ... other integrals
    
    auto lines = VCodeLines();

    // prefactors
    
    if (_find_factor(rgroup, "rpa_x"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_x = -ket_fe[i] * ab_x * fe_0;"});
    }
    
    if (_find_factor(rgroup, "rpa_y"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_y = -ket_fe[i] * ab_y * fe_0;"});
    }
    
    if (_find_factor(rgroup, "rpa_z"))
    {
        lines.push_back({2, 0, 2, "const auto rpa_z = -ket_fe[i] * ab_z * fe_0;"});
    }
    
    if (_find_factor(rgroup, "rpb_x"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_x = bra_exp * ab_x * fe_0;"});
    }
    
    if (_find_factor(rgroup, "rpb_y"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_y = bra_exp * ab_y * fe_0;"});
    }
    
    if (_find_factor(rgroup, "rpb_z"))
    {
        lines.push_back({2, 0, 2, "const auto rpb_z = bra_exp * ab_z * fe_0;"});
    }
    
    // generate loop simd code

    for (size_t i = 0; i < labels.size(); i++)
    {
        const auto rdist = rgroup[i];
        
        const auto nterms = rdist.terms();
        
        auto nbatches = nterms / 5;
        
        if ((nterms % 5) != 0) nbatches++;
        
        for (size_t j = 0; j < nbatches; j++)
        {
            const auto sterm = 5 * j;
            
            const auto eterm = ((sterm + 5) > nterms) ? nterms : sterm + 5;
            
            std::string simd_str;
            
            for (size_t k = sterm; k < eterm; k++)
            {
                simd_str += _get_factor_label(rdist[k], k == sterm);
            }
            
            if ((eterm - sterm) > 1) simd_str = "(" + simd_str + ")";
            
            int shift = 2;
            
            if ((labels.size() == (i + 1)) && (nbatches == (j + 1))) shift = 1;
            
            lines.push_back({2, 0, shift, labels[i] + "[i] += fss * " + simd_str + ";"});
        }
    }
    
    ost::write_code_lines(fstream, lines);
}
