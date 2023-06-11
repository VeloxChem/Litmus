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

std::map<Operator, std::string>
T2CCPUGenerator::_get_integrands_map() const
{
    return std::map<Operator, std::string>({{Operator("1"), ""}, });
}

std::string
T2CCPUGenerator::_get_integrand_label(const Operator& integrand) const
{
    auto labels = std::map<Operator, std::string>({// list of operators
                                                   {Operator("1"), ""},
                                                  });
    
    return labels[integrand];
}

std::string
T2CCPUGenerator::_get_namespace_label(const Operator& integrand) const
{
    auto labels = std::map<Operator, std::string>({ // list of operators
                                                   {Operator("1"), "ovlrec"},
                                                   });
    
    return labels[integrand];
}

std::vector<std::string>
T2CCPUGenerator::_get_operator_components(const Operator&    integrand,
                                          const std::string& label) const
{
    if (const auto op_comps = integrand.components(); op_comps.size() == 1)
    {
        return std::vector<std::string>({label,});
    }
    else
    {
        std::vector<std::string> op_labels;
        
        for (const auto& op_comp : op_comps)
        {
            op_labels.push_back(label + "_" + op_comp.label());
        }
        
        return op_labels;
    }
}

std::vector<std::string>
T2CCPUGenerator::_get_tensor_components(const Tensor&      tensor,
                                        const std::string& label) const
{
    if (const auto tcomps = tensor.components(); tcomps.size() == 1)
    {
        return std::vector<std::string>({label,});
    }
    else
    {
        std::vector<std::string> tlabels;
        
        for (const auto& tcomp : tcomps)
        {
            tlabels.push_back(label + "_" + tcomp.label());
        }
        
        return tlabels;
    }
}

std::map<Operator, std::string>
T2CCPUGenerator::_get_namespaces_map() const
{
    return std::map<Operator, std::string>({{Operator("1"), "ovlrec"},});
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
    
    if (integral[0] == integral[1])
    {
        _write_func_docstr(fstream, integral, true);
        
        _write_func_decl(fstream, integral, true, true);
    }
    
    _write_func_docstr(fstream, integral, false);
    
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
    
    _write_gtos_decl(fstream, false);
    
    _write_ket_data_decl(fstream);
    
    _write_buffers_decl(fstream, integral);
    
    _write_batches_loop_start_decl(fstream, false);
    
    _write_main_call_tree_decl(fstream, integral, false);
    
    _write_batches_loop_end_decl(fstream);
    
    fstream << "}" << std::endl << std::endl;
    
    //_write_prim_funcs_to_cpp_header(fstream, integral);

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
    
    lines.push_back({0, 0, 2, "#include \"SimdTypes.hpp\""});
    
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
T2CCPUGenerator::_write_func_docstr(      std::ofstream& fstream,
                                    const I2CIntegral&   integral,
                                    const bool           diagonal) const
{
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    const auto integrand = integral.integrand();
    
    auto lines = VCodeLines();
    
    auto label = " Evaluates <" + bra.label() + "|";
    
    label += _get_integrand_label(integrand);
    
    label += "|" + ket.label() + ">  integrals for given ";
    
    label += (diagonal) ? "GTOs block." : "pair of GTOs blocks.";
    
    lines.push_back({0, 0, 1, "/**"});
    
    lines.push_back({0, 0, 2, label});
    
    if (const auto labels = _get_operator_components(integrand, "matrix"); labels.size() == 1)
    {
        lines.push_back({0, 1, 1, "@param matrix the pointer to matrix for storage of integrals."});
    }
    else
    {
        for (const auto& label : labels)
        {
            auto lcomp = fstr::upcase(label);
            
            lcomp.erase(0, lcomp.find('_') + 1);
            
            lcomp = "@param " + label + "the pointer to matrix for storage of Cartesian integral component " + lcomp + ".";
            
            lines.push_back({0, 1, 1, lcomp});
        }
    }
    
    if (diagonal)
    {
        lines.push_back({0, 1, 1, "@param gto_block the GTOs block."});
    }
    else
    {
        lines.push_back({0, 1, 1, "@param bra_gto_block the GTOs block on bra side."});
        
        lines.push_back({0, 1, 1, "@param ket_gto_block the GTOs block on ket side."});
    }
    
    lines.push_back({0, 1, 1, "@param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side."});
    
    lines.push_back({0, 1, 1, "@param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side."});
    
    lines.push_back({0, 0, 1, "*/"});
    
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
    
    if (const auto labels = _get_operator_components(integral.integrand(), "matrix"); labels.size() == 1)
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
    
    lines.push_back({0, fsize, 1, "const int64_t     bra_first,"});
    
    fname = (terminus) ? ";" : "";
    
    lines.push_back({0, fsize, 2, "const int64_t     bra_last) -> void" + fname});
        
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_prim_funcs_to_cpp_header(      std::ofstream& fstream,
                                                 const I2CIntegral&   integral) const
{
    if (const auto labels = _get_operator_components(integral.integrand(), "buffer"); labels.size() == 1)
    {
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            _write_prim_func_docstr(fstream, integral);
            
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
T2CCPUGenerator::_write_prim_func_docstr(      std::ofstream& fstream,
                                         const I2CIntegral&   integral) const
{
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    // generate function name
    
    auto fname = "<" + bra.label() + "|" ;
    
    fname += _get_integrand_label(integral.integrand());
    
    fname += "|" + ket.label() + ">";
    
    // write code
    
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "/**"});
    
    lines.push_back({0, 1, 2, "Evaluates block of primitive " + fname + " integrals."});
    
    std::vector<std::string> labels({"buffer", });
    
    if (integral[0] > 0) labels = _get_tensor_components(bra, "buffer");
    
    if (integral[1] > 0) labels = _get_tensor_components(ket, "buffer");
    
    if (labels.empty())
    {
        lines.push_back({0, 1, 1, "@param buffer the integrals buffer."});
    }
    else
    {
        for (const auto& label : labels)
        {
            lines.push_back({0, 1, 1, "@param " + label + " the partial integrals buffer."});
        }
    }
    
    ost::write_code_lines(fstream, lines);
        
    _write_prim_data_docstr(fstream);
}

void
T2CCPUGenerator::_write_prim_func_docstr(      std::ofstream&   fstream,
                                         const TensorComponent& component,
                                         const I2CIntegral&     integral,
                                         const bool             bra_first) const
{
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    // generate function name
    
    std::string fname = "<" + bra.label();
    
    fname += (bra_first) ? "_" + fstr::upcase(component.label()) : "";
    
    fname += "|" + _get_integrand_label(integral.integrand()) + "|";
    
    fname += ket.label();
    
    fname += (bra_first) ? "" : "_" + fstr::upcase(component.label());
    
    fname += ">";
    
    // write code
    
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "/**"});
    
    lines.push_back({0, 1, 2, "Evaluates block of primitive " + fname + " integrals."});
    
    const auto labels = (bra_first) ? _get_tensor_components(ket, "buffer")
                                    : _get_tensor_components(bra, "buffer");
    
    for (const auto& label : labels)
    {
        lines.push_back({0, 1, 1, "@param " + label + " the partial integrals buffer."});
    }
    
    ost::write_code_lines(fstream, lines);
    
    _write_prim_data_docstr(fstream);
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
    
    fname += "|" + _get_integrand_label(integrand) + "|";
    
    fname += ket.label() + "_" + fstr::upcase(ket_component.label());
    
    fname += ">";
    
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "/**"});
    
    lines.push_back({0, 1, 2, "Evaluates block of primitive " + fname + " integrals."});
    
    for (const auto& label : _get_operator_components(integrand, "buffer"))
    {
        lines.push_back({0, 1, 1, "@param " + label + " the partial integrals buffer."});
    }
    
    ost::write_code_lines(fstream, lines);
            
    _write_prim_data_docstr(fstream);
}

void
T2CCPUGenerator::_write_prim_data_docstr(std::ofstream& fstream) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 1, 1, "@param bra_exp the primitive exponent on bra side."});
    
    lines.push_back({0, 1, 1, "@param bra_norm the primitive normalization factor on bra side."});
    
    lines.push_back({0, 1, 1, "@param bra_coord the 3d coordinate of basis function on bra side."});
    
    lines.push_back({0, 1, 1, "@param ket_exps the array of primitive exponents on ket side."});
    
    lines.push_back({0, 1, 1, "@param ket_norms the array of primitive normalization factors on ket side."});
    
    lines.push_back({0, 1, 1, "@param ket_coords_x the array of Cartesian X coordinates on ket side."});
    
    lines.push_back({0, 1, 1, "@param ket_coords_y the array of Cartesian Y coordinates on ket side."});
    
    lines.push_back({0, 1, 1, "@param ket_coords_z the array of Cartesian Z coordinates on ket side."});
    
    lines.push_back({0, 1, 1, "@param ket_dim the end size of ket arrays." });
    
    lines.push_back({0, 1, 1, "@param ket_dim the end size of ket arrays." });
    
    lines.push_back({0, 0, 1, "*/"});
    
    ost::write_code_lines(fstream, lines);
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
    
    if (integral[0] > 0) labels = _get_tensor_components(bra, "buffer");
    
    if (integral[1] > 0) labels = _get_tensor_components(ket, "buffer");
    
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
    
    const auto labels = (bra_first) ? _get_tensor_components(integral[1], "buffer")
                                    : _get_tensor_components(integral[0], "buffer");
    
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
    
    const auto labels = _get_operator_components(integral.integrand(), "buffer");
    
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
    
    auto labels = _get_operator_components(integral.integrand(), "buffer");

    if (labels.size() == 1)
    {
        const auto bra = Tensor(integral[0]);
        
        const auto ket = Tensor(integral[1]);
        
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            if (integral[0] > 0) labels = _get_tensor_components(bra, "buffer");
            
            if (integral[1] > 0) labels = _get_tensor_components(ket, "buffer");
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                labels = _get_tensor_components(ket, "buffer");
            }
            else
            {
                labels = _get_tensor_components(bra, "buffer");
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
    
    lines.push_back({2, 0, 2, "}"});
    
    lines.push_back({1, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CCPUGenerator::_write_main_call_tree_decl(      std::ofstream& fstream,
                                            const I2CIntegral&   integral,
                                            const bool           diagonal) const
{
    if (const auto labels = _get_operator_components(integral.integrand(), "buffer"); labels.size() == 1)
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
    
    if (integral[0] > 0) labels = _get_tensor_components(bra, "buffer");
    
    if (integral[1] > 0) labels = _get_tensor_components(ket, "buffer");
    
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
}

void
T2CCPUGenerator::_write_prim_call_tree_block_decl(      std::ofstream&   fstream,
                                                  const TensorComponent& component,
                                                  const I2CIntegral&     integral,
                                                  const bool             bra_first,
                                                  const bool             diagonal) const
{
    const auto labels = (bra_first) ? _get_tensor_components(integral[1], "buffer")
                                    : _get_tensor_components(integral[0], "buffer");
    
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
}

void
T2CCPUGenerator::_write_prim_call_tree_block_decl(      std::ofstream&   fstream,
                                                  const TensorComponent& bra_component,
                                                  const TensorComponent& ket_component,
                                                  const I2CIntegral&     integral,
                                                  const bool             diagonal) const
{
   
    
    const auto labels = _get_operator_components(integral.integrand(), "buffer");
    
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
    
    lines.push_back({3, 0, 1, "}"});
    
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
