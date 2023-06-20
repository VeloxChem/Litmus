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
#include "t2c_decl.hpp"
#include "t2c_body.hpp"
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



bool
T2CCPUGenerator::_is_available(const std::string& label) const
{
    if (fstr::lowercase(label) == "overlap") return true;
    
    return false;
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
    return t2c::integral_label(integral) + "Rec" + integral.label();
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
    
    T2CDeclDriver decl_drv;
    
    if (integral[0] == integral[1])
    {
        docs_drv.write_doc_str(fstream, integral, true);
        
        decl_drv.write_func_decl(fstream, integral, true, true);
    }
    
    docs_drv.write_doc_str(fstream, integral, true);
    
    decl_drv.write_func_decl(fstream, integral, false, true);
    
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
    
    T2CDeclDriver decl_drv;
    
    T2CFuncBodyDriver func_drv;
    
    if (integral[0] == integral[1])
    {
        decl_drv.write_func_decl(fstream, integral, true, false);
        
        func_drv.write_func_body(fstream, integral, true);
    }
    
    decl_drv.write_func_decl(fstream, integral, false, false);
    
    func_drv.write_func_body(fstream, integral, true);
    
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
T2CCPUGenerator::_write_prim_funcs_to_cpp_header(      std::ofstream& fstream,
                                                 const I2CIntegral&   integral) const
{
    T2CDocuDriver docs_drv;
    
    T2CDeclDriver decl_drv;
    
    if (integral.is_simple_integrand() && integral.is_simple())
    {
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            docs_drv.write_prim_doc_str(fstream, integral);
            
            decl_drv.write_prim_func_decl(fstream, integral, true);
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                for (const auto& bcomp: Tensor(integral[0]).components())
                {
                    docs_drv.write_prim_doc_str(fstream, bcomp, integral, true);
                    
                    decl_drv.write_prim_func_decl(fstream, bcomp, integral, true, true);
                }
            }
            else
            {
                for (const auto& kcomp: Tensor(integral[1]).components())
                {
                    docs_drv.write_prim_doc_str(fstream, kcomp, integral, false);
                    
                    decl_drv.write_prim_func_decl(fstream, kcomp, integral, false, true);
                }
            }
        }
    }
    else
    {
        for (const auto& bcomp: Tensor(integral[0]).components())
        {
            for (const auto& kcomp: Tensor(integral[1]).components())
            {
                docs_drv.write_prim_doc_str(fstream, bcomp, kcomp, integral);
                
                decl_drv.write_prim_func_decl(fstream, bcomp, kcomp, integral, true);
            }
        }
    }
}

void
T2CCPUGenerator::_write_prim_funcs_to_cpp_file(      std::ofstream& fstream,
                                                 const I2CIntegral&   integral) const
{
    T2CDeclDriver decl_drv;

    if (integral.is_simple_integrand() && integral.is_simple())
    {
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            decl_drv.write_prim_func_decl(fstream, integral, false);
            
            _write_prim_func_body(fstream, integral); 
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                for (const auto& bcomp: Tensor(integral[0]).components())
                {
                    decl_drv.write_prim_func_decl(fstream, bcomp, integral, true, false);
                    
                    _write_prim_func_body(fstream, bcomp, integral, true);
                }
            }
            else
            {
                for (const auto& kcomp: Tensor(integral[1]).components())
                {
                    decl_drv.write_prim_func_decl(fstream, kcomp, integral, false, false);
                    
                    _write_prim_func_body(fstream, kcomp, integral, false);
                }
            }
        }
    }
    else
    {
        for (const auto& bcomp: Tensor(integral[0]).components())
        {
            for (const auto& kcomp: Tensor(integral[1]).components())
            {
                decl_drv.write_prim_func_decl(fstream, bcomp, kcomp, integral, false);
                
                _write_prim_func_body(fstream, bcomp, kcomp, integral);
            }
        }
    }
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
