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

#include "operator.hpp"
#include "string_formater.hpp"

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
        
        fstream << "}" << std::endl << std::endl;
    }
    
    _write_func_decl(fstream, integral, false, false);
    
    fstream << "{" << std::endl;
    
    _write_gtos_decl(fstream, false);
    
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
    
    if (start)
    {
        fstream << "#ifndef " << fname << std::endl;
        
        fstream << "#define " << fname  << std::endl;
        
        fstream << std::endl;
    }
    else
    {
        fstream << "#endif /* " << fname << " */";
    }
}

void
T2CCPUGenerator::_write_hpp_includes(      std::ofstream& fstream,
                                     const I2CIntegral&   integral) const
{
    // C/C++ headers
    
    fstream << "#include <cstdint>" << std::endl;
    
    // custom headers
    
    fstream << std::endl;

    fstream << "#include \"GtoBlock.hpp\"" << std::endl;
    
    fstream << "#include \"SubMatrix.hpp\"" << std::endl;

    fstream << "#include \"SimdTypes.hpp\"" << std::endl;
    
    fstream << std::endl;
}

void
T2CCPUGenerator::_write_cpp_includes(      std::ofstream& fstream,
                                     const I2CIntegral&   integral) const
{
    // main header
    
    fstream << "#include \"" + _file_name(integral) +  ".hpp\"" << std::endl;
    
    // C/C++ headers
    
    fstream << std::endl;
    
    fstream << "#include <cmath>" << std::endl;
    
    // custom headers
    
    fstream << std::endl;
    
    fstream << "#include \"BatchFunc.hpp\"" << std::endl;
    
    fstream << "#include \"MathConst.hpp\"" << std::endl;

    fstream << "#include \"T2CDistributor.hpp\"" << std::endl;
    
    fstream << std::endl;
}

void
T2CCPUGenerator::_write_namespace(      std::ofstream& fstream,
                                  const I2CIntegral&   integral,
                                  const bool           start) const
{
    auto tlabels = std::map<Operator, std::string>({{Operator("1"), "ovlrec"},});
    
    const auto label = tlabels[integral.integrand()];
    
    if (start)
    {
        fstream << "namespace "  << label << " { // " << label << " namespace" << std::endl;
    }
    else
    {
        fstream << "} // " << label << " namespace" << std::endl;
    }
    
    fstream << std::endl;
}

void
T2CCPUGenerator::_write_func_docstr(      std::ofstream& fstream,
                                    const I2CIntegral&   integral,
                                    const bool           diagonal) const
{
    auto tlabels = _get_integrands_map();
    
    const auto label = tlabels[integral.integrand()];
    
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    fstream << "/**" << std::endl;
    
    fstream << " Evaluates <" << bra.label() << "|" << label << "|" << ket.label() << ">  integrals for given ";
    
    if (diagonal)
    {
        fstream << "GTOs block." << std::endl;
    }
    else
    {
        fstream << "pair of GTOs blocks." << std::endl;
    }
    
    fstream << std::endl;
    
    fstream << " @param matrix the pointer to matrix for storage of integrals."  << std::endl;
    
    if (diagonal)
    {
        fstream <<  " @param gto_block the GTOs block." << std::endl; 
    }
    else
    {
        fstream << " @param bra_gto_block the GTOs block on bra side." << std::endl;
        
        fstream << " @param ket_gto_block the GTOs block on ket side." << std::endl;
    }
    
    fstream << " @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side." << std::endl;
    
    fstream << " @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side." << std::endl;
    
    fstream << "*/" << std::endl;
}

void
T2CCPUGenerator::_write_matrix_docstr(      std::ofstream& fstream,
                                      const I2CIntegral&   integral) const
{
    const auto op = integral.integrand();
    
    if (const auto op_comps = op.components(); op_comps.size() == 1)
    {
        fstream << " @param matrix " << std::endl;
    }
    else
    {
        for (const auto& op_comp : op_comps)
        {
            fstream << " @param matrix_" << op_comp.label();
            
            fstream << " the pointer to matrix for storage of";
            
            fstream << " Cartesian integral component ";
            
            fstream << fstr::upcase(op_comp.label()) << " ." << std::endl;
        }
    }
}

void
T2CCPUGenerator::_write_func_decl(      std::ofstream& fstream,
                                  const I2CIntegral&   integral,
                                  const bool           diagonal,
                                  const bool           terminus) const
{
    const auto fname = "comp" + _get_label(integral) + integral.label();
    
    const auto spacer = std::string(fname.size() + 1,  ' ');
    
    fstream << "auto" << std::endl;
    
    fstream << fname << "(";
    
    _write_matrix_decl(fstream, integral, spacer);
    
    if (diagonal)
    {
        fstream << spacer << "const CGtoBlock&  gto_block," << std::endl;
    }
    else
    {
        fstream << spacer << "const CGtoBlock&  bra_gto_block," << std::endl;
        
        fstream << spacer << "const CGtoBlock&  ket_gto_block," << std::endl;
    }
    
    fstream << spacer << "const int64_t     bra_first," << std::endl;
    
    fstream << spacer << "const int64_t     bra_last) -> void";
    
    if (terminus) fstream << ";" << std::endl;
    
    fstream << std::endl;
}

void
T2CCPUGenerator::_write_matrix_decl(      std::ofstream& fstream,
                                    const I2CIntegral&   integral,
                                    const std::string&   spacer) const
{
    const auto op = integral.integrand();
    
    const auto padding = std::string(6, ' ');
    
    if (const auto op_comps = op.components(); op_comps.size() == 1)
    {
        fstream << padding << "CSubMatrix* matrix," << std::endl;
    }
    else
    {
        for (size_t i = 0; i < op_comps.size(); i++)
        {
            if (i > 0) fstream << spacer;
            
            fstream << padding << "CSubMatrix* matrix_";
            
            fstream << op_comps[i].label() << "," << std::endl;
        }
    }
}

void
T2CCPUGenerator::_write_prim_funcs_to_cpp_header(      std::ofstream& fstream,
                                                 const I2CIntegral&   integral) const
{
    const auto op = integral.integrand();
    
    if (const auto op_comps = op.components(); op_comps.size() == 1)
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
    auto tlabels = _get_integrands_map();
    
    const auto label = tlabels[integral.integrand()];
    
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    fstream << "/**" << std::endl;
        
    fstream <<  " Evaluates block of primitive <" << bra.label() << "|";
        
    fstream << label << "|" << ket.label() << ">  integrals.";
        
    fstream << std::endl << std::endl;
        
    if ((integral[0] + integral[1]) == 0)
    {
        fstream << " @param buffer the integrals buffer." << std::endl;
    }
        
    if (integral[0] > 0)
    {
        for (const auto& bcomp : bra.components())
        {
            fstream << " @param buffer_" << bcomp.label();
                
            fstream << " the partial integrals buffer." << std::endl;
        }
    }
        
    if (integral[1] > 0)
    {
        
        for (const auto& kcomp : ket.components())
        {
            fstream << " @param buffer_" << kcomp.label();
                
            fstream << " the partial integrals buffer." << std::endl;
        }
    }
        
    _write_prim_data_docstr(fstream);
}

void
T2CCPUGenerator::_write_prim_func_docstr(      std::ofstream&   fstream,
                                         const TensorComponent& component,
                                         const I2CIntegral&     integral,
                                         const bool             bra_first) const
{
    auto tlabels = _get_integrands_map();
    
    const auto label = tlabels[integral.integrand()];
    
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    fstream << "/**" << std::endl;
    
    if (bra_first)
    {
        fstream <<  " Evaluates block of primitive <";
        
        fstream << bra.label() << "_" << fstr::upcase(component.label());
        
        fstream << "|" << label << "|" << ket.label();
        
        fstream << ">  integrals." <<  std::endl << std::endl;
        
        for (const auto& kcomp : ket.components())
        {
            fstream << " @param buffer_" << kcomp.label();
        
            fstream << " the partial integrals buffer." << std::endl;
        }
    }
    else
    {
        fstream <<  " Evaluates block of primitive <";
        
        fstream  << bra.label() << "|" << label << "|";
        
        fstream << ket.label() << "_" << fstr::upcase(component.label());
        
        fstream << ">  integrals." << std::endl << std::endl;
        
        for (const auto& bcomp : bra.components())
        {
            fstream << " @param buffer_" << bcomp.label();
        
            fstream << " the partial integrals buffer." << std::endl;
        }
    }
    
    _write_prim_data_docstr(fstream);
}

void
T2CCPUGenerator::_write_prim_func_docstr(      std::ofstream&   fstream,
                                         const TensorComponent& bra_component,
                                         const TensorComponent& ket_component,
                                         const I2CIntegral&     integral) const
{
    const auto op = integral.integrand();
    
    auto tlabels = _get_integrands_map();
    
    const auto label = tlabels[op];
    
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    fstream << "/**" << std::endl;
            
    fstream <<  " Evaluates block of primitive <";
            
    fstream << bra.label() << "_" << fstr::upcase(bra_component.label());
            
    fstream << "|" << label << "|";
            
    fstream << ket.label() << "_" << fstr::upcase(ket_component.label());
            
    fstream << ">  integrals." << std::endl << std::endl;
            
    for (const auto& op_comp : op.components())
    {
        fstream << " @param buffer_" << op_comp.label();
                
        fstream << " the partial integrals buffer." << std::endl;
    }
            
    _write_prim_data_docstr(fstream);
}

void
T2CCPUGenerator::_write_prim_data_docstr(std::ofstream& fstream) const
{
    fstream << " @param bra_exp the primitive exponent on bra side." << std::endl;
    
    fstream << " @param bra_norm the primitive normalization factor on bra side." << std::endl;
    
    fstream << " @param bra_coord the 3d coordinate of basis function on bra side." << std::endl;
    
    fstream << " @param ket_exps the array of primitive exponents on ket side." << std::endl;
    
    fstream << " @param ket_norms the array of primitive normalization factors on ket side." << std::endl;
    
    fstream << " @param ket_coords_x the array of Cartesian X coordinates on ket side." << std::endl;
    
    fstream << " @param ket_coords_y the array of Cartesian Y coordinates on ket side." << std::endl;
    
    fstream << " @param ket_coords_z the array of Cartesian Z coordinates on ket side." << std::endl;
    
    fstream << " @param ket_dim the end size of ket arrays." << std::endl;
    
    fstream << " @param ket_dim the end size of ket arrays." << std::endl;
    
    fstream << "*/" << std::endl;
    
}

void
T2CCPUGenerator::_write_prim_func_decl(      std::ofstream& fstream,
                                       const I2CIntegral&   integral,
                                       const bool           terminus) const
{
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    const auto fname = "compPrimitive" + _get_label(integral) + integral.label();
    
    const auto spacer = std::string(fname.size() + 1,  ' ');
    
    const auto padding = std::string(6, ' ');
    
    fstream << "auto" << std::endl;
    
    fstream << fname << "(";
    
    if ((integral[0] + integral[1]) == 0)
    {
        fstream << padding << "TDoubleArray& buffer," << std::endl;
    }
        
    if (integral[0] > 0)
    {
        const auto bra_comps = bra.components();
        
        for (size_t i = 0; i < bra_comps.size(); i++)
        {
            if (i > 0) fstream << spacer;
            
            fstream << padding << "TDoubleArray& buffer_";
            
            fstream << bra_comps[i].label() << "," << std::endl;
        }
    }
        
    if (integral[1] > 0)
    {
        const auto ket_comps = ket.components();
        
        for (size_t i = 0; i < ket_comps.size(); i++)
        {
            if (i > 0) fstream << spacer;
            
            fstream << padding << "TDoubleArray& buffer_";
            
            fstream << ket_comps[i].label() << "," << std::endl;
        }
    }
    
    _write_prim_data_decl(fstream, spacer, terminus);
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
    
    const auto spacer = std::string(fname.size() + 1,  ' ');
    
    const auto padding = std::string(6, ' ');
    
    fstream << "auto" << std::endl;
    
    fstream << fname << "(";
    
    if (bra_first)
    {
        const auto ket = Tensor(integral[1]);
        
        const auto ket_comps = ket.components();
        
        for (size_t i = 0; i < ket_comps.size(); i++)
        {
            if (i > 0) fstream << spacer;
            
            fstream << padding << "TDoubleArray& buffer_";
            
            fstream << ket_comps[i].label() << "," << std::endl;
        }
    }
    else
    {
        const auto bra = Tensor(integral[0]);
        
        const auto bra_comps = bra.components();
        
        for (size_t i = 0; i < bra_comps.size(); i++)
        {
            if (i > 0) fstream << spacer;
            
            fstream << padding << "TDoubleArray& buffer_";
            
            fstream << bra_comps[i].label() << "," << std::endl;
        }
    }
    
    _write_prim_data_decl(fstream, spacer, terminus);
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
   
    const auto spacer = std::string(fname.size() + 1,  ' ');
    
    const auto padding = std::string(6, ' ');
    
    fstream << "auto" << std::endl;
    
    fstream << fname << "(";
    
    const auto op = integral.integrand();
    
    const auto op_comps = op.components();
    
    for (size_t i = 0; i < op_comps.size(); i++)
    {
        if (i > 0) fstream << spacer;
        
        fstream << padding << "TDoubleArray& buffer_";
        
        fstream << op_comps[i].label() << "," << std::endl;
    }
    
    _write_prim_data_decl(fstream, spacer, terminus);
}

void
T2CCPUGenerator::_write_prim_data_decl(      std::ofstream& fstream,
                                       const std::string&   spacer, 
                                       const bool           terminus) const
{
    fstream << spacer << "const double        bra_exp," << std::endl;
    
    fstream << spacer << "const double        bra_norm," << std::endl;
    
    fstream << spacer << "const TPoint3D&     bra_coord," << std::endl;
    
    fstream << spacer << "const TDoubleArray& ket_exps," << std::endl;
    
    fstream << spacer << "const TDoubleArray& ket_norms," << std::endl;
    
    fstream << spacer << "const TDoubleArray& ket_coords_x," << std::endl;
    
    fstream << spacer << "const TDoubleArray& ket_coords_y," << std::endl;
    
    fstream << spacer << "const TDoubleArray& ket_coords_z," << std::endl;
    
    fstream << spacer << "const int64_t       ket_dim) -> void";
    
    if (terminus) fstream << ";" << std::endl;
    
    fstream << std::endl;
}


void
T2CCPUGenerator::_write_gtos_decl(      std::ofstream& fstream,
                                  const bool           diagonal) const
{
    const auto spacer = std::string(4, ' ');
    
    if (diagonal)
    {
        fstream << spacer << "// intialize GTOs data" << std::endl << std::endl;
        
        fstream << spacer << "const auto gto_coords = gto_block.getCoordinates();" << std::endl << std::endl;
       
        fstream << spacer << "const auto gto_exps = gto_block.getExponents();" << std::endl << std::endl;
       
        fstream << spacer << "const auto gto_norms = gto_block.getNormalizationFactors();" << std::endl << std::endl;
        
        fstream << spacer << "const auto gto_indexes = gto_block.getOrbitalIndexes();" << std::endl << std::endl;
       
        fstream << spacer << "const auto ncgtos = gto_block.getNumberOfBasisFunctions();" << std::endl << std::endl;

        fstream << spacer << "const auto npgtos = gto_block.getNumberOfPrimitives();" << std::endl << std::endl;
    }
    else
    {
        fstream << spacer << "// intialize GTOs data on bra side" << std::endl << std::endl;
        
        fstream << spacer << "const auto bra_gto_coords = bra_gto_block.getCoordinates();" << std::endl << std::endl;
       
        fstream << spacer << "const auto bra_gto_exps = bra_gto_block.getExponents();" << std::endl << std::endl;
       
        fstream << spacer << "const auto bra_gto_norms = bra_gto_block.getNormalizationFactors();" << std::endl << std::endl;
        
        fstream << spacer << "const auto bra_gto_indexes = bra_gto_block.getOrbitalIndexes();" << std::endl << std::endl;
       
        fstream << spacer << "const auto bra_ncgtos = bra_gto_block.getNumberOfBasisFunctions();" << std::endl << std::endl;

        fstream << spacer << "const auto bra_npgtos = bra_gto_block.getNumberOfPrimitives();" << std::endl << std::endl;
        
        fstream << spacer << "// intialize GTOs data on ket side" << std::endl << std::endl;
        
        fstream << spacer << "const auto ket_gto_coords = ket_gto_block.getCoordinates();" << std::endl << std::endl;
       
        fstream << spacer << "const auto ket_gto_exps = ket_gto_block.getExponents();" << std::endl << std::endl;
       
        fstream << spacer << "const auto ket_gto_norms = ket_gto_block.getNormalizationFactors();" << std::endl << std::endl;
        
        fstream << spacer << "const auto ket_gto_indexes = ket_gto_block.getOrbitalIndexes();" << std::endl << std::endl;
       
        fstream << spacer << "const auto ket_ncgtos = ket_gto_block.getNumberOfBasisFunctions();" << std::endl << std::endl;

        fstream << spacer << "const auto ket_npgtos = ket_gto_block.getNumberOfPrimitives();" << std::endl << std::endl;
    }
}
