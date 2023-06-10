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

#ifndef t2c_cpu_generator_hpp
#define t2c_cpu_generator_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>


#include "operator.hpp"
#include "tensor_component.hpp"
#include "t2c_defs.hpp"

// Two-center integrals code generator for CPU.
class T2CCPUGenerator
{
    /// Checks if recursion is available for two-center inetgral with given label.
    /// @param label The label of requested two-center integral.
    bool _is_available(const std::string& label) const;
    
    /// Gets standart capitalized label of integral.
    /// @param integral The base two center integral.
    /// @return The standart capitalized label of integral.
    std::string _get_label(const I2CIntegral& integral) const;
    
    /// Gets map of standart integrand labels.
    /// @return The map of standart integrad labels.
    std::map<Operator, std::string> _get_integrands_map() const;
    
    /// Gets two-center inetgral with requested label.
    /// @param label The label of requested two-center integral.
    /// @param ang_a The angular momentum of center A.
    /// @param ang_b The angular momentum of center B.
    /// @return The two-center integral.
    I2CIntegral _get_integral(const std::string& label,
                              const int          ang_a,
                              const int          ang_b) const;
    
    /// Gets file name of file with recursion functions for two center integral.
    /// @param integral The base two center integral.
    /// @return The file name.
    std::string _file_name(const I2CIntegral& integral) const;
    
    /// Writes header file for recursion.
    /// @param integral The base four center integral.
    void _write_cpp_header(const I2CIntegral& integral) const;
    
    /// Writes C++ code file for recursion.
    /// @param integral The base four center integral.
    void _write_cpp_file(const I2CIntegral& integral) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_defines(      std::ofstream& fstream,
                            const I2CIntegral&   integral,
                            const bool           start) const;
    
    /// Writes definitions of includes for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_hpp_includes(      std::ofstream& fstream,
                             const I2CIntegral&   integral) const;
    
    /// Writes definitions of includes for C++ code file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_cpp_includes(      std::ofstream& fstream,
                             const I2CIntegral&   integral) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const I2CIntegral&   integral,
                          const bool           start) const;
    
    /// Writes documentation string for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_func_docstr(      std::ofstream& fstream,
                            const I2CIntegral&   integral,
                            const bool           diagonal) const;
    
    /// Writes documentation strings for all matrices associates with integral.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_matrix_docstr(      std::ofstream& fstream,
                              const I2CIntegral&   integral) const;
    
    /// Writes declaration for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @param terminus The flag to add termination symbol.
    void _write_func_decl(      std::ofstream& fstream,
                          const I2CIntegral&   integral,
                          const bool           diagonal,
                          const bool           terminus) const;
    
    /// Writes documentation strings for all matrices associates with integral.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param spacer The spacer for formatting matrices.
    void _write_matrix_decl(      std::ofstream& fstream,
                            const I2CIntegral&   integral,
                            const std::string&   spacer) const;
    
    /// Writes primitive functions documentation and declaration to header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_prim_funcs_to_cpp_header(      std::ofstream& fstream,
                                         const I2CIntegral&   integral) const;
    
    /// Writes documentation string for primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_prim_func_docstr(      std::ofstream& fstream,
                                 const I2CIntegral&   integral) const;
    
    /// Writes documentation string for primitive compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    void _write_prim_func_docstr(      std::ofstream&   fstream,
                                 const TensorComponent& component,
                                 const I2CIntegral&     integral,
                                 const bool             bra_first) const;
    
    /// Writes documentation string for primitive compute function.
    /// @param fstream the file stream.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    void _write_prim_func_docstr(      std::ofstream&   fstream,
                                 const TensorComponent& bra_component,
                                 const TensorComponent& ket_component,
                                 const I2CIntegral&     integral) const;
    
    /// Writes documentation string for primitive data.
    /// @param fstream the file stream.
    void _write_prim_data_docstr(std::ofstream& fstream) const;
    
    /// Writes declaration of primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param terminus The flag to add termination symbol.
    void _write_prim_func_decl(      std::ofstream& fstream,
                               const I2CIntegral&   integral,
                               const bool           terminus) const;
    
    /// Writes declaration of primitive compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    /// @param terminus The flag to add termination symbol.
    void _write_prim_func_decl(      std::ofstream&   fstream,
                               const TensorComponent& component,
                               const I2CIntegral&     integral,
                               const bool             bra_first,
                               const bool             terminus) const;
    
    /// Writes declatation of primitive compute function.
    /// @param fstream the file stream.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    /// @param terminus The flag to add termination symbol.
    void _write_prim_func_decl(      std::ofstream&   fstream,
                               const TensorComponent& bra_component,
                               const TensorComponent& ket_component,
                               const I2CIntegral&     integral,
                               const bool             terminus) const;
    
    
    /// Writes declaration of primitive data.
    /// @param fstream the file stream.
    /// @param spacer The spacer for formatting variables. 
    /// @param terminus The flag to add termination symbol.
    void _write_prim_data_decl(      std::ofstream& fstream,
                               const std::string&   spacer,
                               const bool           terminus) const;
    
    /// Writes declaration for GTOs dta in compute function.
    /// @param fstream the file stream.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_gtos_decl(      std::ofstream& fstream,
                          const bool           diagonal) const;
        
public:
    /// Creates an electron repulsion integrals CPU code generator.
    T2CCPUGenerator() = default;
     
    /// Generates selected one-electron integrals up to given angular momentum (inclusive) )on A and B centers.
    /// @param label The label of requested two-center integral.
    /// @param angmom The maximum angular momentum of A and B centers.
    void generate(const std::string& label,
                  const int          angmom) const;
};

#endif /* t2c_cpu_generator_hpp */
