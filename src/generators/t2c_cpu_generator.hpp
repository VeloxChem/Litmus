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
#include "file_stream.hpp"

// Two-center integrals code generator for CPU.
class T2CCPUGenerator
{
    /// Checks if recursion is available for two-center inetgral with given label.
    /// @param label The label of requested two-center integral.
    bool _is_available(const std::string& label) const;
    
    /// Gets two-center inetgral with requested label.
    /// @param label The label of requested two-center integral.
    /// @param ang_a The angular momentum of center A.
    /// @param ang_b The angular momentum of center B.
    /// @param bra_gdrv The geometrical derivative of bra side.
    /// @param ket_gdrv The geometrical derivative of ket side.
    /// @param op_gdrv The geometrical derivative of operator. 
    /// @return The two-center integral.
    I2CIntegral _get_integral(const std::string& label,
                              const int          ang_a,
                              const int          ang_b,
                              const int          bra_gdrv = 0,
                              const int          ket_gdrv = 0,
                              const int          op_gdrv  = 0) const;
    
    /// Gets file name of file with recursion functions for two center integral.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @return The file name.
    std::string _file_name(const I2CIntegral& integral,
                           const bool         sum_form) const;
    
    /// Writes header file for recursion.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _write_cpp_header(const I2CIntegral& integral,
                           const bool         sum_form) const;
    
    /// Writes C++ code file for recursion.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _write_cpp_file(const I2CIntegral& integral,
                         const bool         sum_form) const;
    
    /// Writes header files for primitive recursion.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _write_cpp_prim_headers(const I2CIntegral& integral,
                                 const bool         sum_form) const;
    
    /// Writes C++ code files for primitive recursion.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _write_cpp_prim_files(const I2CIntegral& integral,
                               const bool         sum_form) const;
    
    /// Writes header files for specific primitive recursion.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _write_cpp_prim_header(const I2CIntegral& integral,
                                const bool         sum_form) const;
    
    /// Writes header files for specific primitive recursion.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @param bra_first The flag to set bra as expansion point.
    void _write_cpp_prim_header(const TensorComponent& component,
                                const I2CIntegral&     integral,
                                const bool             sum_form,
                                const bool             bra_first) const;
    
    /// Writes header files for specific primitive recursion.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _write_cpp_prim_header(const TensorComponent& bra_component,
                                const TensorComponent& ket_component,
                                const I2CIntegral&     integral,
                                const bool             sum_form) const;
    
    /// Writes C++ code files for specific primitive recursion.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _write_cpp_prim_file(const I2CIntegral& integral,
                              const bool         sum_form) const;
    
    /// Writes C++ code files for specific primitive recursion.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @param bra_first The flag to set bra as expansion point.
    void _write_cpp_prim_file(const TensorComponent& component,
                              const I2CIntegral&     integral,
                              const bool             sum_form,
                              const bool             bra_first) const;
    
    /// Writes C++ code files for specific primitive recursion.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _write_cpp_prim_file(const TensorComponent& bra_component,
                              const TensorComponent& ket_component,
                              const I2CIntegral&     integral,
                              const bool             sum_form) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_defines(      std::ofstream& fstream,
                            const I2CIntegral&   integral,
                            const bool           sum_form,
                            const bool           start) const;
    
    /// Writes definitions of define for primitive header file.
    /// @param fstream the file stream.
    /// @param fname The base file name.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_prim_defines(      std::ofstream& fstream,
                                 const std::string&   fname,
                                 const bool           start) const;
    
    /// Writes definitions of includes for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _write_hpp_includes(      std::ofstream& fstream,
                             const I2CIntegral&   integral,
                             const bool           sum_form) const;
    
    /// Writes definitions of includes for primitives header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _write_hpp_prim_includes(      std::ofstream& fstream,
                                  const I2CIntegral&   integral,
                                  const bool           sum_form) const;
    
    /// Writes definitions of includes for C++ code file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _write_cpp_includes(      std::ofstream& fstream,
                             const I2CIntegral&   integral,
                             const bool           sum_form) const;
    
    /// Adds primitive functions headers to code lines container.
    /// @param lines The code lines container to which primitives loop start definition are added.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _add_prim_call_includes(      VCodeLines&  lines,
                                 const I2CIntegral& integral,
                                 const bool         sum_form) const;
    
    /// Writes definitions of includes for C++ code file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_cpp_prim_includes(      std::ofstream& fstream,
                                  const I2CIntegral&   integral) const;
    
    /// Writes definitions of includes for C++ code file.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    void _write_cpp_prim_includes(      std::ofstream&   fstream,
                                  const TensorComponent& component,
                                  const I2CIntegral&     integral,
                                  const bool             bra_first) const;
    
    /// Writes definitions of includes for C++ code file.
    /// @param fstream the file stream.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    void _write_cpp_prim_includes(      std::ofstream&   fstream,
                                  const TensorComponent& bra_component,
                                  const TensorComponent& ket_component,
                                  const I2CIntegral&     integral) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const I2CIntegral&   integral,
                          const bool           start) const;
        
    /// Writes primitive functions documentation and declaration to header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _write_prim_funcs_to_cpp_header(      std::ofstream& fstream,
                                         const I2CIntegral&   integral,
                                         const bool           sum_form) const;
public:
    /// Creates an electron repulsion integrals CPU code generator.
    T2CCPUGenerator() = default;
     
    /// Generates selected one-electron integrals up to given angular momentum (inclusive) )on A and B centers.
    /// @param label The label of requested two-center integral.
    /// @param angmom The maximum angular momentum of A and B centers.
    /// @param bra_gdrv The geometrical derivative of bra side.
    /// @param ket_gdrv The geometrical derivative of ket side.
    /// @param op_gdrv The geometrical derivative of operator.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void generate(const std::string& label,
                  const int          angmom,
                  const int          bra_gdrv,
                  const int          ket_gdrv,
                  const int          op_gdrv,
                  const bool         sum_form) const;
};

#endif /* t2c_cpu_generator_hpp */
