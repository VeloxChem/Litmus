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
    /// Gets index of tensor componnent in it's component wwise expansion.
    /// @param component The tensor component to find index.
    /// @return The index of tensor component.
    int _get_tensor_component_index(const TensorComponent& component) const;
    
    /// Combines two symbolic factors into  one.
    /// @param bra_factor The  symbolic factor on bra side.
    /// @param ket_factor The  symbolic factor on ket side. 
    /// @return The combined  factor.
    std::string _combine_factors(const std::string& bra_factor,
                                 const std::string& ket_factor) const;
    
    /// Checks if recursion is available for two-center inetgral with given label.
    /// @param label The label of requested two-center integral.
    bool _is_available(const std::string& label) const;
    
    /// Gets label of standart namespaces for integrand.
    /// @param integrand the integrand operator.
    /// @return The label of namespace.
    std::string _get_namespace_label(const Operator& integrand) const;
    
    /// Gets symmetry  of matrix  for integrand.
    /// @param integrand the integrand operator.
    /// @return The symmetry of integrand.
    std::string _get_matrix_symmetry(const Operator& integrand) const;
    
    /// Select integrals components with predefined bra or ket side components.
    /// @param component the tensor component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    VT2CIntegrals _select_integral_components(const TensorComponent& component,
                                              const I2CIntegral&     integral,
                                              const bool             bra_first) const;
    /// Writes documentation string for primitive compute function.
    /// @param bra_component the tensor component.
    /// @param ket_component the tensor component.
    /// @param integral The base two center integral.
    VT2CIntegrals _select_integral_components(const TensorComponent& bra_component,
                                              const TensorComponent& ket_component,
                                              const I2CIntegral&     integral) const;
    
    /// Gets two-center inetgral with requested label.
    /// @param label The label of requested two-center integral.
    /// @param ang_a The angular momentum of center A.
    /// @param ang_b The angular momentum of center B.
    /// @return The two-center integral.
    I2CIntegral _get_integral(const std::string& label,
                              const int          ang_a,
                              const int          ang_b) const;
    
    /// Gets recursion factors label.
    /// @param rterm The recursion term.
    /// @param first The position of recursion term in code line.
    std::string _get_factor_label(const R2CTerm& rterm,
                                  const bool     first) const;
    /// Checks if factors is needed by recursion group.
    /// @param rgroup The recursion group.
    /// @param label The label of factor to find.
    /// @return True if factor is found, False otherwise.
    bool _find_factor(const R2Group&     rgroup,
                      const std::string& label) const;
    
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
        
    /// Writes primitive functions documentation and declaration to header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_prim_funcs_to_cpp_header(      std::ofstream& fstream,
                                         const I2CIntegral&   integral) const;
    
    /// Writes primitive function implementation to C++ code file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_prim_funcs_to_cpp_file(      std::ofstream& fstream,
                                       const I2CIntegral&   integral) const;
    
    /// Writes documentation string for primitive data.
    /// @param fstream the file stream.
    void _write_prim_data_docstr(std::ofstream& fstream) const;
    
    /// Writes declaration for GTOs data in compute function.
    /// @param fstream the file stream.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_gtos_decl(      std::ofstream& fstream,
                          const bool           diagonal) const;
    
    /// Writes declaration for  spherical momentum factors in compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_angmom_decl(      std::ofstream& fstream,
                            const I2CIntegral&   integral) const;
    
    /// Writes declaration for ket data in compute function.
    /// @param fstream the file stream.
    void _write_ket_data_decl(std::ofstream& fstream) const;
    
    /// Writes declaration of buffers for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_buffers_decl(      std::ofstream&   fstream,
                             const I2CIntegral&     integral) const;
    
    /// Writes declaration for main batches loop start in compute function.
    /// @param fstream the file stream.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_batches_loop_start_decl(      std::ofstream& fstream,
                                        const bool           diagonal) const;
    
    /// Writes declaration for main batches loop end in compute function.
    /// @param fstream the file stream.
    void _write_batches_loop_end_decl(std::ofstream& fstream) const;
    
    /// Writes declaration of main call tree for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_main_call_tree_decl(      std::ofstream& fstream,
                                    const I2CIntegral&   integral,
                                    const bool           diagonal) const;
    
    /// Writes declaration of call tree block for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_prim_call_tree_block_decl(      std::ofstream& fstream,
                                          const I2CIntegral&   integral,
                                          const bool           diagonal) const;
    
    /// Writes declaration of call tree block for compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_prim_call_tree_block_decl(      std::ofstream&   fstream,
                                          const TensorComponent& component,
                                          const I2CIntegral&     integral,
                                          const bool             bra_first,
                                          const bool             diagonal) const;
    
    /// Writes declaration of call tree block for compute function.
    /// @param fstream the file stream.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_prim_call_tree_block_decl(      std::ofstream&   fstream,
                                          const TensorComponent& bra_component,
                                          const TensorComponent& ket_component,
                                          const I2CIntegral&     integral,
                                          const bool             diagonal) const;
    
    /// Writes declaration for primitives loop start in compute function.
    /// @param fstream the file stream.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_primitives_loop_start_decl(      std::ofstream& fstream,
                                          const bool           diagonal) const;
    
    /// Writes declaration for primitives loop end in compute function.
    /// @param fstream the file stream.
    void _write_primitives_loop_end_decl(std::ofstream& fstream) const;
    
    /// Writes declaration of common variables for primitive computation function.
    /// @param fstream the file stream.
    /// @param spacer The size of spacer for formatting variables.
    void _write_primitives_call_data_decl(      std::ofstream& fstream,
                                          const size_t         spacer) const;
    
    /// Writes declaration of block distribution call tree for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    ///@param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_block_distributor_decl(      std::ofstream& fstream,
                                       const I2CIntegral&   integral,
                                       const bool           diagonal) const;
    
    /// Writes declaration of block distribution call tree for compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_block_distributor_decl(      std::ofstream&   fstream,
                                       const TensorComponent& component,
                                       const I2CIntegral&     integral,
                                       const bool             bra_first,
                                       const bool             diagonal) const;
    
    /// Writes declaration of of block distribution call tree  for compute function.
    /// @param fstream the file stream.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_block_distributor_decl(      std::ofstream&   fstream,
                                       const TensorComponent& bra_component,
                                       const TensorComponent& ket_component,
                                       const I2CIntegral&     integral,
                                       const bool             diagonal) const;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_prim_func_body(      std::ofstream& fstream,
                                const I2CIntegral&   integral) const;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    void _write_prim_func_body(      std::ofstream&   fstream,
                               const TensorComponent& component,
                               const I2CIntegral&     integral,
                               const bool             bra_first) const;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    void _write_prim_func_body(      std::ofstream&   fstream,
                               const TensorComponent& bra_component,
                               const TensorComponent& ket_component,
                               const I2CIntegral&     integral) const;
    
    /// Writes common data definitions in primitive compute function.
    /// @param fstream the file stream.
    void _write_prim_func_common_data(std::ofstream& fstream) const;
    
    /// Writes buffers for primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_prim_func_buffers(      std::ofstream& fstream,
                                  const I2CIntegral&   integral) const;
    
    /// Writes buffers for primitive compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    void _write_prim_func_buffers(      std::ofstream&   fstream,
                                  const TensorComponent& component,
                                  const I2CIntegral&     integral,
                                  const bool             bra_first) const;
    
    /// Writes buffers for primitive compute function.
    /// @param fstream the file stream.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    void _write_prim_func_buffers(      std::ofstream&   fstream,
                                  const TensorComponent& bra_component,
                                  const TensorComponent& ket_component,
                                  const I2CIntegral&     integral) const;
    
    /// Writes main loop pragma for primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_prim_func_pragma(      std::ofstream& fstream,
                                 const I2CIntegral&   integral) const;
    
    /// Writes main loop pragma for primitive compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    void _write_prim_func_pragma(      std::ofstream&   fstream,
                                 const TensorComponent& component,
                                 const I2CIntegral&     integral,
                                 const bool             bra_first) const;
    
    /// Writes main loop pragma for primitive compute function.
    /// @param fstream the file stream.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    void _write_prim_func_pragma(      std::ofstream&   fstream,
                                 const TensorComponent& bra_component,
                                 const TensorComponent& ket_component,
                                 const I2CIntegral&     integral) const;
    
    /// Writes common pragma variables for primitive compute function.
    /// @param fstream the file stream.
    void _write_prim_func_common_pragma(std::ofstream& fstream) const;
    
    /// Writes main loop start for primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_prim_func_loop_start(      std::ofstream& fstream,
                                     const I2CIntegral&    integral) const;
    
    /// Writes main loop end for primitive compute function.
    /// @param fstream the file stream.
    void _write_prim_func_loop_end(std::ofstream& fstream) const;
    
    /// Writes simd code generated for selected set of integral components.
    /// @param fstream the file stream.
    /// @param labels the vector of labels for integral components.
    /// @param components the vector of integral component.
    /// @param integral the base integral.
    void _write_simd_code(      std::ofstream&            fstream,
                          const std::vector<std::string>& labels,
                          const VT2CIntegrals&            components,
                          const I2CIntegral&              integral) const;
    
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
