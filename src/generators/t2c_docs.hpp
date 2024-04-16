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

#ifndef t2c_docs_hpp
#define t2c_docs_hpp

#include <string>
#include <vector>
#include <fstream>

#include "t2c_defs.hpp"

// Two-center documentation generator for CPU.
class T2CDocuDriver
{
    /// Generates compute string.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The compute string.
    std::string _get_compute_str(const I2CIntegral& integral,
                                 const bool         diagonal) const;
    
    /// Generates auxilary compute string.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The auxilary compute string.
    std::string _get_auxilary_compute_str(const I2CIntegral& integral,
                                          const bool         diagonal) const;
    
    /// Generates primitive compute string.
    /// @param integral The base two center integral.
    /// @return The primitive compute string.
    std::string _get_prim_compute_str(const I2CIntegral& integral) const;
    
    /// Generates primitive compute string.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    /// @return The primitive compute string.
    std::string _get_prim_compute_str(const TensorComponent& component,
                                      const I2CIntegral&     integral,
                                      const bool             bra_first) const;
    
    /// Generates primitive compute string.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    /// @return The primitive compute string.
    std::string _get_prim_compute_str(const TensorComponent& bra_component,
                                      const TensorComponent& ket_component,
                                      const I2CIntegral&     integral) const;
    
    /// Generates vector of matrix strings.
    /// @param integral The base two center integral.
    /// @return The vector of matrix strings.
    std::vector<std::string> _get_matrix_str(const I2CIntegral& integral) const;
    
    /// Generates vector of special variable strings.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @return The vector of special variable strings.
    std::vector<std::string> _get_special_vars_str(const I2CIntegral& integral,
                                                   const bool         sum_form) const;
    
    /// Generates vector of GTOs block strings.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of GTOs block strings
    std::vector<std::string> _get_gto_blocks_str(const I2CIntegral& integral,
                                                 const bool         is_auxilary,
                                                 const bool         diagonal) const;
    
    /// Generates vector of index strings.
    /// @return The vector of index strings.
    std::vector<std::string> _get_indexes_str() const;
    
    /// Generates vector of auxilary index strings.
    /// @return The vector of index strings.
    std::vector<std::string> _get_auxilary_indexes_str() const;
    
    /// Generates matrix type string.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The matrix type string.
    std::string _get_matrix_type_str(const I2CIntegral& integral,
                                     const bool         diagonal) const;
    
    /// Generates vector of primitive buffer strings.
    /// @param integral The base two center integral.
    /// @return The vector of primitive buffer strings.
    std::vector<std::string> _get_prim_buffer_str(const I2CIntegral& integral) const;
    
    /// Generates vector of primitive buffer strings.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    /// @return The vector of primitive buffer strings.
    std::vector<std::string> _get_prim_buffer_str(const I2CIntegral& integral,
                                                  const bool         bra_first) const;
    
    /// Generates vector of primitive common variable strings.
    /// @return The vector of primitive  common variable strings.
    std::vector<std::string> _get_prim_variables_str() const;
    
    /// Generates vector of primitive common variable strings.
    /// @return The vector of primitive  common variable strings.
    std::vector<std::string> _get_prim_variables_str(const I2CIntegral& integral) const;
    
public:
    /// Creates a two-center documentation generator.
    T2CDocuDriver() = default;
    
    /// Writes documentation string for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void write_doc_str(      std::ofstream& fstream,
                       const I2CIntegral&   integral,
                       const bool           sum_form,
                       const bool           diagonal) const;
    
    /// Writes documentation string for primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void write_prim_doc_str(      std::ofstream& fstream,
                            const I2CIntegral&    integral) const;
    
    /// Writes documentation string for auxilary compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void write_auxilary_doc_str(      std::ofstream& fstream,
                                const I2CIntegral&   integral,
                                const bool           diagonal) const;
    
    /// Writes documentation string for primtive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void write_prim_doc_str(      std::ofstream& fstream,
                            const I2CIntegral&   integral,
                            const bool           sum_form) const;
    
    /// Writes documentation string for primitive compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @param bra_first The flag to set bra as expansion point.
    void write_prim_doc_str(      std::ofstream&   fstream,
                            const TensorComponent& component,
                            const I2CIntegral&     integral,
                            const bool             sum_form,
                            const bool             bra_first) const;
    
    /// Writes documentation string for primitive compute function.
    /// @param fstream the file stream.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void write_prim_doc_str(      std::ofstream&   fstream,
                            const TensorComponent& bra_component,
                            const TensorComponent& ket_component,
                            const I2CIntegral&     integral,
                            const bool             sum_form) const;
};

#endif /* t2c_docs_hpp */
