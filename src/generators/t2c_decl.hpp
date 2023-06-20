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

#ifndef t2c_decl_hpp
#define t2c_decl_hpp

#include <string>
#include <vector>
#include <fstream>
#include <utility>

#include "t2c_defs.hpp"

// Two-center functions declaration generator for CPU.
class T2CDeclDriver
{
    /// Generates vector of matrix strings.
    /// @param integral The base two center integral.
    /// @return The vector of matrix strings.
    std::vector<std::string> _get_matrix_str(const I2CIntegral& integral) const;
    
    /// Generates vector of GTOs block strings.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of GTOs block strings
    std::vector<std::string> _get_gto_blocks_str(const I2CIntegral& integral,
                                                 const bool         diagonal) const;
    
    /// Generates vector of index strings.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @param terminus The flag to add termination symbol.
    /// @return The vector of index strings.
    std::vector<std::string> _get_indexes_str(const I2CIntegral& integral,
                                              const bool         diagonal,
                                              const bool         terminus) const;
    
    /// Generates vector of primitive buffer strings.
    /// @param integral The base two center integral.
    /// @param terminus The flag to add termination symbol.
    /// @return The vector of primitive buffer strings.
    std::vector<std::string> _get_prim_buffer_str(const I2CIntegral& integral,
                                                  const bool         terminus) const;
    
    /// Generates vector of primitive buffer strings.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    /// @param terminus The flag to add termination symbol.
    /// @return The vector of primitive buffer strings.
    std::vector<std::string> _get_prim_buffer_str(const TensorComponent& component,
                                                  const I2CIntegral&     integral,
                                                  const bool             bra_first,
                                                  const bool             terminus) const;
    
    /// Generates vector of primitive buffer strings.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    /// @param terminus The flag to add termination symbol.
    /// @return The vector of primitive buffer strings.
    std::vector<std::string> _get_prim_buffer_str(const TensorComponent& bra_component,
                                                  const TensorComponent& ket_component,
                                                  const I2CIntegral&     integral,
                                                  const bool             terminus) const;
    
    /// Adds common variable strings to vector of strings.
    /// @param vstrings  The vector of strings to addd common variable strings.
    /// @param spacer The padding spacer for  common  variables.
    /// @param terminus The flag to add termination symbol.
    void _add_prim_variables(      std::vector<std::string>& vstrings,
                             const std::string&              spacer,
                             const bool                      terminus) const;
    
public:
    /// Creates a two-center functions declaration generator.
    T2CDeclDriver() = default;
    
    /// Writes declaration for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @param terminus The flag to add termination symbol.
    void write_func_decl(      std::ofstream& fstream,
                         const I2CIntegral&   integral,
                         const bool           diagonal,
                         const bool           terminus) const;
    
    /// Writes declaration of primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param terminus The flag to add termination symbol.
    void write_prim_func_decl(      std::ofstream& fstream,
                              const I2CIntegral&   integral,
                              const bool           terminus) const;
    
    /// Writes declaration of primitive compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    /// @param terminus The flag to add termination symbol.
    void write_prim_func_decl(      std::ofstream&   fstream,
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
    void write_prim_func_decl(      std::ofstream&   fstream,
                              const TensorComponent& bra_component,
                              const TensorComponent& ket_component,
                              const I2CIntegral&     integral,
                              const bool             terminus) const;
    
};

#endif /* t2c_decl_hpp */
