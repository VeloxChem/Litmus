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
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @return The vector of matrix strings.
    std::vector<std::string> _get_matrices_str(const I2CIntegral&           integral,
                                               const std::pair<bool, bool>& rec_form) const;
    
    /// Generates vector of special variables strings.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @return The vector of special variables strings.
    std::vector<std::string> _get_special_variables_str(const I2CIntegral& integral,
                                                        const std::pair<bool, bool>& rec_form) const;
    
    /// Generates vector of GTOs block strings.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of GTOs block strings,
    std::vector<std::string> _get_gto_blocks_str(const I2CIntegral&           integral,
                                                 const std::pair<bool, bool>& rec_form,
                                                 const bool                   diagonal) const;
    
    /// Generates vector of distributor variables strings.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of distributor variables strings
    std::vector<std::string> _get_distributor_variables_str(const I2CIntegral&           integral,
                                                            const std::pair<bool, bool>& rec_form,
                                                            const bool                   diagonal) const;
    
    /// Generates vector of indices strings.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @param terminus The flag to add termination symbol.
    /// @return The vector of indices strings.
    std::vector<std::string> _get_indices_str(const I2CIntegral&           integral,
                                              const std::pair<bool, bool>& rec_form,
                                              const bool                   diagonal,
                                              const bool                   terminus) const;

public:
    /// Creates a two-center functions declaration generator.
    T2CDeclDriver() = default;
    
    /// Writes declaration for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @param terminus The flag to add termination symbol.
    void write_func_decl(      std::ofstream&         fstream,
                         const I2CIntegral&           integral,
                         const std::pair<bool, bool>& rec_form,
                         const bool                   diagonal,
                         const bool                   terminus) const;
};

#endif /* t2c_decl_hpp */
