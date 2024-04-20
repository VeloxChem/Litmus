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
    
    /// Generates vector of matrix strings.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @return The vector of matrix strings.
    std::vector<std::string> _get_matrices_str(const I2CIntegral& integral,
                                               const std::pair<bool, bool>& rec_form) const;
    
    /// Generates vector of GTOs block strings.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of GTOs block strings,
    std::vector<std::string> _get_gto_blocks_str(const I2CIntegral& integral,
                                                 const bool         diagonal) const;
    /// Generates vector of distributor variables strings.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of distributor variables strings
    std::vector<std::string> _get_distributor_variables_str(const I2CIntegral& integral,
                                                            const bool         diagonal) const;
    
    /// Generates vector of indices strings.
    /// @return The vector of indices strings.
    std::vector<std::string> _get_indices_str() const;
    
public:
    /// Creates a two-center documentation generator.
    T2CDocuDriver() = default;
    
    /// Writes documentation string for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void write_doc_str(      std::ofstream&         fstream,
                       const I2CIntegral&           integral,
                       const std::pair<bool, bool>& rec_form,
                       const bool                   diagonal) const;
    
};

#endif /* t2c_docs_hpp */