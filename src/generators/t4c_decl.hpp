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

#ifndef t4c_decl_hpp
#define t4c_decl_hpp

#include <string>
#include <vector>
#include <fstream>
#include <utility>

#include "t4c_defs.hpp"

// Four-center functions declaration generator for CPU.
class T4CDeclDriver
{
    /// Generates vector of matrix strings.
    /// @param integral The base two center integral.
    /// @return The vector of matrix strings.
    std::vector<std::string> _get_matrices_str(const I4CIntegral& integral) const;
    
    /// Generates vector of GTOs block strings.
    /// @param integral The base two center integral.
    /// @return The vector of GTOs block strings,
    std::vector<std::string> _get_gto_pair_blocks_str(const I4CIntegral& integral) const;
        
    /// Generates vector of indices strings.
    /// @param integral The base two center integral.
    /// @param terminus The flag to add termination symbol.
    /// @return The vector of indices strings.
    std::vector<std::string> _get_indices_str(const I4CIntegral& integral,
                                              const bool         terminus) const;
    
    /// Generates vector of matrix strings.
    /// @param integral The base two center integral.
    /// @return The vector of matrix strings.
    std::vector<std::string> _get_diag_matrices_str(const I4CIntegral& integral) const;
    
    /// Generates vector of GTOs block strings.
    /// @param integral The base two center integral.
    /// @return The vector of GTOs block strings,
    std::vector<std::string> _get_diag_gto_pair_blocks_str(const I4CIntegral& integral) const;
        
    /// Generates vector of indices strings.
    /// @param integral The base two center integral.
    /// @param terminus The flag to add termination symbol.
    /// @return The vector of indices strings.
    std::vector<std::string> _get_diag_indices_str(const I4CIntegral& integral,
                                                   const bool         terminus) const;

public:
    /// Creates a four-center functions declaration generator.
    T4CDeclDriver() = default;
    
    /// Writes declaration for compute function.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    /// @param terminus The flag to add termination symbol.
    void write_func_decl(      std::ofstream& fstream,
                         const I4CIntegral&   integral,
                         const bool           terminus) const;
    
    /// Writes declaration for compute function.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    /// @param terminus The flag to add termination symbol.
    void write_diag_func_decl(      std::ofstream& fstream,
                              const I4CIntegral&   integral,
                              const bool           terminus) const;
};


#endif /* t4c_decl_hpp */
