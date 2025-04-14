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

#ifndef g2c_decl_hpp
#define g2c_decl_hpp

#include <string>
#include <vector>
#include <fstream>
#include <utility>

#include "t2c_defs.hpp"

// Two-center functions declaration generator for CPU.
class G2CDeclDriver
{
    /// Generates vector of matrix strings.
    /// @param integral The base two center integral.
    /// @param use_rs The flag for use of range-separated Coulomb interactions.
    /// @return The vector of matrix strings.
    std::vector<std::string> _get_distributor_str(const I2CIntegral& integral,
                                                  const bool         use_rs) const;
        
    /// Generates vector of GTOs block strings.
    /// @param integral The base two center integral.
    /// @param use_rs The flag for use of range-separated Coulomb interactions.
    /// @return The vector of GTOs block strings,
    std::vector<std::string> _get_gto_blocks_str(const I2CIntegral& integral,
                                                 const bool         use_rs) const;
        
    /// Generates vector of indices strings.
    /// @param integral The base two center integral.
    /// @param use_rs The flag for use of range-separated Coulomb interactions.
    /// @param terminus The flag to add termination symbol.
    /// @return The vector of indices strings.
    std::vector<std::string> _get_indices_str(const I2CIntegral& integral,
                                              const bool         use_rs,
                                              const bool         terminus) const;

public:
    /// Creates a two-center functions declaration generator.
    G2CDeclDriver() = default;
    
    /// Writes declaration for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param use_rs The flag for use of range-separated Coulomb interactions.
    /// @param terminus The flag to add termination symbol.
    void write_func_decl(      std::ofstream&         fstream,
                         const I2CIntegral&           integral,
                         const bool                   use_rs,
                         const bool                   terminus) const;
};

#endif /* g2c_decl_hpp */
