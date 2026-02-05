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

#ifndef t2c_ecp_prim_body_hpp
#define t2c_ecp_prim_body_hpp

#include <string>
#include <array>
#include <vector>
#include <utility>
#include <fstream>

#include "t2c_defs.hpp"
#include "file_stream.hpp"

// Two-center local ECP compute function body generators for CPU.
class T2CECPPrimFuncBodyDriver
{
    /// Generates vector of factor strings.
    /// @param integral The base two center integral.
    /// @return The vector of factor strings.
    std::vector<std::string> _get_factors_str(const I2CIntegral& integral) const;
    
    /// Computes VRR recursion for integral component.
    /// @param integral The base two center integral component.
    /// @return The recursion expansion of integral component.
    R2CDist _get_vrr_recursion(const T2CIntegral& integral) const;
    
public:
    /// Creates a two-center compute function body generator.
    T2CECPPrimFuncBodyDriver() = default;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void write_func_body(      std::ofstream& fstream,
                         const I2CIntegral&   integral) const;
};

#endif /* t2c_ecp_prim_body_hpp */
