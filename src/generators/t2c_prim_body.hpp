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

#ifndef t2c_prim_body_hpp
#define t2c_prim_body_hpp

#include <string>
#include <array>
#include <vector>
#include <utility>
#include <fstream>

#include "t2c_defs.hpp"
#include "file_stream.hpp"

// Two-center compute function body generators for CPU.
class T2CPrimFuncBodyDriver
{
    /// Generates vector of buffer strings.
    /// @param integral The base two center integral.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_buffers_str(const I2CIntegral& integral) const;
    
    /// Gets tensor label for integral.
    /// @param integral The base two center integral.
    /// @return The tensorial label.
    std::string _get_tensor_label(const I2CIntegral& integral) const;
    
    /// Adds single loop computation of primitive integrals.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    /// @param components The vector of integral components.
    /// @param rec_range The recursion range [first, last) in integral components space.
    void _add_recursion_loop(      VCodeLines&         lines,
                             const I2CIntegral&        integral,
                             const VT2CIntegrals&      components,
                             const std::array<int, 2>& rec_range) const;
    
    std::string _get_pragma(const I2CIntegral&   integral,
                            const VT2CIntegrals& components) const;

public:
    /// Creates a two-center compute function body generator.
    T2CPrimFuncBodyDriver() = default;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void write_func_body(      std::ofstream& fstream,
                         const I2CIntegral&   integral) const;
};


#endif /* t2c_prim_body_hpp */
