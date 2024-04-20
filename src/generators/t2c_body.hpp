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

#ifndef t2c_body_hpp
#define t2c_body_hpp

#include <string>
#include <vector>
#include <utility>
#include <fstream>

#include "t2c_defs.hpp"
#include "file_stream.hpp"

// Two-center compute function body generators for CPU.
class T2CFuncBodyDriver
{
    /// Generates vector of strings with GTOs definitions in compute function.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of strings with GTOS definitions in compute function.
    std::vector<std::string> _get_gtos_def(const bool diagonal) const;
    
    /// Generates vector of ket factors in compute function.
    /// @return The vector of ket factors in compute function.
    std::vector<std::string> _get_ket_variables_def(const bool diagonal) const;
    
    /// Generates vector of distances in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of distances in compute function.
    std::vector<std::string> _get_coordinates_def(const I2CIntegral& integral) const;
    
    /// Generates vector of distances in compute function.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    /// @return The vector of distances in compute function.
    std::vector<std::string> _get_buffers_def(const SI2CIntegrals& integrals,
                                              const I2CIntegral&   integral) const;
    
    /// Adds loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_loop_start(      VCodeLines&  lines,
                         const I2CIntegral& integral,
                         const bool         diagonal) const;
    
    /// Adds loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_loop_end(      VCodeLines&  lines,
                       const I2CIntegral& integral,
                       const bool         diagonal) const;
    
    /// Adds ket loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_ket_loop_start(      VCodeLines&  lines,
                             const I2CIntegral& integral,
                             const bool         diagonal) const;
    
    /// Adds ket loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_ket_loop_end(      VCodeLines&  lines,
                           const I2CIntegral& integral,
                           const bool         diagonal) const;
    
    /// Adds auxilary integrals.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    void _add_auxilary_integrals(      VCodeLines&  lines,
                                 const SI2CIntegrals& integrals) const;
    
    /// Adds call tree for recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    void _add_call_tree(      VCodeLines&  lines,
                        const SI2CIntegrals& integrals) const;
    
    /// Gets arguments list for primitive function call.
    /// @param integral The base two center integral.
    std::string _get_arguments(const I2CIntegral& integral) const;
    
public:
    /// Creates a two-center compute function body generator.
    T2CFuncBodyDriver() = default;
    
    /// Writes body of compute function.
    /// @param fstream the file stream.
    /// @param integrals The set of integrals.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void write_func_body(      std::ofstream&         fstream,
                         const SI2CIntegrals&         integrals,
                         const I2CIntegral&           integral,
                         const std::pair<bool, bool>& rec_form,
                         const bool                   diagonal) const;
};

#endif /* t2c_body_hpp */
