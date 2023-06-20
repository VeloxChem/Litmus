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
#include <fstream>

#include "t2c_defs.hpp"
#include "file_stream.hpp"

// Two-center compute function body generators for CPU.
class T2CFuncBodyDriver
{
    /// Generates vector of strings with spherical momentum factors in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of strings with spherical momentum factors in compute function.
    std::vector<std::string> _get_angmom_def(const I2CIntegral& integral) const;
    
    /// Generates vector of strings with GTOs definitions in compute function.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of strings with GTOS definitions in compute function.
    std::vector<std::string> _get_gtos_def(const bool diagonal) const;
    
    /// Generates vector of ket factors in compute function.
    /// @return The vector of ket factors in compute function.
    std::vector<std::string> _get_ket_variables_def() const;
    
    /// Generates vector of buffer definitions in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of buffer definitions in compute function.
    std::vector<std::string> _get_buffers_def(const I2CIntegral& integral) const;
    
    /// Generates vector of strings with main loop definition in compute function.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of strings with main loop definition in compute function.
    std::vector<std::string> _get_batches_def(const bool diagonal) const;
    
    /// Adds loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    void _add_batches_loop_start(VCodeLines& lines) const;
    
    /// Adds loop body definitions to code lines container.
    /// @param lines The code lines container to which loop body definition are added.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_batches_loop_body(      VCodeLines& lines,
                                const bool        diagonal) const;
    
    /// Adds loop end definitions to code lines container.
    /// @param lines The code lines container to which loop end definition are added.
    void _add_batches_loop_end(VCodeLines& lines) const;
    
    /// Adds bra loop start definitions to code lines container.
    /// @param lines The code lines container to which bra loop start definition are added.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_bra_loop_start(      VCodeLines& lines,
                             const bool        diagonal) const;
    
    /// Adds bra loop body definitions to code lines container.
    /// @param lines The code lines container to which bra loop body definition are added.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_bra_loop_body(      VCodeLines&  lines,
                            const I2CIntegral& integral,
                            const bool         diagonal) const;
    
    /// Adds bra loop end definitions to code lines container.
    /// @param lines The code lines container to which bra loop end definition are added.
    void _add_bra_loop_end(VCodeLines& lines) const;
    
    /// Adds bra loop call tree definitions to code lines container.
    /// @param lines The code lines container to which bra loop call tree definition are added.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_loop_call_tree(      VCodeLines&  lines,
                             const I2CIntegral& integral,
                             const bool         diagonal) const;
    
    /// Adds bra loop call tree definitions to code lines container.
    /// @param lines The code lines container to which bra loop call tree definition are added.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_loop_call_tree(      VCodeLines&      lines,
                             const TensorComponent& component,
                             const I2CIntegral&     integral,
                             const bool             bra_first,
                             const bool             diagonal) const;
    
    /// Adds bra loop call tree definitions to code lines container.
    /// @param lines The code lines container to which bra loop start definition are added.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_loop_call_tree(      VCodeLines&      lines,
                             const TensorComponent& bra_component,
                             const TensorComponent& ket_component,
                             const I2CIntegral&     integral,
                             const bool             diagonal) const;
    
    /// Adds primitives loop start definitions to code lines container.
    /// @param lines The code lines container to which primitives loop start definition are added.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_prim_loop_start(      VCodeLines& lines,
                              const bool        diagonal) const;
    
    /// Adds primitives loop end definitions to code lines container.
    /// @param lines The code lines container to which primitives loop end definition are added.
    void _add_prim_loop_end(VCodeLines& lines) const;
    
    /// Adds primitive function call definitions to code lines container.
    /// @param lines The code lines container to which primitives loop start definition are added.
    /// @param spacer The size of formatting shift.
    void _add_prim_call_data(      VCodeLines& lines,
                             const size_t      spacer) const;
    
    /// Adds definition of block distribution call tree for compute function.
    /// @param lines The code lines container to which primitives loop start definition are added.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_block_distributor(      VCodeLines&  lines,
                                  const I2CIntegral& integral,
                                  const bool         diagonal) const;
    
    /// Adds definition of block distribution call tree for compute function.
    /// @param lines The code lines container to which primitives loop start definition are added.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_block_distributor(      VCodeLines&      lines,
                                  const TensorComponent& component,
                                  const I2CIntegral&     integral,
                                  const bool             bra_first,
                                  const bool             diagonal) const;
    
    /// Adds definition  of of block distribution call tree  for compute function.
    /// @param lines The code lines container to which primitives loop start definition are added.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _write_block_distributor(      VCodeLines&      lines,
                                  const TensorComponent& bra_component,
                                  const TensorComponent& ket_component,
                                  const I2CIntegral&     integral,
                                  const bool             diagonal) const;
    
public:
    /// Creates a two-center compute function body generator.
    T2CFuncBodyDriver() = default;
    
    /// Writes body of compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void write_func_body(      std::ofstream& fstream,
                         const I2CIntegral&   integral,
                         const bool           diagonal) const;
};

#endif /* t2c_body_hpp */
