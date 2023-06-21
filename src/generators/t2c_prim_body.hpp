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
#include <vector>
#include <fstream>

#include "t2c_defs.hpp"
#include "file_stream.hpp"

// Two-center primitive compute function body generators for CPU.
class T2CPrimFuncBodyDriver
{
    /// Generates the common data definitions in primitive compute function.
    /// @return The vector of common primitives  data.
    std::vector<std::string> _get_common_data_str() const;
    
    /// Generates the buffer definitions in primitive compute function.
    /// @return The vector of common primitive buffers.
    std::vector<std::string> _get_buffers_str(const I2CIntegral& integral) const;
    
    /// Generates the buffer definitions in primitive compute function.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    /// @return The vector of common primitive buffers.
    std::vector<std::string> _get_buffers_str(const TensorComponent& component,
                                              const I2CIntegral&     integral,
                                              const bool             bra_first) const;
    
    /// Generates the buffer definitions in primitive compute function.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    /// @return The vector of common primitive buffers.
    std::vector<std::string> _get_buffers_str(const TensorComponent& bra_component,
                                              const TensorComponent& ket_component,
                                              const I2CIntegral&     integral) const;
    
    
    /// Adds pragmas for primitive compute function.
    /// @param lines The code lines container to which  pragma definitions are added.
    /// @param integral The base two center integral.
    void _add_func_pragma(      VCodeLines&    lines,
                          const I2CIntegral&   integral) const;
    
    /// Adds pragmas for primitive compute function.
    /// @param lines The code lines container to which  pragma definitions are added.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    void _add_func_pragma(      VCodeLines&      lines,
                          const TensorComponent& component,
                          const I2CIntegral&     integral,
                          const bool             bra_first) const;
    
    /// Adds pragmas for primitive compute function.
    /// @param lines The code lines container to which  pragma definitions are added.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    void _add_func_pragma(      VCodeLines&      lines,
                          const TensorComponent& bra_component,
                          const TensorComponent& ket_component,
                          const I2CIntegral&     integral) const;
    
    /// Adds pragmas for primitive compute function.
    /// @param lines The code lines container to which common pragma definitions are added.
    void _add_common_pragma(VCodeLines& lines) const;
    
    /// Adds main loop start for primitive compute function.
    /// @param lines The code lines container to which loop start are added.
    /// @param integral The base two center integral.
    void _add_loop_start(      VCodeLines&  lines,
                         const I2CIntegral& integral) const;
    
    /// Adds main loop end for primitive compute function.
    /// @param lines The code lines container to which loop end are added.
    void _add_loop_end(VCodeLines&  lines) const;
    
    /// Select integrals components with predefined bra or ket side components.
    /// @param component the tensor component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    VT2CIntegrals _select_integral_components(const TensorComponent& component,
                                              const I2CIntegral&     integral,
                                              const bool             bra_first) const;
    /// Writes documentation string for primitive compute function.
    /// @param bra_component the tensor component.
    /// @param ket_component the tensor component.
    /// @param integral The base two center integral.
    VT2CIntegrals _select_integral_components(const TensorComponent& bra_component,
                                              const TensorComponent& ket_component,
                                              const I2CIntegral&     integral) const;
    
    /// Adds simd code generated for selected set of integral components.
    /// @param lines The code lines container to which simd code are added.
    /// @param labels The vector of labels for integral components.
    /// @param components The vector of integral component.
    /// @param integral The base integral.
    void _add_simd_code(      VCodeLines&               lines,
                        const std::vector<std::string>& labels,
                        const VT2CIntegrals&            components,
                        const I2CIntegral&              integral) const;
    
    /// Generates  the recursion  group for given vector of integral components.
    /// @param components The vector of integral component.
    /// @param integral The base integral.
    /// @return The recursion group.
    R2Group _generate_integral_group(const VT2CIntegrals& components,
                                     const I2CIntegral&   integral) const;
    
    /// Adds conditional prefactors for given recursion group.
    /// @param lines The code lines container to which simd code are added.
    /// @param rgroup The recursion group.
    void _add_prefactors(      VCodeLines& lines,
                         const R2Group&    rgroup) const;
    
    /// Adds simd code lines for given reccursion group.
    /// @param lines The code lines container to which simd code are added.
    /// @param labels The vector of labels for integral components.
    /// @param rgroup The recursion group.
    void _add_simd_lines(      VCodeLines&               lines,
                         const std::vector<std::string>& labels,
                         const R2Group&                  rgroup) const;

public:
    /// Creates a two-center compute function body generator.
    T2CPrimFuncBodyDriver() = default;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void write_prim_func_body(      std::ofstream& fstream,
                              const I2CIntegral&   integral) const;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param bra_first The flag to set bra as expansion point.
    void write_prim_func_body(      std::ofstream&   fstream,
                              const TensorComponent& component,
                              const I2CIntegral&     integral,
                              const bool             bra_first) const;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    void write_prim_func_body(      std::ofstream&   fstream,
                              const TensorComponent& bra_component,
                              const TensorComponent& ket_component,
                              const I2CIntegral&     integral) const;
};

#endif /* t2c_prim_body_hpp */
