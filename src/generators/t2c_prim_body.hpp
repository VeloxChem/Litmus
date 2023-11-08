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
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _add_func_pragma(      VCodeLines&    lines,
                          const I2CIntegral&   integral,
                          const bool           sum_form) const;
    
    /// Adds pragmas for primitive compute function.
    /// @param lines The code lines container to which  pragma definitions are added.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @param bra_first The flag to set bra as expansion point.
    void _add_func_pragma(      VCodeLines&      lines,
                          const TensorComponent& component,
                          const I2CIntegral&     integral,
                          const bool             sum_form,
                          const bool             bra_first) const;
    
    /// Adds pragmas for primitive compute function.
    /// @param lines The code lines container to which  pragma definitions are added.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _add_func_pragma(      VCodeLines&      lines,
                          const TensorComponent& bra_component,
                          const TensorComponent& ket_component,
                          const I2CIntegral&     integral,
                          const bool             sum_form) const;
    
    /// Adds pragmas for primitive compute function.
    /// @param lines The code lines container to which common pragma definitions are added.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _add_common_pragma(      VCodeLines& lines,
                            const bool        sum_form) const;
    
    /// Adds main loop start for primitive compute function.
    /// @param lines The code lines container to which loop start are added.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _add_loop_start(      VCodeLines&  lines,
                         const I2CIntegral& integral,
                         const bool         sum_form) const;
    
    /// Adds index loop start for primitive compute function.
    /// @param lines The code lines container to which loop start are added.
    /// @param integral The base two center integral.
    void _add_index_loop_start(      VCodeLines&  lines,
                               const I2CIntegral& integral) const;
    
    /// Adds overlap integral variables after the main loop start for primitive compute function.
    /// @param lines The code lines container to which loop start are added.
    /// @param integral The base two center integral.
    void _add_overlap_vars(      VCodeLines&  lines,
                           const I2CIntegral& integral) const;
    
    /// Adds kinetic energy integral variables after the main loop start for primitive compute function.
    /// @param lines The code lines container to which loop start are added.
    /// @param integral The base two center integral.
    void _add_kinetic_energy_vars(      VCodeLines&  lines,
                                  const I2CIntegral& integral) const;
    
    /// Adds nuclear potential integral variables after the main loop start for primitive compute function.
    /// @param lines The code lines container to which loop start are added.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _add_nuclear_potential_vars(      VCodeLines&  lines,
                                     const I2CIntegral& integral,
                                     const bool         sum_form) const;
    
    /// Adds nuclear potential geometrical derivatives integral variables after the main loop start for primitive compute function.
    /// @param lines The code lines container to which loop start are added.
    /// @param integral The base two center integral.
    void _add_nuclear_potential_geom_vars(      VCodeLines&  lines,
                                          const I2CIntegral& integral) const;
    
    /// Adds multipole integral variables after the main loop start for primitive compute function.
    /// @param lines The code lines container to which loop start are added.
    /// @param integral The base two center integral.
    void _add_multipole_vars(      VCodeLines&  lines,
                             const I2CIntegral& integral) const;
    
    /// Adds three center overlap integral variables after the main loop start for primitive compute function.
    /// @param lines The code lines container to which loop start are added.
    /// @param integral The base two center integral.
    void _add_three_center_overlap_vars(      VCodeLines&  lines,
                                        const I2CIntegral& integral) const;
    
    /// Adds main loop end for primitive compute function.
    /// @param lines The code lines container to which loop end are added.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _add_loop_end(      VCodeLines&  lines,
                       const bool         sum_form) const;
    
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
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _add_simd_code(      VCodeLines&               lines,
                        const std::vector<std::string>& labels,
                        const VT2CIntegrals&            components,
                        const I2CIntegral&              integral,
                        const bool                      sum_form) const;
    
    /// Generates  the recursion  group for given vector of integral components.
    /// @param components The vector of integral component.
    /// @param integral The base integral.
    /// @return The recursion group.
    R2Group _generate_integral_group(const VT2CIntegrals& components,
                                     const I2CIntegral&   integral) const;
    
    /// Adds conditional prefactors for given recursion group.
    /// @param lines The code lines container to which simd code are added.
    /// @param rgroup The recursion group.
    /// @param integral The base integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _add_prefactors(      VCodeLines&  lines,
                         const R2Group&     rgroup,
                         const I2CIntegral& integral,
                         const bool         sum_form) const;
    
    /// Adds simd code lines for given reccursion group.
    /// @param lines The code lines container to which simd code are added.
    /// @param labels The vector of labels for integral components.
    /// @param rgroup The recursion group.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _add_simd_lines(      VCodeLines&               lines,
                         const std::vector<std::string>& labels,
                         const R2Group&                  rgroup,
                         const bool                      sum_form) const;
    
    /// Adds block of  simd code lines for given reccursion distribution.
    /// @param lines The code lines container to which simd code are added.
    /// @param label The label of integral components.
    /// @param integral The integral  component.
    /// @param rdist  The recursion distribution.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _add_simd_lines_block(      VCodeLines&  lines,
                               const std::string& label,
                               const T2CIntegral& integral,
                               const R2CDist&     rdist,
                               const bool         sum_form) const;
    
    /// Generates standard label  of auxilary integral.
    /// @param integral The integral  component.
    /// @param base The base integral.
    /// @return The label of auxilary integral.
    std::string _get_aux_label(const T2CIntegral& integral,
                               const T2CIntegral& base) const;
    
    /// Generates vector of special variable strings.
    /// @param integral The base two center integral.
    /// @param geom_form The geometrical form of special parameters.
    /// @return The vector of special variable strings.
    std::vector<std::string> _get_special_vars_str(const I2CIntegral& integral,
                                                   const bool         geom_form) const;
    
    /// Generates vector of Boys function variables strings.
    /// @param integral The base two center integral.
    /// @return The vector of special variable strings.
    std::vector<std::string> _get_boys_vars_str(const I2CIntegral& integral) const;
    
    /// Adds Boys function computation code lines.
    /// @param lines The code lines container to which simd code are added.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _add_boys_compute_lines(      VCodeLines&  lines,
                                 const I2CIntegral& integral,
                                 const bool         sum_form) const;

public:
    /// Creates a two-center compute function body generator.
    T2CPrimFuncBodyDriver() = default;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void write_prim_func_body(      std::ofstream& fstream,
                              const I2CIntegral&   integral,
                              const bool           sum_form) const;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @param bra_first The flag to set bra as expansion point.
    void write_prim_func_body(      std::ofstream&   fstream,
                              const TensorComponent& component,
                              const I2CIntegral&     integral,
                              const bool             sum_form,
                              const bool             bra_first) const;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param bra_component the integral component on bra side.
    /// @param ket_component the integral component on ket side.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void write_prim_func_body(      std::ofstream&   fstream,
                              const TensorComponent& bra_component,
                              const TensorComponent& ket_component,
                              const I2CIntegral&     integral,
                              const bool             sum_form) const;
    
    
};

#endif /* t2c_prim_body_hpp */
