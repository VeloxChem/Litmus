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
    /// Generates vector of strings with external data definitions in compute function.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @return The vector of strings with external data definitions in compute function.
    std::vector<std::string> _get_external_data_def(const I2CIntegral&           integral,
                                                    const std::pair<bool, bool>& rec_form) const;
    
    /// Generates vector of strings with GTOs definitions in compute function.
    /// @return The vector of strings with GTOS definitions in compute function.
    std::vector<std::string> _get_gtos_def() const;
    
    /// Generates vector of ket factors in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of ket factors in compute function.
    std::vector<std::string> _get_ket_variables_def(const I2CIntegral& integral) const;
    
    /// Generates vector of buffers in compute function.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    /// @param geom_drvs The geometrical derivative of bra side, integrand, and  ket side.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_buffers_def(const SI2CIntegrals& integrals,
                                              const I2CIntegral&   integral,
                                              const std::array<int, 3>& geom_drvs) const;
    
    /// Generates vector of Boys function definitions in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of Boys function definitions in compute function.
    std::vector<std::string> _get_boys_function_def(const I2CIntegral& integral) const;
    
    /// Adds loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_loop_start(      VCodeLines&  lines,
                         const I2CIntegral& integral) const;
    
    /// Adds loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    void _add_loop_end(      VCodeLines&  lines,
                       const I2CIntegral& integral,
                       const std::pair<bool, bool>& rec_form) const;
    
    /// Adds ket loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    void _add_ket_loop_start(      VCodeLines&            lines,
                             const I2CIntegral&           integral,
                             const std::pair<bool, bool>& rec_form) const;
    
    /// Adds ket loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    void _add_ket_loop_end(      VCodeLines&  lines,
                           const SI2CIntegrals& integrals,
                           const I2CIntegral& integral,
                           const std::pair<bool, bool>& rec_form) const;
    
    /// Adds sum loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @param use_rs The flag for use of range-separated Coulomb interactions.
    void _add_sum_loop_start(      VCodeLines&            lines,
                             const I2CIntegral&           integral,
                             const std::pair<bool, bool>& rec_form,
                             const bool                   use_rs) const;
    
    /// Adds sum loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    void _add_sum_loop_end(      VCodeLines&            lines,
                           const SI2CIntegrals&         integrals,
                           const I2CIntegral&           integral,
                           const std::pair<bool, bool>& rec_form) const;
    
    /// Adds auxilary integrals.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @param in_sum_loop The flag indicating call from inside summation loop.
    void _add_auxilary_integrals(      VCodeLines&            lines,
                                 const SI2CIntegrals&         integrals,
                                 const I2CIntegral&           integral,
                                 const std::pair<bool, bool>& rec_form,
                                 const bool                   in_sum_loop) const;
    
    /// Adds call tree for recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    void _add_call_tree(      VCodeLines&            lines,
                        const SI2CIntegrals&         integrals,
                        const I2CIntegral&           integral,
                        const std::pair<bool, bool>& rec_form) const;
    
    
    /// Adds call tree for recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param geom_integrals The set of inetrgals in geometrical recursion.
    /// @param vrr_integrals The set of inetrgals in vertical recursion.
    ///
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    void _add_geom_call_tree(      VCodeLines&            lines,
                             const SI2CIntegrals&         geom_integrals,
                             const SI2CIntegrals&         vrr_integrals,
                             const I2CIntegral&           integral,
                             const std::array<int, 3>&    geom_drvs,
                             const std::pair<bool, bool>& rec_form) const;
    
    /// Gets arguments list for primitive function call.
    /// @param integral The base two center integral.
    std::string _get_arguments(const I2CIntegral& integral) const;
    
    /// Gets arguments list for primitive function call.
    /// @param integral The base two center integral.
    /// @param integrals The set of inetrgals.
    std::string _get_arguments(const I2CIntegral&   integral,
                               const SI2CIntegrals& integrals) const;
    
    /// Gets position of integral in integrals buffer.
    /// @param integral The base two center integral.
    /// @param integrals The set of inetrgals.
    size_t _get_position(const I2CIntegral&   integral,
                         const SI2CIntegrals& integrals) const;
    
    /// Checks if coordinates of center P are required for integration.
    /// @param integral The base two center integral.
    bool _need_center_p(const I2CIntegral& integral) const;
    
    /// Checks if distances of (P-C) are required for integration.
    /// @param integral The base two center integral.
    bool _need_distances_pc_in_call_tree(const I2CIntegral& integral) const;
    
    /// Checks if distances of (P-C) are required for integration.
    /// @param integral The base two center integral.
    bool _need_distances_pc(const I2CIntegral& integral) const;
    
    /// Checks if distances of (P-A) are required for integration.
    /// @param integral The base two center integral.
    bool _need_distances_pa(const I2CIntegral& integral) const;
    
    /// Checks if distances of (P-B) are required for integration.
    /// @param integral The base two center integral.
    bool _need_distances_pb(const I2CIntegral& integral) const;
    
    /// Checks if exponents on center A  are required for integration.
    /// @param integral The base two center integral.
    bool _need_exponents(const I2CIntegral& integral) const;
    
    /// Checks if Boys function data are required for integration.
    /// @param integral The base two center integral.
    bool _need_boys_func(const I2CIntegral& integral) const;
    
    /// Checks if external coordinates are needed.
    /// @param integral The base two center integral.
    bool _need_external_coords(const I2CIntegral& integral) const;
    
    /// Checks if geometrical derivatives are needed.
    /// @param geom_drvs The geometrical derivative of bra side, integrand, and  ket side.
    bool _need_geom_drvs(const std::array<int, 3>& geom_drvs) const;
    
    /// Gets index of R(PA) distances in factors list.
    /// @param integral The base two center integral.
    /// @return The string with index label.
    int _get_index_pa(const I2CIntegral& integral) const;
    
    /// Gets index of R(PB) distances in factors list.
    /// @param integral The base two center integral.
    /// @return The string with index label.
    int _get_index_pb(const I2CIntegral& integral) const;
    
    /// Gets index of R(PC) distances in factors list.
    /// @param integral The base two center integral.
    /// @return The string with index label.
    int _get_index_pc(const I2CIntegral& integral) const;
    
public:
    /// Creates a two-center compute function body generator.
    T2CFuncBodyDriver() = default;
    
    /// Writes body of compute function.
    /// @param fstream the file stream.
    /// @param geom_integrals The set of inetrgals in geometrical recursion.
    /// @param vrr_integrals The set of inetrgals in vertical recursion.
    /// @param integral The base two center integral.
    /// @param geom_drvs The geometrical derivative of bra side, integrand, and  ket side.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @param use_rs The flag for use of range-separated Coulomb interactions.
    void write_func_body(      std::ofstream&         fstream,
                         const SI2CIntegrals&         geom_integrals,
                         const SI2CIntegrals&         vrr_integrals,
                         const I2CIntegral&           integral,
                         const std::array<int, 3>& geom_drvs, 
                         const std::pair<bool, bool>& rec_form,
                         const bool                   use_rs) const;
};

#endif /* t2c_body_hpp */
