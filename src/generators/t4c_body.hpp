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

#ifndef t4c_body_hpp
#define t4c_body_hpp

#include <string>
#include <vector>
#include <utility>
#include <fstream>

#include "t4c_defs.hpp"
#include "file_stream.hpp"

// Four-center compute function body generators for CPU.
class T4CFuncBodyDriver
{
    /// Generates vector of strings with GTOs definitions in compute function.
    /// @return The vector of strings with GTOS definitions in compute function.
    std::vector<std::string> _get_gto_pairs_def() const;
    
    /// Generates vector of ket factors in compute function.
    /// @param integral The base four center integral.
    /// @return The vector of ket factors in compute function.
    std::vector<std::string> _get_ket_variables_def(const I4CIntegral& integral) const;
    
    /// Generates vector of ket factors in compute function.
    /// @return The vector of ket factors in compute function.
    std::vector<std::string> _get_diag_ket_variables_def() const;
    
    /// Generates vector of distances in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of distances in compute function.
    std::vector<std::string> _get_coordinates_def(const I4CIntegral& integral) const;
    
    /// Generates vector of distances in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of distances in compute function.
    std::vector<std::string> _get_diag_coordinates_def(const I4CIntegral& integral) const;
    
    /// Generates vector of distances in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of distances in compute function.
    std::vector<std::string> _get_full_coordinates_def(const I4CIntegral& integral) const;
    
    /// Generates vector of Cartesian buffer integrals.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @return The vector of Cartesian integrals in compute function.
    SI4CIntegrals _get_cart_buffer_integrals(const SI4CIntegrals& bra_integrals,
                                             const SI4CIntegrals& ket_integrals) const;
    
    /// Generates vector of contracted  buffer integrals.
    /// @param integrals The set of unique integrals.
    /// @return The vector of contracted integrals in compute function.
    SI4CIntegrals _get_contr_buffers_integrals(const SI4CIntegrals& integrals) const;
    
    /// Generates vector of half spherical buffer integrals.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base four center integral.
    /// @return The vector of half spherical integrals in compute function.
    SI4CIntegrals _get_half_spher_buffers_integrals(const SI4CIntegrals& bra_integrals,
                                                    const SI4CIntegrals& ket_integrals,
                                                    const I4CIntegral&   integral) const;
    
    /// Generates vector of primitive buffers in compute function.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_prim_buffers_def(const SI4CIntegrals& integrals,
                                                   const I4CIntegral&   integral) const;
    
    /// Generates vector of primitive buffers in compute function.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_diag_prim_buffers_def(const SI4CIntegrals& integrals,
                                                        const I4CIntegral&   integral) const;
    
    /// Generates vector of primitive buffers in compute function.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_full_prim_buffers_def(const SI4CIntegrals& integrals,
                                                        const I4CIntegral&   integral) const;
    
    /// Generates vector of Cartesian buffers in compute function.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_cart_buffers_def(const SI4CIntegrals& bra_integrals,
                                                   const SI4CIntegrals& ket_integrals,
                                                   const I4CIntegral&   integral) const;
    
    /// Generates vector of Cartesian buffers in compute function.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_diag_cart_buffers_def(const SI4CIntegrals& bra_integrals,
                                                        const SI4CIntegrals& ket_integrals,
                                                        const I4CIntegral&   integral) const;
    
    /// Generates vector of Cartesian buffer integrals.
    /// @param integrals The set of unique integrals.
    /// @param integral The base four center integral.
    /// @return The vector of Cartesian integrals in compute function.
    std::vector<std::string> _get_full_cart_buffers_def(const SI4CIntegrals& integrals,
                                                        const I4CIntegral&   integral) const;
    
    /// Generates vector of contracted buffers in compute function.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_contr_buffers_def(const SI4CIntegrals& bra_integrals,
                                                    const SI4CIntegrals& ket_integrals,
                                                    const I4CIntegral&   integral) const;
    
    /// Generates vector of contracted buffers in compute function.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_diag_contr_buffers_def(const SI4CIntegrals& bra_integrals,
                                                         const SI4CIntegrals& ket_integrals,
                                                         const I4CIntegral&   integral) const;
    
    /// Generates vector of half transformed buffers in compute function.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_half_spher_buffers_def(const SI4CIntegrals& bra_integrals,
                                                         const SI4CIntegrals& ket_integrals,
                                                         const I4CIntegral&   integral) const;
    
    /// Generates vector of half transformed buffers in compute function.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_diag_half_spher_buffers_def(const SI4CIntegrals& bra_integrals,
                                                              const SI4CIntegrals& ket_integrals,
                                                              const I4CIntegral&   integral) const;
    
    /// Generates vector of half transformed buffers in compute function.
    /// @param integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_bra_half_spher_buffers_def(const SI4CIntegrals& integrals,
                                                             const I4CIntegral&   integral) const;
    
    /// Generates vector of half transformed buffers in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_spher_buffers_def(const I4CIntegral& integral) const;
    
    /// Generates vector of half transformed buffers in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_diag_spher_buffers_def(const I4CIntegral& integral) const;
    
    /// Generates vector of half transformed buffers in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_full_spher_buffers_def(const I4CIntegral& integral) const;
    
    /// Generates vector of Boys function definitions in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of Boys function definitions in compute function.
    std::vector<std::string> _get_boys_function_def(const I4CIntegral& integral) const;
    
    /// Generates vector of Boys function definitions in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of Boys function definitions in compute function.
    std::vector<std::string> _get_diag_boys_function_def(const I4CIntegral& integral) const;
    
    /// Generates vector of maximum integral values in compute function.
    /// @return The vector of maximum integral values  in compute function.
    std::vector<std::string> _get_max_array_def() const;
    
    /// Adds loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    void _add_loop_start(      VCodeLines&  lines,
                         const SI4CIntegrals& bra_integrals,
                         const SI4CIntegrals& ket_integrals,
                         const I4CIntegral& integral) const;
    
    /// Adds loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    void _add_diag_loop_start(      VCodeLines&  lines,
                              const SI4CIntegrals& bra_integrals,
                              const SI4CIntegrals& ket_integrals,
                              const I4CIntegral& integral) const;
    
    /// Adds loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of unique integrals.
    /// @param integral The base two center integral.
    void _add_full_loop_start(      VCodeLines&    lines,
                              const SI4CIntegrals& integrals,
                              const I4CIntegral&   integral) const;
    
    /// Adds loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_loop_end(      VCodeLines&  lines,
                       const I4CIntegral& integral) const;
    
    /// Adds loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_diag_loop_end(      VCodeLines&  lines,
                            const I4CIntegral& integral) const;
    
    /// Adds loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_full_loop_end(      VCodeLines&  lines,
                            const I4CIntegral& integral) const;
    
    /// Adds ket loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_ket_loop_start(      VCodeLines&  lines,
                             const I4CIntegral& integral) const;
    
    /// Adds ket loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_diag_ket_loop_start(      VCodeLines&  lines,
                                  const I4CIntegral& integral) const;
    
    /// Adds ket loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_full_ket_loop_start(      VCodeLines&  lines,
                                  const I4CIntegral& integral) const;
    
    /// Adds ket loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param vrr_integrals The set of unique integrals for vertical recursion.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    void _add_ket_loop_end(      VCodeLines&    lines,
                           const SI4CIntegrals& vrr_integrals,
                           const SI4CIntegrals& bra_integrals,
                           const SI4CIntegrals& ket_integrals,
                           const I4CIntegral&   integral) const;
    
    /// Adds ket loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    void _add_diag_ket_loop_end(      VCodeLines&    lines,
                                const SI4CIntegrals& bra_integrals,
                                const SI4CIntegrals& ket_integrals,
                                const I4CIntegral&   integral) const;
    
    /// Adds ket loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_full_ket_loop_end(      VCodeLines&    lines,
                                const I4CIntegral&   integral) const;
    
    /// Adds auxilary integrals.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    void _add_auxilary_integrals(      VCodeLines&    lines,
                                 const SI4CIntegrals& integrals,
                                 const I4CIntegral&   integral) const;
        
    /// Adds call tree for vertical recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    void _add_vrr_call_tree(      VCodeLines&  lines,
                            const SI4CIntegrals& integrals,
                            const I4CIntegral&   integral) const;
    
    /// Adds call tree for vertical recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    void _add_full_vrr_call_tree(      VCodeLines&  lines,
                                 const SI4CIntegrals& integrals) const;
    
    /// Gets arguments list for primitive vertical recursion function call.
    /// @param start The indexes starting position.
    /// @param integrals The set of inetrgals.
    /// @param integral The base four center integral.
    std::string _get_vrr_arguments(const size_t start,
                                   const SI4CIntegrals& integrals,
                                   const I4CIntegral&   integral) const;
    
    /// Gets arguments list for primitive vertical recursion function call.
    /// @param integral The base four center integral.
    std::string _get_full_vrr_arguments(const I4CIntegral& integral) const;
    
    /// Adds call tree for vertical recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    /// @param integral The base four center integral.
    void _add_geom_call_tree(      VCodeLines&    lines,
                             const SI4CIntegrals& integrals,
                             const I4CIntegral&   integral) const;
    
    /// Gets arguments list for primitive vertical recursion function call.
    /// @param integrals The set of inetrgals.
    /// @param integral The base four center integral.
    std::string _get_geom_arguments(const SI4CIntegrals& integrals,
                                    const I4CIntegral&   integral) const;
    
    /// Adds call tree for ket horizontal recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    void _add_ket_hrr_call_tree(      VCodeLines&  lines,
                                const SI4CIntegrals& bra_integrals,
                                const SI4CIntegrals& ket_integrals) const;
    
    /// Gets arguments list for ket horizontal recursion function call.
    /// @param start The starting index of arguments list. 
    /// @param integral The base four center integral.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    std::string _get_ket_hrr_arguments(const size_t       start,
                                       const I4CIntegral& integral,
                                       const SI4CIntegrals& bra_integrals,
                                       const SI4CIntegrals& ket_integrals) const;
    
    /// Adds call tree for ket side transformation.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    void _add_ket_trafo_call_tree(      VCodeLines&  lines,
                                  const SI4CIntegrals& bra_integrals,
                                  const SI4CIntegrals& ket_integrals,
                                  const I4CIntegral&   integral) const;
    
    /// Adds call tree for bra horizontal recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    void _add_bra_hrr_call_tree(      VCodeLines&  lines,
                                const SI4CIntegrals& bra_integrals,
                                const SI4CIntegrals& ket_integrals,
                                const I4CIntegral&   integral) const;
    
    /// Gets arguments list for bra horizontal recursion function call.
    /// @param integral The base four center integral.
    /// @param integrals The set of inetrgals.
    std::string _get_bra_hrr_arguments(const size_t start,
                                       const I4CIntegral& integral,
                                       const SI4CIntegrals& integrals) const;
    
    /// Adds call tree for bra side transformation.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    void _add_bra_trafo_call_tree(      VCodeLines&  lines,
                                  const SI4CIntegrals& bra_integrals,
                                  const SI4CIntegrals& ket_integrals,
                                  const I4CIntegral&   integral) const;
    
    /// Adds call for full transformation.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_full_trafo(      VCodeLines&  lines,
                         const I4CIntegral& integral) const;
    
    /// Checks if coordinates of center W are required for integration.
    /// @param integral The base four center integral.
    bool _need_center_w(const I4CIntegral& integral) const;
    
    /// Checks if distances of (Q-D) are required for integration.
    /// @param integral The base four center integral.
    bool _need_distances_qd(const I4CIntegral& integral) const;
    
    /// Checks if distances of (W-Q) are required for integration.
    /// @param integral The base four center integral.
    bool _need_distances_wq(const I4CIntegral& integral) const;
    
    /// Checks if distances of (W-P) are required for integration.
    /// @param integral The base four center integral.
    bool _need_distances_wp(const I4CIntegral& integral) const;
    
    /// Checks if horizontal recursion on ket side is required for integration.
    /// @param integral The base four center integral.
    bool _need_hrr_for_ket(const I4CIntegral& integral) const;
    
    /// Checks if horizontal recursion on bra side is required for integration.
    /// @param integral The base four center integral.
    bool _need_hrr_for_bra(const I4CIntegral& integral) const;
    
    /// Gets index of Cartesian center W in factors buffer.
    /// @param integral The base four center integral.
    size_t _get_index_w(const I4CIntegral& integral) const;
    
    /// Gets index of distances of (Q-D)  in factors buffer.
    /// @param integral The base four center integral.
    size_t _get_index_qd(const I4CIntegral& integral) const;
    
    /// Gets index of distances of (W-Q)  in factors buffer.
    /// @param integral The base four center integral.
    size_t _get_index_wq(const I4CIntegral& integral) const;
    
    /// Gets index of distances of (W-P)  in factors buffer.
    /// @param integral The base four center integral.
    size_t _get_index_wp(const I4CIntegral& integral) const;
    
    /// Gets index of requested integral in set of integrals.
    /// @param start The initial index.
    /// @param integral The base four center integral.
    /// @param integrals The set of inetrgals.
    size_t _get_index(const size_t         start,
                      const I4CIntegral&   integral,
                      const SI4CIntegrals& integrals) const;
    
    /// Gets index of requested integral in set of half transformed integrals.
    /// @param start The initial index.
    /// @param integral The base four center integral.
    /// @param integrals The set of inetrgals.
    size_t _get_half_spher_index(const size_t         start,
                                 const I4CIntegral&   integral,
                                 const SI4CIntegrals& integrals) const;
    
    /// Gets total number of components in set of integrals.
    /// @param integrals The set of inetrgals.
    size_t _get_all_components(const SI4CIntegrals& integrals) const;
    
    /// Gets total number of half transformed components in set of integrals.
    /// @param integrals The set of inetrgals.
    size_t _get_all_half_spher_components(const SI4CIntegrals& integrals) const;
    
    /// Gets total number of spherical components in set of integrals.
    /// @param integral The base four center integral.
    size_t _get_all_spher_components(const I4CIntegral& integral) const;

public:
    /// Creates a four-center compute function body generator.
    T4CFuncBodyDriver() = default;
    
    /// Writes body of compute function.
    /// @param fstream the file stream.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param vrr_integrals The set of unique integrals for vertical recursion.
    /// @param integral The base two center integral.
    void write_func_body(      std::ofstream& fstream,
                         const SI4CIntegrals& bra_integrals,
                         const SI4CIntegrals& ket_integrals,
                         const SI4CIntegrals& vrr_integrals,
                         const I4CIntegral&   integral) const;
    
    /// Writes body of compute function.
    /// @param fstream the file stream.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param vrr_integrals The set of unique integrals for vertical recursion.
    /// @param integral The base two center integral.
    void write_diag_func_body(      std::ofstream& fstream,
                              const SI4CIntegrals& bra_integrals,
                              const SI4CIntegrals& ket_integrals,
                              const SI4CIntegrals& vrr_integrals,
                              const I4CIntegral&   integral) const;
    
    /// Writes body of compute function.
    /// @param fstream the file stream.
    /// @param geom_integrals The set of unique integrals for geometrical recursion.
    /// @param vrr_integrals The set of unique integrals for vertical recursion.
    /// @param integral The base two center integral.
    void write_geom_func_body(      std::ofstream& fstream,
                              const SI4CIntegrals& geom_integrals,
                              const SI4CIntegrals& vrr_integrals,
                              const I4CIntegral&   integral) const;
};


#endif /* t4c_body_hpp */
