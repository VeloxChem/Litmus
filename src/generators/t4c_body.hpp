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
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of strings with GTOS definitions in compute function.
    std::vector<std::string> _get_gto_pairs_def(const bool diagonal) const;
    
    /// Generates vector of ket factors in compute function.
    /// @return The vector of ket factors in compute function.
    std::vector<std::string> _get_ket_variables_def(const bool diagonal) const;
    
    /// Generates vector of distances in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of distances in compute function.
    std::vector<std::string> _get_coordinates_def(const I4CIntegral& integral) const;
    
    /// Generates vector of Cartesian buffer integrals.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @return The vector of Cartesian integrals in compute function.
    SI4CIntegrals _get_cart_buffer_integrals(const SI4CIntegrals& bra_integrals,
                                             const SI4CIntegrals& ket_integrals) const;
    
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
    
    /// Generates vector of contracted buffers in compute function.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_contr_buffers_def(const SI4CIntegrals& bra_integrals,
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
    /// @param integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_bra_half_spher_buffers_def(const SI4CIntegrals& integrals,
                                                             const I4CIntegral&   integral) const;
    
    /// Generates vector of half transformed buffers in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_spher_buffers_def(const I4CIntegral& integral) const;
    
    /// Generates vector of Boys function definitions in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of Boys function definitions in compute function.
    std::vector<std::string> _get_boys_function_def(const I4CIntegral& integral) const;
    
    /// Adds loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_loop_start(      VCodeLines&  lines,
                         const SI4CIntegrals& bra_integrals,
                         const SI4CIntegrals& ket_integrals,
                         const I4CIntegral& integral,
                         const bool         diagonal) const;
    
    /// Adds loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_loop_end(      VCodeLines&  lines,
                       const I4CIntegral& integral,
                       const bool         diagonal) const;
    
    /// Adds ket loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_ket_loop_start(      VCodeLines&  lines,
                             const I4CIntegral& integral,
                             const bool         diagonal) const;
    
    /// Adds ket loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_ket_loop_end(      VCodeLines&    lines,
                           const SI4CIntegrals& bra_integrals,
                           const SI4CIntegrals& ket_integrals,
                           const I4CIntegral&   integral,
                           const bool           diagonal) const;
    
    /// Adds auxilary integrals.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    void _add_auxilary_integrals(      VCodeLines&    lines,
                                 const SI4CIntegrals& integrals) const;
    
    /// Adds call tree for vertical recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    void _add_vrr_call_tree(      VCodeLines&  lines,
                            const SI4CIntegrals& integrals) const;
    
    /// Gets arguments list for primitive vertical recursion function call.
    /// @param integral The base four center integral.
    std::string _get_vrr_arguments(const I4CIntegral& integral) const;
    
    /// Adds call tree for ket horizontal recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    void _add_ket_hrr_call_tree(      VCodeLines&  lines,
                                const SI4CIntegrals& integrals) const;
    
    /// Gets arguments list for ket horizontal recursion function call.
    /// @param integral The base four center integral.
    std::string _get_ket_hrr_arguments(const I4CIntegral& integral) const;
    
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
    /// @param integrals The set of inetrgals.
    void _add_bra_hrr_call_tree(      VCodeLines&  lines,
                                const SI4CIntegrals& integrals) const;
    
    /// Gets arguments list for bra horizontal recursion function call.
    /// @param integral The base four center integral.
    std::string _get_bra_hrr_arguments(const I4CIntegral& integral) const;
    
    /// Adds call tree for bra side transformation.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_bra_trafo_call_tree(      VCodeLines&  lines,
                                  const I4CIntegral& integral) const;
    
public:
    /// Creates a four-center compute function body generator.
    T4CFuncBodyDriver() = default;
    
    /// Writes body of compute function.
    /// @param fstream the file stream.
    /// @param bra_integrals The set of unique integrals for bra horizontal recursion.
    /// @param ket_integrals The set of unique integrals for ket horizontal recursion.
    /// @param vrr_integrals The set of unique integrals for vertical recursion.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void write_func_body(      std::ofstream& fstream,
                         const SI4CIntegrals& bra_integrals,
                         const SI4CIntegrals& ket_integrals,
                         const SI4CIntegrals& vrr_integrals,
                         const I4CIntegral&   integral,
                         const bool           diagonal) const;
    
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
