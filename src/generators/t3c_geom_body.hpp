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

#ifndef t3c_geom_body_hpp
#define t3c_geom_body_hpp

#include <string>
#include <array>
#include <vector>
#include <utility>
#include <fstream>

#include "t3c_defs.hpp"
#include "file_stream.hpp"

// Four-center compute function body generators for CPU.
class T3CGeomFuncBodyDriver
{
    /// Generates vector of strings with GTOs definitions in compute function.
    /// @return The vector of strings with GTOS definitions in compute function.
    std::vector<std::string> _get_gto_pairs_def() const;
    
    /// Generates vector of ket factors in compute function.
    /// @param integral The base four center integral.
    /// @return The vector of ket factors in compute function.
    std::vector<std::string> _get_ket_variables_def(const I3CIntegral& integral) const;
    
    /// Checks if coordinates of center W are required for integration.
    /// @param integral The base four center integral.
    bool _need_center_w(const I3CIntegral& integral) const;
    
    /// Checks if distances of (Q-D) are required for integration.
    /// @param integral The base four center integral.
    bool _need_distances_qd(const I3CIntegral& integral) const;
    
    /// Checks if distances of (W-Q) are required for integration.
    /// @param integral The base four center integral.
    bool _need_distances_wq(const I3CIntegral& integral) const;
    
    /// Checks if distances of (W-A) are required for integration.
    /// @param integral The base four center integral.
    bool _need_distances_wa(const I3CIntegral& integral) const;
    
    /// Checks if horizontal recursion is required for integration.
    /// @param integral The base four center integral.
    bool _need_hrr(const I3CIntegral& integral) const;
    
    /// Gets total number of components in set of integrals.
    /// @param integrals The set of inetrgals.
    size_t _get_all_components(const SI3CIntegrals& integrals) const;
    
    /// Generates vector of primitive buffers in compute function.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_prim_buffers_def(const SI3CIntegrals& integrals,
                                                   const I3CIntegral&   integral) const;
    
    /// Generates vector of Cartesian buffers in compute function.
    /// @param cterms The set of filtered geometrical terms.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_cart_buffers_def(const SG3Terms&    cterms,
                                                   const I3CIntegral& integral) const;
    
    /// Generates vector of half transformed buffers in compute function.
    /// @param skterms The set of filtered geometrical terms.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_half_spher_buffers_def(const SG3Terms&    skterms,
                                                         const I3CIntegral& integral) const;
    
    /// Generates vector of half transformed buffers in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_spher_buffers_def(const I3CIntegral& integral) const;
    
    /// Generates vector of Boys function definitions in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of Boys function definitions in compute function.
    std::vector<std::string> _get_boys_function_def(const I3CIntegral& integral) const;
    
    /// Adds loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_loop_start(      VCodeLines&  lines,
                         const I3CIntegral& integral) const;
    
    /// Adds loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_loop_end(      VCodeLines&  lines,
                       const I3CIntegral& integral) const;
    
    /// Adds ket loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_ket_loop_start(      VCodeLines&  lines,
                             const I3CIntegral& integral) const;
    
    /// Adds ket loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param cterms The set of filtered geometrical terms.
    /// @param vrr_integrals The set of unique integrals for vertical recursion.
    /// @param integral The base two center integral.
    void _add_ket_loop_end(      VCodeLines&    lines,
                           const SG3Terms&      cterms,
                           const SI3CIntegrals& vrr_integrals,
                           const I3CIntegral&   integral) const;
    
    /// Gets index of Cartesian center W in factors buffer.
    /// @param integral The base four center integral.
    size_t _get_index_w(const I3CIntegral& integral) const;
    
    /// Gets index of distances of (Q-D)  in factors buffer.
    /// @param integral The base four center integral.
    size_t _get_index_qd(const I3CIntegral& integral) const;
    
    /// Gets index of distances of (W-Q)  in factors buffer.
    /// @param integral The base four center integral.
    size_t _get_index_wq(const I3CIntegral& integral) const;
    
    /// Gets index of distances of (W-P)  in factors buffer.
    /// @param integral The base four center integral.
    size_t _get_index_wa(const I3CIntegral& integral) const;
    
    /// Adds auxilary integrals.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    /// @param spacer The tabulation spacer.
    void _add_auxilary_integrals(      VCodeLines&    lines,
                                 const SI3CIntegrals& integrals,
                                 const I3CIntegral&   integral,
                                 const size_t         spacer) const;
    
    /// Gets index of requested integral in set of integrals.
    /// @param start The initial index.
    /// @param integral The base four center integral.
    /// @param integrals The set of inetrgals.
    size_t _get_index(const size_t         start,
                      const I3CIntegral&   integral,
                      const SI3CIntegrals& integrals) const;
    
    /// Adds call tree for vertical recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    /// @param spacer The tabulation spacer.
    void _add_vrr_call_tree(      VCodeLines&    lines,
                            const SI3CIntegrals& integrals,
                            const I3CIntegral&   integral,
                            const size_t         spacer) const;
    
    /// Adds call tree for ket horizontal recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param cterms The set of filtered geometrical terms.
    /// @param skterms The set of filtered geometrical terms.
    /// @param integral The base two center integral.
    /// @param spacer The tabulation spacer.
    void _add_ket_hrr_call_tree(      VCodeLines&  lines,
                                const SG3Terms&    cterms,
                                const SG3Terms&    skterms,
                                const I3CIntegral& integral,
                                const size_t       spacer) const;
    
    /// Gets arguments list for primitive vertical recursion function call.
    /// @param start The indexes starting position.
    /// @param integrals The set of inetrgals.
    /// @param integral The base four center integral.
    std::string _get_vrr_arguments(const size_t start,
                                   const SI3CIntegrals& integrals,
                                   const I3CIntegral&   integral) const;
    
    /// Gets index of requested integral in set of integrals.
    /// @param term The base four center term.
    /// @param terms The set of four center terms.
    size_t _get_index(const G3Term&   term,
                      const SG3Terms& terms) const;
    
    /// Adds call tree for ket side transformation.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param skterms The set of filtered geometrical terms.
    /// @param integral The base two center integral.
    void _add_ket_trafo_call_tree(      VCodeLines&  lines,
                                  const SG3Terms&    skterms,
                                  const I3CIntegral& integral) const;
    
    /// Gets index of requested integral in set of half transformed integrals.
    /// @param term The base four center term.
    /// @param terms The set of four center terms.
    size_t _get_half_spher_index(const G3Term&   term,
                                 const SG3Terms& terms) const;
    
    /// Adds call tree for bra side transformation.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param cterms The set of filtered geometrical terms.
    /// @param skterms The set of filtered geometrical terms.
    /// @param integral The base two center integral.
    void _add_bra_trafo_call_tree(      VCodeLines&  lines,
                                  const SG3Terms&    cterms,
                                  const SG3Terms&    skterms,
                                  const I3CIntegral& integral) const;
    
    /// Adds call tree for ket horizontal recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param skterms The set of filtered geometrical terms.
    /// @param integral The base two center integral.
    void _add_hrr_call_tree(      VCodeLines&  lines,
                            const SG3Terms&    skterms,
                            const I3CIntegral& integral) const;
    
    /// Gets arguments list for ket horizontal recursion function call.
    /// @param skterms The set of filtered geometrical terms.
    /// @param term The base integral term.
    std::string _get_hrr_arguments(const SG3Terms& skterms,
                                   const G3Term&   term) const;
    
    /// Gets arguments list for ket horizontal recursion function call.
    /// @param skterms The set of filtered geometrical terms.
    /// @param term The base integral term.
    /// @param icomponent The index of geometrical derivative on bra side.
    std::string _get_hrr_arguments(const SG3Terms& skterms,
                                   const G3Term&   term,
                                   const int       icomponent) const; 
    
    /// Adds call tree for bra horizontal recursion.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param cterms The set of filtered geometrical terms.
    /// @param integral The base two center integral.
    void _add_bra_geom_call_tree(      VCodeLines&  lines,
                                 const SG3Terms&    cterms,
                                 const I3CIntegral& integral) const;
    
    /// Gets arguments list for bra recursion function call.
    /// @param term The recursion term.
    /// @param cterms The set of filtered geometrical terms.
    std::string _get_bra_geom_arguments(const G3Term&  term,
                                        const SG3Terms& cterms) const;
    
public:
    /// Creates a two-center compute function body generator.
    T3CGeomFuncBodyDriver() = default;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param cterms The set of filtered geometrical terms.
    /// @param skterms The set of filtered geometrical terms.
    /// @param vrr_integrals The set of unique integrals for vertical recursion.
    /// @param integral The base four center integral.
    void write_func_body(      std::ofstream& fstream,
                         const SG3Terms&      cterms,
                         const SG3Terms&      skterms,
                         const SI3CIntegrals& vrr_integrals,
                         const I3CIntegral&   integral) const;
};

#endif /* t3c_geom_body_hpp */
