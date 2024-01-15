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

#ifndef c2c_body_hpp
#define c2c_body_hpp

#include <string>
#include <vector>
#include <array>
#include <set>
#include <fstream>

#include "t2c_defs.hpp"
#include "t2c_utils.hpp"
#include "file_stream.hpp"

// Two-center compute function body generators for CPU.
class C2CFuncBodyDriver
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
    /// @param rgroup The recursion group.
    /// @param integral The base two center integral.
    /// @return The vector of buffer definitions in compute function.
    std::vector<std::string> _get_buffers_def(const R2Group&     rgroup,
                                              const I2CIntegral& integral) const;
    
    /// Generates vector of auxilaries definitions in compute function.
    /// @param rgroup The recursion group.
    /// @return The vector of auxilaries definitions in compute function.
    std::vector<std::string> _get_auxilaries_def(const R2Group& rgroup) const;
    
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
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_bra_loop_start(      VCodeLines&  lines,
                             const I2CIntegral& integral,
                             const bool         diagonal) const;
    
    /// Adds bra loop body definitions to code lines container.
    /// @param lines The code lines container to which bra loop body definition are added.
    /// @param rgroup The recursion group.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_bra_loop_body(      VCodeLines&  lines,
                            const R2Group&     rgroup,
                            const I2CIntegral& integral,
                            const bool         sum_form,
                            const bool         diagonal) const;
    
    /// Adds bra loop end definitions to code lines container.
    /// @param lines The code lines container to which bra loop end definition are added.
    void _add_bra_loop_end(VCodeLines& lines) const;
    
    /// Adds code block in bra loop definition to code lines container.
    /// @param lines The code lines container to which bra loop body definition are added.
    /// @param rgroup The recursion group.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @param first The first recursion expansion in code block.
    /// @param last The last recursion expansion in code block.
    void _add_bra_loop_block(      VCodeLines&  lines,
                             const R2Group&     rgroup,
                             const I2CIntegral& integral,
                             const bool         sum_form,
                             const bool         diagonal,
                             const size_t       first,
                             const size_t       last) const;
    
    /// Adds code line in bra loop definition to code lines container.
    /// @param lines The code lines container to which bra loop body definition are added.
    /// @param rdist The recursion distribution.
    /// @param integral The base two center integral.
    /// @param auxilaries The set of unique auxilaries (n,mt).
    /// @param index The index of integrals buffer.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    void _add_bra_loop_line(      VCodeLines&   lines,
                            const R2CDist&      rdist,
                            const I2CIntegral&  integral,
                            const V4Auxilaries& auxilaries,
                            const size_t        index,
                            const bool          sum_form) const;
    
    /// Converts recursion term into it's code representation.
    /// @param rterm The recursion term.
    /// @param auxilaries The set of unique auxilaries (n,mt).
    /// @param is_first The flag indixating if recursion term is first in recursion expansion.
    std::string _get_rterm_code(const R2CTerm&      rterm,
                                const V4Auxilaries& auxilaries,
                                const bool          is_first) const;
    
    /// Adds loop prefactors for range of expansions in recursion group.
    /// @param lines The code lines container to which bra loop body definition are added.
    /// @param rgroup The recursion group.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @param first The first recursion expansion in code block.
    /// @param last The last recursion expansion in code block.
    void _add_loop_prefactors(      VCodeLines& lines,
                              const R2Group&    rgroup,
                              const bool        sum_form,
                              const size_t      first,
                              const size_t      last) const;
    
    /// Adds distribution code in bra loop definition to code lines container.
    /// @param lines The code lines container to which bra loop body definition are added.
    /// @param rgroup The recursion group.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @param first The first recursion expansion in code block.
    /// @param last The last recursion expansion in code block.
    void _write_block_distributor(      VCodeLines&  lines,
                                  const R2Group&     rgroup,
                                  const I2CIntegral& integral,
                                  const bool         sum_form,
                                  const bool         diagonal,
                                  const size_t       first,
                                  const size_t       last) const;
    
    
    /// Gets matrix label for integral component.
    /// @param integral The two center integral component.
    /// @return The label of integral component.
    std::string _get_matrix_label(const T2CIntegral& integral) const;
    
    /// Gets dimensions of code block.
    /// @return The size of code block.
    size_t _get_block_size() const;
    
    /// Checks if auxilaries are needed for this integral.
    /// @param integral The two center integral component.
    /// @return True if auxilaries are needed, false otherwise..
    bool _need_auxilaries(const I2CIntegral& integral) const;
   
public:
    /// Creates a two-center compute function body generator.
    C2CFuncBodyDriver() = default;
    
    /// Writes body of compute function.
    /// @param fstream the file stream.
    /// @param rgroup the recursion group.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void write_func_body(      std::ofstream& fstream,
                         const R2Group&       rgroup,
                         const I2CIntegral&   integral,
                         const bool           sum_form,
                         const bool           diagonal) const;
};

#endif /* c2c_body_hpp */
