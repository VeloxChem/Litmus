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

#ifndef c2c_auxilary_body_hpp
#define c2c_auxilary_body_hpp

#include <string>
#include <vector>
#include <array>
#include <set>
#include <fstream>

#include "t2c_defs.hpp"
#include "t2c_utils.hpp"
#include "file_stream.hpp"

// Two-center auxilary compute function body generators for CPU.
class C2CAuxilaryBodyDriver
{
    /// Generates vector of strings with math constant definitions in compute function.
    /// @return The vector of strings with math constant definitions in compute function.
    std::vector<std::string> _get_math_def() const;
    
    /// Generates vector of strings with GTOs definitions in compute function.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of strings with GTOS definitions in compute function.
    std::vector<std::string> _get_gtos_def(const bool diagonal) const;
    
    /// Generates vector of ket factors in compute function.
    /// @return The vector of ket factors in compute function.
    std::vector<std::string> _get_ket_variables_def() const;
    
    /// Generates vector of auxilaries definitions in compute function.
    /// @param rgroup The recursion group.
    /// @return The vector of auxilaries definitions in compute function.
    std::vector<std::string> _get_auxilaries_def(const R2Group& rgroup) const;
    
    /// Generates vector of strings with coordinates definitions in bra side.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of strings with coordinates definitions.
    std::vector<std::string> _get_bra_coords(const bool diagonal) const;
    
    /// Generates vector of strings with coordinates definitions in ket side.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of strings with coordinates definitions.
    std::vector<std::string> _get_ket_coords(const bool diagonal) const;
    
    /// Generates vector of pointers to GTOs data on ket side.
    /// @return The vector of ket factors in compute function.
    std::vector<std::string> _get_ket_pointers_def() const;
    
    /// Adds primitives loop start definitions to code lines container.
    /// @param lines The code lines container to which bra loop start definition are added.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void _add_prim_loop_start(      VCodeLines& lines,
                             const bool         diagonal) const;
    
    /// Adds auxilary loop body to compute function.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param rgroup The recursion group.
    /// @param integral The base two center integral.
    void _add_aux_loop_body(      VCodeLines&  lines,
                            const R2Group&     rgroup,
                            const I2CIntegral& integral) const;
    
    /// Adds primitives loop end definitions to code lines container.
    /// @param lines The code lines container to which bra loop end definition are added.
    void _add_prim_loop_end(VCodeLines& lines) const;
    
    /// Adds overlap factors to code lines container.
    /// @param lines The code lines container to which bra loop end definition are added.
    /// @param auxilaries The set of unique auxilaries (n,m,t).
    void _add_aux_overlap_factor(      VCodeLines&   lines,
                                 const V4Auxilaries& auxilaries) const;
    
    /// Adds polynomial factors to code lines container.
    /// @param lines The code lines container to which bra loop end definition are added.
    /// @param auxilaries The set of unique auxilaries (n,m,t).
    void _add_aux_polynomial_factors(      VCodeLines&   lines,
                                     const V4Auxilaries& auxilaries) const;
    
    /// Adds auxilaries values computation to code lines container.
    /// @param lines The code lines container to which bra loop end definition are added.
    /// @param auxilaries The set of unique auxilaries (n,m,t).
    void _add_aux_values(      VCodeLines&   lines,
                         const V4Auxilaries& auxilaries) const;
    
public:
    /// Creates a two-center auxilary compute function body generator.
    C2CAuxilaryBodyDriver() = default;
    
    /// Writes body of auxilary compute function.
    /// @param fstream the file stream.
    /// @param rgroup the recursion group.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void write_aux_body(      std::ofstream& fstream,
                        const R2Group&       rgroup,
                        const I2CIntegral&   integral,
                        const bool           diagonal) const;
};

#endif /* c2c_auxilary_body_hpp */
