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
// limitations under the License. by Zilvinas Rinkevicius on 2025-01-22.
//

#ifndef t3c_hrr_body_hpp
#define t3c_hrr_body_hpp

#include <string>
#include <array>
#include <vector>
#include <utility>
#include <fstream>

#include "t3c_defs.hpp"
#include "file_stream.hpp"

// Three-center compute function body generators for CPU.
class T3CHrrFuncBodyDriver
{
    /// Generates vector of buffer strings.
    /// @param integral The base four center integral.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_ket_buffers_str(const std::vector<R3CDist>& rec_dists,
                                                  const I3CIntegral&          integral) const;
    
    /// Generates vector of buffer strings.
    /// @param integral The base four center integral.
    /// @param components The vector of integral components.
    /// @param rec_range The recursion range [first, last) in integral components space.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_ket_buffers_str(const I3CIntegral&        integral,
                                                  const VT3CIntegrals&      components,
                                                  const std::array<int, 2>& rec_range) const;
    
    /// Checks if integral is needed by recursion.
    /// @param rec_dists The vector of recursion distributions.
    /// @param integral The base two center integral.
    /// @return True if integral component is actually used in simplified recursion, False otherwise.
    bool _find_integral(const std::vector<R3CDist>& rec_dists,
                        const T3CIntegral&          integral) const;
    
    /// Gets tensor label for integral.
    /// @param integral The base four center integral.
    /// @return The tensorial label.
    std::string _get_tensor_label(const I3CIntegral& integral) const;
    
    /// Gets tensor label for integral.
    /// @param integral The base four center integral.
    /// @return The tensorial label.
    std::string _get_tensor_label(const T3CIntegral& integral) const;
    
    /// Adds single loop computation of primitive integrals.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base four center integral.
    /// @param components The vector of integral components.
    /// @param rec_range The recursion range [first, last) in integral components space.
    void _add_ket_recursion_loop(      VCodeLines&         lines,
                                 const I3CIntegral&        integral,
                                 const VT3CIntegrals&      components,
                                 const std::array<int, 2>& rec_range) const;
    
    /// Gets pragma string for vector of recursion distributions.
    /// @param integral The base four center integral.
    std::string _get_ket_pragma_str(const I3CIntegral& integral,
                                    const std::vector<R3CDist>& rec_distributions) const;
    
    /// Computes ket horizontal recursion for integral component.
    /// @param integral The base four center integral component.
    /// @return The recursion expansion of integral component.
    R3CDist _get_ket_hrr_recursion(const T3CIntegral& integral) const;
    
    /// Creates code line for recursion expansion.
    /// @param rec_distribution The recursion distribution
    /// @return The string with code line.
    std::string _get_ket_code_line(const R3CDist& rec_distribution) const;
    
    /// Creates code string for recursion term.
    /// @param rec_term The recursion distribution.
    /// @param is_first The flag to indicate first term in recursion expnasion.
    /// @return The string with code term.
    std::string _get_ket_rterm_code(const R3CTerm& rec_term,
                                    const bool     is_first) const;
    
    /// Gets integral component label.
    /// @param integral The base four center integral component.
    /// @return The string with integral component label.
    std::string _get_component_label(const T3CIntegral& integral) const;
    
    /// Gets integral component label.
    /// @param integral The base four center integral component.
    /// @return The string with integral component label.
    std::string _get_ket_component_label(const T3CIntegral& integral) const;
    
    /// Gets integral offset definition.
    /// @param integral The base four center integral component.
    /// @return The string with integral offset definition.
    std::string _get_ket_offset_def(const I3CIntegral& integral) const;
    
    /// Gets integral offset label.
    /// @param integral The base four center integral component.
    /// @return The string with integral offset label.
    std::string _get_ket_offset_label(const I3CIntegral& integral) const;
    
public:
    /// Creates a three-center compute function body generator.
    T3CHrrFuncBodyDriver() = default;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void write_func_body(      std::ofstream& fstream,
                         const I3CIntegral&   integral) const;
};

#endif /* t3c_hrr_body_hpp */
