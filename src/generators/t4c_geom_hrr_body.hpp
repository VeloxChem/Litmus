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

#ifndef t4c_geom_hrr_body_hpp
#define t4c_geom_hrr_body_hpp

#include <string>
#include <array>
#include <vector>
#include <utility>
#include <fstream>

#include "t4c_defs.hpp"
#include "file_stream.hpp"

// Four-center compute function body generators for CPU.
class T4CGeomHrrFuncBodyDriver
{
    /// Computes ket horizontal recursion for integral component.
    /// @param integral The base four center integral component.
    /// @return The recursion expansion of integral component.
    R4CDist _get_ket_hrr_recursion(const T4CIntegral& integral) const;
    
    /// Computes ket horizontal recursion for integral component.
    /// @param integral The base four center integral component.
    /// @return The recursion expansion of integral component.
    R4CDist _get_ket_geom_hrr_recursion(const T4CIntegral& integral) const;
    
    /// Computes ket horizontal recursion for integral component.
    /// @param integral The base four center integral component.
    /// @return The recursion expansion of integral component.
    R4CDist _get_bra_hrr_recursion(const T4CIntegral& integral) const;
    
    /// Generates vector of buffer strings.
    /// @param integral The base four center integral.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_ket_buffers_str(const std::vector<R4CDist>& rec_dists,
                                                  const I4CIntegral&          integral) const;
    
    /// Generates vector of buffer strings.
    /// @param integral The base four center integral.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_ket_geom_buffers_str(const std::vector<R4CDist>& rec_dists,
                                                       const I4CIntegral&          integral) const;
    
    /// Generates vector of buffer strings.
    /// @param integral The base four center integral.
    /// @param components The vector of integral components.
    /// @param rec_range The recursion range [first, last) in integral components space.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_ket_buffers_str(const I4CIntegral&        integral,
                                                  const VT4CIntegrals&      components,
                                                  const std::array<int, 2>& rec_range) const;
    
    /// Generates vector of buffer strings.
    /// @param integral The base four center integral.
    /// @param components The vector of integral components.
    /// @param rec_range The recursion range [first, last) in integral components space.
    /// @param ket_index The index of geometrical derivative.
    /// @param ket_components The numbed of components in integral.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_ket_geom_buffers_str(const I4CIntegral&        integral,
                                                       const VT4CIntegrals&      components,
                                                       const std::array<int, 2>& rec_range,
                                                       const int                 ket_index,
                                                       const int                 ket_components) const;
    
    /// Generates vector of buffer strings.
    /// @param integral The base four center integral.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_bra_buffers_str(const std::vector<R4CDist>& rec_dists,
                                                  const I4CIntegral&          integral) const;
    
    /// Generates vector of buffer strings.
    /// @param integral The base four center integral.
    /// @param components The vector of integral components.
    /// @param rec_range The recursion range [first, last) in integral components space.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_bra_buffers_str(const I4CIntegral&        integral,
                                                  const VT4CIntegrals&      components,
                                                  const std::array<int, 2>& rec_range) const;
    
    /// Checks if integral is needed by recursion.
    /// @param rec_dists The vector of recursion distributions.
    /// @param integral The base two center integral.
    /// @return True if integral component is actually used in simplified recursion, False otherwise.
    bool _find_integral(const std::vector<R4CDist>& rec_dists,
                        const T4CIntegral&          integral) const;
    
    /// Gets tensor label for integral.
    /// @param integral The base four center integral.
    /// @return The tensorial label.
    std::string _get_tensor_label(const T4CIntegral& integral) const;
    
    /// Gets integral component label.
    /// @param integral The base four center integral component.
    /// @return The string with integral component label.
    std::string _get_bra_component_label(const T4CIntegral& integral) const;
    
    /// Gets integral component label.
    /// @param integral The base four center integral component.
    /// @return The string with integral component label.
    std::string _get_full_bra_component_label(const T4CIntegral& integral) const;
    
    /// Gets integral component label.
    /// @param integral The base four center integral component.
    /// @return The string with integral component label.
    std::string _get_ket_component_label(const T4CIntegral& integral) const;
    
    /// Gets integral offset definition.
    /// @param integral The base four center integral component.
    /// @return The string with integral offset definition.
    std::string _get_bra_offset_def(const I4CIntegral& integral) const;
    
    /// Gets integral offset definition.
    /// @param integral The base four center integral component.
    /// @return The string with integral offset definition.
    std::string _get_full_bra_offset_def(const I4CIntegral& integral) const;
    
    /// Gets integral offset definition.
    /// @param integral The base four center integral component.
    /// @return The string with integral offset definition.
    std::string _get_ket_offset_def(const I4CIntegral& integral) const;
    
    /// Gets integral offset label.
    /// @param integral The base four center integral component.
    /// @return The string with integral offset label.
    std::string _get_bra_offset_label(const I4CIntegral& integral) const;
    
    /// Gets integral offset label.
    /// @param integral The base four center integral component.
    /// @return The string with integral offset label.
    std::string _get_full_bra_offset_label(const I4CIntegral& integral) const;
    
    /// Gets integral offset label.
    /// @param integral The base four center integral component.
    /// @return The string with integral offset label.
    std::string _get_ket_offset_label(const I4CIntegral& integral) const;
    
    /// Adds single loop computation of primitive integrals.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base four center integral.
    /// @param components The vector of integral components.
    /// @param rec_range The recursion range [first, last) in integral components space.
    /// @param ket_index The index of geometrical derivative.
    /// @param ket_components The numbed of components in integral.
    void _add_ket_recursion_loop(      VCodeLines&         lines,
                                 const I4CIntegral&        integral,
                                 const VT4CIntegrals&      components,
                                 const std::array<int, 2>& rec_range,
                                 const int                 ket_index,
                                 const int                 ket_components) const;
    
    /// Gets pragma string for vector of recursion distributions.
    /// @param integral The base four center integral.
    std::string _get_ket_pragma_str(const I4CIntegral& integral,
                                    const std::vector<R4CDist>& rec_distributions) const;
    
    /// Creates code line for recursion expansion.
    /// @param rec_distribution The recursion distribution
    /// @return The string with code line.
    std::string _get_ket_code_line(const R4CDist& rec_distribution) const;
    
    /// Creates code string for recursion term.
    /// @param rec_term The recursion distribution.
    /// @param is_first The flag to indicate first term in recursion expnasion.
    /// @return The string with code term.
    std::string _get_ket_rterm_code(const R4CTerm& rec_term,
                                    const bool     is_first) const;
    
    /// Adds single loop computation of primitive integrals.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base four center integral.
    /// @param components The vector of integral components.
    /// @param rec_range The recursion range [first, last) in integral components space.
    void _add_bra_recursion_loop(      VCodeLines&         lines,
                                 const I4CIntegral&        integral,
                                 const VT4CIntegrals&      components,
                                 const std::array<int, 2>& rec_range) const;
    
    /// Gets pragma string for vector of recursion distributions.
    /// @param integral The base four center integral.
    std::string _get_bra_pragma_str(const I4CIntegral& integral,
                                    const std::vector<R4CDist>& rec_distributions) const;
    
    /// Creates code line for recursion expansion.
    /// @param rec_distribution The recursion distribution
    /// @return The string with code line.
    std::string _get_bra_code_line(const R4CDist& rec_distribution) const;
    
    /// Creates code string for recursion term.
    /// @param rec_term The recursion distribution.
    /// @param is_first The flag to indicate first term in recursion expnasion.
    /// @return The string with code term.
    std::string _get_bra_rterm_code(const R4CTerm& rec_term,
                                    const bool     is_first) const;

public:
    
    /// Creates a two-center compute function body generator.
    T4CGeomHrrFuncBodyDriver() = default;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void write_ket_func_body(      std::ofstream& fstream,
                             const I4CIntegral&   integral) const;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void write_ket_geom_func_body(      std::ofstream& fstream,
                                  const I4CIntegral&   integral) const;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void write_bra_func_body(      std::ofstream& fstream,
                             const I4CIntegral&   integral) const;
};

#endif /* t4c_geom_hrr_body_hpp */
