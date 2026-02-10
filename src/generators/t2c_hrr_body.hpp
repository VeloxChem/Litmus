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

#ifndef t2c_hrr_body_hpp
#define t2c_hrr_body_hpp

#include <string>
#include <array>
#include <vector>
#include <utility>
#include <fstream>

#include "t2c_defs.hpp"
#include "file_stream.hpp"

// Two-center compute function body generators for CPU.
class T2CHRRFuncBodyDriver
{
    /// Generates vector of factor strings.
    /// @param integral The base two center integral.
    /// @return The vector of factor strings.
    std::vector<std::string> _get_factors_str(const I2CIntegral& integral) const;
    
    /// Generates vector of buffer strings.
    /// @param rec_dists The vector of recursion distributions.
    /// @param integral The base two center integral.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_buffers_str(const std::vector<R2CDist>& rec_dists,
                                              const I2CIntegral&          integral) const;
    
    /// Checks if integral is needed by recursion.
    /// @param rec_dists The vector of recursion distributions.
    /// @param integral The base two center integral.
    /// @return True if integral component is actually used in simplified recursion, False otherwise.
    bool _find_integral(const std::vector<R2CDist>& rec_dists,
                        const T2CIntegral&          integral) const;
    
    /// Generates vector of buffer strings.
    /// @param integral The base two center integral.
    /// @param components The vector of integral components.
    /// @param rec_range The recursion range [first, last) in integral components space.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_buffers_str(const I2CIntegral&        integral,
                                              const VT2CIntegrals&      components,
                                              const std::array<int, 2>& rec_range) const;
    
    /// Gets tensor label for integral.
    /// @param integral The base two center integral.
    /// @return The tensorial label.
    std::string _get_tensor_label(const I2CIntegral& integral) const;
    
    /// Gets tensor label for integral.
    /// @param integral The base two center integral.
    /// @return The tensorial label.
    std::string _get_tensor_label(const T2CIntegral& integral) const;
    
    /// Adds single loop computation of primitive integrals.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    /// @param components The vector of integral components.
    /// @param rec_range The recursion range [first, last) in integral components space.
    void _add_recursion_loop(      VCodeLines&         lines,
                             const I2CIntegral&        integral,
                             const VT2CIntegrals&      components,
                             const std::array<int, 2>& rec_range) const;
    
    
    /// Adds single loop computation of primitive integrals.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param rec_distributions The recursion distributions.
    void _get_factor_lines(                VCodeLines& lines,
                           const std::vector<R2CDist>& rec_distributions) const;
    
    /// Gets pragma string for vector of recursion distributions.
    /// @param integral The base two center integral.
    /// @param rec_distributions The recursion distributions.
    std::string _get_pragma_str(const I2CIntegral&          integral,
                                const std::vector<R2CDist>& rec_distributions) const;
    
    /// Computes HRR recursion for integral component.
    /// @param integral The base two center integral component.
    /// @return The recursion expansion of integral component.
    R2CDist _get_hrr_recursion(const T2CIntegral& integral) const;
    
    /// Creates code line for recursion expansion.
    /// @param rec_distribution The recursion distribution
    /// @return The string with code line.
    std::string _get_code_line(const R2CDist& rec_distribution) const;
    
    /// Creates code string for recursion term.
    /// @param rec_term The recursion distribution.
    /// @param is_first The flag to indicate first term in recursion expnasion.
    /// @return The string with code term.
    std::string _get_rterm_code(const R2CTerm& rec_term,
                                const bool     is_first) const;
    
    /// Gets integral component label.
    /// @param integral The base two center integral component.
    /// @return The string with integral component label.
    std::string _get_component_label(const T2CIntegral& integral) const;
    
    /// Checks if R(PA) distances are needed for recursion implementation.
    /// @param integral The base two center integral component.
    /// @return True if R(PA) distances are needed, false otherwise.
    bool _need_distances_pa(const I2CIntegral& integral) const;
    
    /// Checks if R(PB) distances are needed for recursion implementation.
    /// @param integral The base two center integral component.
    /// @return True if R(PB) distances are needed, false otherwise.
    bool _need_distances_pb(const I2CIntegral& integral) const;
    
    /// Checks if R(PC) distances are needed for recursion implementation.
    /// @param integral The base two center integral component.
    /// @return True if R(PC) distances are needed, false otherwise.
    bool _need_distances_pc(const I2CIntegral& integral) const;
    
    /// Checks if GTOs exponents are needed for recursion implementation.
    /// @param integral The base two center integral component.
    /// @return True if GTOs exponents are needed, false otherwise.
    bool _need_exponents(const I2CIntegral& integral) const;

public:
    /// Creates a two-center compute function body generator.
    T2CHRRFuncBodyDriver() = default;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void write_func_body(      std::ofstream& fstream,
                         const I2CIntegral&   integral) const;
};

#endif /* t2c_hrr_body_hpp */
