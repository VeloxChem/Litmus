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

#ifndef t2c_ecp_body_hpp
#define t2c_ecp_body_hpp

#include <string>
#include <vector>
#include <utility>
#include <fstream>

#include "t2c_defs.hpp"
#include "file_stream.hpp"

// Two-center ECP compute function body generators for CPU.
class T2CECPFuncBodyDriver
{
    /// Generates vector of strings with GTOs definitions in compute function.
    /// @return The vector of strings with GTOS definitions in compute function.
    std::vector<std::string> _get_gtos_def() const;
    
    /// Generates vector of ket factors in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of ket factors in compute function.
    std::vector<std::string> _get_ket_variables_def(const I2CIntegral& integral) const;
    
    /// Generates vector of buffers in compute function.
    /// @param hrr_integrals The set of inetrgals in horizontal recursion.
    /// @param vrr_integrals The set of inetrgals in vertical recursion.
    /// @param integral The base two center integral.
    /// @return The vector of buffers in compute function.
    std::vector<std::string> _get_buffers_def(const SI2CIntegrals& hrr_integrals,
                                              const SI2CIntegrals& vrr_integrals,
                                              const I2CIntegral&   integral) const;
    
    /// Filters set of VRR integrals requiring contraction.
    /// @param integrals The set of inetrgals in vertical recursion.
    /// @param integral The base two center integral.
    /// @return The set of VRR integrals requiring contraction.
    SI2CIntegrals _filter_contracted(const SI2CIntegrals& integrals,
                                     const I2CIntegral&   integral) const;
    
    /// Adds loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_loop_start(      VCodeLines&  lines,
                         const I2CIntegral& integral) const;
    
    /// Adds loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_loop_end(      VCodeLines&  lines,
                       const I2CIntegral& integral) const;
    
    /// Adds ket loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integral The base two center integral.
    void _add_ket_loop_start(      VCodeLines&  lines,
                             const I2CIntegral& integral) const;
    
    /// Adds ket loop end definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param integrals The set of inetrgals.
    /// @param integral The base two center integral.
    void _add_ket_loop_end(      VCodeLines&    lines,
                           const SI2CIntegrals& integrals,
                           const I2CIntegral&   integral) const;
    
    /// Checks if distances of (R-A) are required for integration.
    /// @param integral The base two center integral.
    bool _need_distances_ra(const I2CIntegral& integral) const;
    
    /// Checks if distances of (R-B) are required for integration.
    /// @param integral The base two center integral.
    bool _need_distances_rb(const I2CIntegral& integral) const;
    
public:
    /// Creates a two-center ECP compute function body generator.
    T2CECPFuncBodyDriver() = default;
    
    /// Writes body of local ECP compute function.
    /// @param fstream the file stream.
    /// @param hrr_integrals The set of inetrgals in horizontal recursion.
    /// @param vrr_integrals The set of inetrgals in vertical recursion.
    /// @param integral The base two center integral.
    void write_func_body(      std::ofstream&         fstream,
                         const SI2CIntegrals&         hrr_integrals,
                         const SI2CIntegrals&         vrr_integrals,
                         const I2CIntegral&           integral) const;
};

#endif /* t2c_ecp_body_hpp */
