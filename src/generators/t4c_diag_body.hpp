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

#ifndef t4c_diag_body_hpp
#define t4c_diag_body_hpp

#include <string>
#include <vector>
#include <fstream>

#include "t4c_defs.hpp"
#include "file_stream.hpp"

// Diagonal four-center compute function body generators for CPU.
class T4CDiagFuncBodyDriver
{
    /// Generates vector of strings with GTOs pairs definitions in compute function.
    /// @return The vector of strings with GTOS definitions in compute function.
    std::vector<std::string> _get_gtos_def() const;
    
    /// Generates vector of strings with variables definitions in compute function.
    /// @param integral The base four center integral.
    /// @return The vector of strings with variables definitions in compute function.
    std::vector<std::string> _get_vars_def(const I4CIntegral& integral) const;
    
    /// Generates vector of strings with main loop definition in compute function.
    /// @return The vector of strings with main loop definition in compute function.
    std::vector<std::string> _get_batches_def() const;
    
    /// Adds loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    void _add_batches_loop_start(VCodeLines& lines) const;
    
    /// Adds loop body definitions to code lines container.
    /// @param lines The code lines container to which loop body definition are added.
    /// @param integral The base four center integral.
    void _add_batches_loop_body(      VCodeLines&  lines,
                                const I4CIntegral& integral) const;
    
    /// Adds loop end definitions to code lines container.
    /// @param lines The code lines container to which loop end definition are added.
    void _add_batches_loop_end(VCodeLines& lines) const;
    
    /// Adds loop body definitions to code lines container.
    /// @param lines The code lines container to which loop body definition are added.
    /// @param integral The base four center integral.
    /// @param component The base four center integral component.
    void _add_component_body(      VCodeLines&  lines,
                             const I4CIntegral& integral,
                             const T4CIntegral& component) const;
    
public:
    /// Creates a diagonal four-center compute function body generator.
    T4CDiagFuncBodyDriver() = default;
    
    /// Writes body of compute function.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void write_func_body(      std::ofstream& fstream,
                         const I4CIntegral&   integral) const;
};

#endif /* t4c_diag_body_hpp */