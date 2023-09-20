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

#ifndef t4c_diag_docs_hpp
#define t4c_diag_docs_hpp

#include <string>
#include <vector>
#include <fstream>

#include "t4c_defs.hpp"

// Diagonal four-center documentation generator for CPU.
class T4CDiagDocuDriver
{
    /// Generates compute string.
    /// @param integral The base two center integral.
    /// @return The compute string.
    std::string _get_compute_str(const I4CIntegral& integral) const;
    
    /// Generates primitive compute string.
    /// @param component the integral component.
    /// @param integral The base four center integral.
    /// @return The primitive compute string.
    std::string _get_prim_compute_str(const T4CIntegral& component,
                                      const I4CIntegral& integral) const;
    
    /// Generates vector of variable strings.
    /// @return The vector of variable strings.
    std::vector<std::string> _get_vars_str() const;
    
    /// Generates vector of variable strings.
    /// @param diagonal The flag to set diagonal form.
    /// @return The vector of variable strings.
    std::vector<std::string> _get_prim_vars_str(const bool diagonal) const;

public:
    /// Creates a four-center documentation generator.
    T4CDiagDocuDriver() = default;
    
    /// Writes documentation string for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void write_doc_str(      std::ofstream& fstream,
                       const I4CIntegral&   integral) const;
    
    /// Writes documentation string for primitive compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param diagonal The flag to set diagonal form.
    void write_prim_doc_str(      std::ofstream& fstream,
                            const T4CIntegral&   component,
                            const I4CIntegral&   integral,
                            const bool           diagonal) const;
};

#endif /* t4c_diag_docs_hpp */
