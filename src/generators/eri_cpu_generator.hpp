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

#ifndef eri_cpu_generator_hpp
#define eri_cpu_generator_hpp

#include "eri_driver.hpp"
#include "repository.hpp"

/// Electron repulsion integrals code generator for CPU class.
class EriCPUGenerator
{
    // The diagonal form flag.
    bool _diag_form;
    
    /// Writes header file for VRR recursion.
    /// @param integral The base four center integral.
    /// @param rectype The recursion type.
    /// @return The file name.
     std::string _file_name(const I4CIntegral& integral,
                            const std::string& rectype) const;
    
    /// Writes header file for VRR recursion.
    /// @param integral The base four center integral.
    void _write_vrr_cpp_header(const I4CIntegral& integral) const;
    
public:
    /// Creates an electron repulsion integrals CPU code generator.
    EriCPUGenerator();
 
    /// Sets diagonal form of generated integrals.
    void set_diag_form(); 
    
    /// Generates electron repulsion integrals code for the given repository.
    /// @param repo The repository of two-electron integrals.
    void generate(const Repository<R4Group, T4CIntegral>& repo) const;
};

#endif /* eri_cpu_generator_hpp */
