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

#ifndef t4c_eri_tree_generators_hpp
#define t4c_eri_tree_generators_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <utility>

#include "t4c_defs.hpp"

// Geometrical derivatives of four-center integrals code generator for CPU.
class T4CCallTreeGenerator
{
    /// Checks if recursion is available for four-center inetgral with given label.
    /// @param label The label of requested four-center integral.
    bool _is_available(const std::string& label) const;
    
    /// Gets four-center inetgral with requested label.
    /// @param label The label of requested four-center integral.
    /// @param ang_moms The angular momentum of  A, B, C, and D centers.
    /// @return The four-center integral.
    I4CIntegral _get_integral(const std::string&        label,
                              const std::array<int, 4>& ang_moms) const;
    
    /// Gets file name of file with recursion functions for two center integral.
    /// @param integral The base two center integral.
    /// @return The file name.
    std::string _file_name(const I4CIntegral& integral) const;
        
public:
    /// Creates a geometrical derivatives of four-center integrals CPU code generator.
    T4CCallTreeGenerator() = default;
     
    /// Generates call tree for selected four-center integrals up to given angular momentum (inclusive)  on A, B, C, and D centers.
    /// @param label The label of requested two-center integral.
    /// @param max_ang_mom The maximum angular momentum of A, B, C and D centers.
    void generate(const std::string& label,
                  const int          max_ang_mom) const;
};

#endif /* t4c_eri_tree_generators_hpp */
