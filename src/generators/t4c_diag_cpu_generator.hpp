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

#ifndef t4c_diag_cpu_generator_hpp
#define t4c_diag_cpu_generator_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>

#include "operator.hpp"
#include "tensor_component.hpp"
#include "t2c_defs.hpp"
#include "file_stream.hpp"

// Four-center diagonal integrals code generator for CPU.
class T4CDiagCPUGenerator
{
    /// Checks if recursion is available for four-center diagonal inetgral with given label.
    /// @param label The label of requested two-center integral.
    bool _is_available(const std::string& label) const;
    
public:
    /// Creates an electron repulsion integrals CPU code generator.
    T4CDiagCPUGenerator() = default;
     
    /// Generates selected one-electron integrals up to given angular momentum (inclusive) )on A and B centers.
    /// @param label The label of requested two-center integral.
    /// @param angmom The maximum angular momentum of A and B centers.
    void generate(const std::string& label,
                  const int          angmom) const;
};

#endif /* t4c_diag_cpu_generator_hpp */
