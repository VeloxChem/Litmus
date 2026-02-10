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

#ifndef t2c_proj_ecp_cpu_generators_hpp
#define t2c_proj_ecp_cpu_generators_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <utility>

#include "t2c_defs.hpp"

// Two-center integrals code generator for CPU.
class T2CProjECPCPUGenerator
{
    /// Checks if recursion is available for two-center inetgral with given label.
    /// @param label The label of requested two-center integral.
    bool _is_available(const std::string& label) const;
    
    /// Gets two-center inetgral with requested label.
    /// @param label The label of requested two-center integral.
    /// @param ang_moms The angular momentum of  A and B centers.
    /// @param proj_ang_mom The maximum angular momentum of projector on center C.
    /// @return The two-center integral.
    M2Integral _get_integral(const std::string&        label,
                             const std::array<int, 2>& ang_moms,
                             const int                 proj_ang_mom) const;
    
    /// Generates set of integrals required for vertical Obara-Saika recursion.
    /// @param integral The base two center integral.
    /// @return The set of integrals.
    SM2Integrals _generate_integral_group(const M2Integral& integral) const;
    
public:
    /// Creates a two-center integrals CPU code generator.
    T2CProjECPCPUGenerator() = default;
     
    /// Generates selected two-center integrals up to given angular momentum (inclusive)  on A and B centers.
    /// @param label The label of requested two-center integral.
    /// @param max_ang_mom The maximum angular momentum of A and B centers.
    /// @param proj_ang_mom The maximum angular momentum of projector on center C.
    void generate(const std::string& label,
                  const int          max_ang_mom,
                  const int          proj_ang_mom) const;
};

#endif /* t2c_proj_ecp_cpu_generators_hpp */
