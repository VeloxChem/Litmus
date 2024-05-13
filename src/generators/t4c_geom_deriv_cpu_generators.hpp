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

#ifndef t4c_geom_deriv_cpu_generators_hpp
#define t4c_geom_deriv_cpu_generators_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <utility>

#include "t4c_defs.hpp"

// Geometrical derivatives of four-center integrals code generator for CPU.
class T4CGeomDerivCPUGenerator
{
   
    /// Gets four-center inetgral with requested label.
    /// @param ang_moms The angular momentum of  A, B, C, and D centers.
    /// @param geom_drvs The geometrical derivative of bra side, integrand, and  ket side.
    /// @return The four-center integral.
    I4CIntegral _get_integral(const std::array<int, 4>& ang_moms,
                              const std::array<int, 4>& geom_drvs) const;
    
    /// Creates recursion group for all components of given integral.
    /// @param components The vector of integral components.
    /// @param integral The diff. integral.
    /// @return The recursion group.
    R4Group _generate_integral_group(const VT4CIntegrals& components,
                                     const I4CIntegral&   integral) const;


    /// Writes header file for recursion.
    /// @param integral The base two center integral.
    void _write_cpp_header(const I4CIntegral& integral) const;
    
public:
    /// Creates a geometrical derivatives of four-center integrals CPU code generator.
    T4CGeomDerivCPUGenerator() = default;
     
    /// Generates selected four-center integrals up to given angular momentum (inclusive)  on A, B, C, and D centers.
    /// @param max_ang_mom The maximum angular momentum of A, B, C and D centers.
    /// @param max_geom_order The maximum order of geometrical derrivatives.
    void generate(const int max_ang_mom,
                  const int max_geom_order) const;
};


#endif /* t4c_geom_deriv_cpu_generators_hpp */
