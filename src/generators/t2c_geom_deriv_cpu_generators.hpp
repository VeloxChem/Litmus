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

#ifndef t2c_geom_deriv_cpu_generators_hpp
#define t2c_geom_deriv_cpu_generators_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <utility>

#include "t2c_defs.hpp"

// Geometrical derivatives of two-center integrals code generator for CPU.
class T2CGeomDerivCPUGenerator
{
    /// Gets four-center inetgral with requested label.
    /// @param ang_moms The angular momentum of  A, B, C, and D centers.
    /// @param geom_drvs The geometrical derivative of bra and  ket sides.
    /// @return The four-center integral.
    I2CIntegral _get_integral(const std::array<int, 2>& ang_moms,
                              const std::array<int, 3>& geom_drvs) const;
    
    /// Generates set of integrals required for geometrical derivatives.
    /// @param integral The base four center integral.
    /// @return The set of integrals.
//    SI4CIntegrals _generate_geom_integral_group(const I4CIntegral& integral) const;

    /// Writes header file for recursion.
    /// @param geom_integrals The set of unique integrals for geometrical recursion.
    /// @param integral The base two center integral.
    /// @param geom_drvs The geometrical derivative of bra and  ket sides.
    void _write_cpp_header(const SI2CIntegrals&      geom_integrals,
                           const I2CIntegral&        integral,
                           const std::array<int, 3>& geom_drvs) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param geom_drvs The geometrical derivative of bra and  ket sides.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_defines(      std::ofstream&      fstream,
                            const I2CIntegral&        integral,
                            const std::array<int, 3>& geom_drvs,
                            const bool                start) const;
    
    /// Writes definitions of includes for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param geom_drvs The geometrical derivative of bra and  ket sides.
    void _write_hpp_includes(      std::ofstream& fstream,
                             const I2CIntegral&   integral,
                             const std::array<int, 3>& geom_drvs) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const bool           start) const;
    
    /// Writes C++ code file for primtive recursion.
    /// @param geom_integrals The set of unique integrals for geometrical recursion.
    /// @param integral The base four center integral.
    /// @param geom_drvs The geometrical derivative of bra and  ket sides.
    void _write_cpp_file(const SI2CIntegrals&      geom_integrals,
                         const I2CIntegral&        integral,
                         const std::array<int, 3>& geom_drvs) const;
    
    /// Writes definitions of includes for primitive header file.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
//    void _write_cpp_includes(      std::ofstream& fstream,
//                             const I4CIntegral&  integral) const;
    
public:
    /// Creates a geometrical derivatives of two-center integrals CPU code generator.
    T2CGeomDerivCPUGenerator() = default;
     
    /// Generates selected two-center integrals up to given angular momentum (inclusive)  on A and B centers.
    /// @param max_ang_mom The maximum angular momentum of A and B centers.
    /// @param geom_drvs The geometrical derivative of bra and  ket sides.
    void generate(const int                 max_ang_mom,
                  const std::array<int, 3>& geom_drvs) const;
};

#endif /* t2c_geom_deriv_cpu_generators_hpp */
