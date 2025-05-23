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

#ifndef g2c_cpu_generators_hpp
#define g2c_cpu_generators_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <utility>

#include "t2c_defs.hpp"

// Two-center integrals on grid code generator for CPU.
class G2CCPUGenerator
{
    /// Checks if recursion is available for two-center inetgral with given label.
    /// @param label The label of requested two-center integral.
    bool _is_available(const std::string& label) const;
    
    /// Gets two-center inetgral with requested label.
    /// @param label The label of requested two-center integral.
    /// @param ang_moms The angular momentum of  A and B centers.
    /// @param geom_drvs The geometrical derivative of bra side, integrand, and  ket side.
    /// @return The two-center integral.
    I2CIntegral _get_integral(const std::string&        label,
                              const std::array<int, 2>& ang_moms,
                              const std::array<int, 3>& geom_drvs) const;
    
    /// Generates set of integrals required for vertical Obara-Saika recursion.
    /// @param integral The base two center integral.
    /// @return The set of integrals.
    SI2CIntegrals _generate_integral_group(const I2CIntegral& integral,
                                           const std::array<int, 3>& geom_drvs) const;
    
    /// Gets file name of file with recursion functions for two center integral.
    /// @param integral The base two center integral.
    /// @param use_rs The flag for use of range-separated Coulomb interactions.
    /// @return The file name.
    std::string _file_name(const I2CIntegral& integral,
                           const bool         use_rs) const;
    
    /// Writes header file for recursion.
    /// @param integrals The set of unique integrals.
    /// @param integral The base two center integral.
    /// @param use_rs The flag for use of range-separated Coulomb interactions.
    void _write_cpp_header(const SI2CIntegrals& integrals,
                           const I2CIntegral&   integral,
                           const bool           use_rs) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param use_rs The flag for use of range-separated Coulomb interactions.
    /// @param is_prim_rec The flag to indicate primitive recurion.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_defines(      std::ofstream&         fstream,
                            const I2CIntegral&           integral,
                            const bool                   use_rs,
                            const bool                   is_prim_rec,
                            const bool                   start) const;
    
    /// Writes definitions of includes for header file.
    /// @param fstream the file stream.
    /// @param integrals The set of unique integrals.
    /// @param integral The base two center integral.
    /// @param use_rs The flag for use of range-separated Coulomb interactions.
    void _write_hpp_includes(      std::ofstream&         fstream,
                             const SI2CIntegrals&         integrals,
                             const I2CIntegral&           integral,
                             const bool                   use_rs) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const I2CIntegral&   integral,
                          const bool           start) const;
    
    /// Writes primitive header file for recursion.
    /// @param integral The base two center integral.
    void _write_prim_cpp_header(const I2CIntegral& integral) const;
    
    /// Writes definitions of includes for primitive header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_prim_hpp_includes(      std::ofstream& fstream,
                                  const I2CIntegral&   integral) const;
    
    /// Writes C++ code file for primtive recursion.
    /// @param integral The base two center integral.
    void _write_prim_cpp_file(const I2CIntegral& integral) const;
    
    /// Writes definitions of includes for primitive header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_prim_cpp_includes(      std::ofstream& fstream,
                                  const I2CIntegral&  integral) const;
    
public:
    /// Creates a two-center integrals on grid CPU code generator.
    G2CCPUGenerator() = default;
     
    /// Generates selected two-center integrals up to given angular momentum (inclusive)  on A and B centers.
    /// @param label The label of requested two-center integral.
    /// @param max_ang_mom The maximum angular momentum of A and B centers.
    /// @param geom_drvs The geometrical derivative of bra side, integrand, and  ket side.
    /// @param use_rs The flag for use of range-separated Coulomb interactions.
    void generate(const std::string&           label,
                  const int                    max_ang_mom,
                  const std::array<int, 3>&    geom_drvs,
                  const bool                   use_rs) const;
};

#endif /* g2c_cpu_generators_hpp */
