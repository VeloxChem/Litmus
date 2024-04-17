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

#ifndef t2c_cpu_generators_hpp
#define t2c_cpu_generators_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <utility>

#include "t2c_defs.hpp"

// Two-center integrals code generator for CPU.
class T2CCPUGenerator
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
    
    /// Gets file name of file with recursion functions for two center integral.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @return The file name.
    std::string _file_name(const I2CIntegral&           integral,
                           const std::pair<bool, bool>& rec_form) const;
    
    /// Writes header file for recursion.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    void _write_cpp_header(const I2CIntegral&           integral,
                           const std::pair<bool, bool>& rec_form) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_defines(      std::ofstream&         fstream,
                            const I2CIntegral&           integral,
                            const std::pair<bool, bool>& rec_form,
                            const bool                   start) const;
    
    /// Writes definitions of includes for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    void _write_hpp_includes(      std::ofstream&         fstream,
                             const I2CIntegral&           integral,
                             const std::pair<bool, bool>& rec_form) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const I2CIntegral&   integral,
                          const bool           start) const;
    
public:
    /// Creates a two-center integrals CPU code generator.
    T2CCPUGenerator() = default;
     
    /// Generates selected two-center integrals up to given angular momentum (inclusive)  on A and B centers.
    /// @param label The label of requested two-center integral.
    /// @param max_ang_mom The maximum angular momentum of A and B centers.
    /// @param geom_drvs The geometrical derivative of bra side, integrand, and  ket side.
    /// @param rec_form The recursion form for two center integrals (summation, convolution flags).
    void generate(const std::string&           label,
                  const int                    max_ang_mom,
                  const std::array<int, 3>&    geom_drvs,
                  const std::pair<bool, bool>& rec_form) const;
};

#endif /* t2c_cpu_generators_hpp */