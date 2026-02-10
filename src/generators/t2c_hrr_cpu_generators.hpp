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

#ifndef t2c_hrr_cpu_generators_hpp
#define t2c_hrr_cpu_generators_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <utility>

#include "t2c_defs.hpp"

// Horizontal recursion of two-center integrals code generator for CPU.
class T2CHRRCPUGenerator
{
    
    /// Gets two-center inetgral with requested label.
    /// @param ang_moms The angular momentum of  A and B centers.
    /// @return The two-center integral.
    I2CIntegral _get_integral(const std::array<int, 2>& ang_moms) const;
    
    /// Writes horizontal header file for recursion.
    /// @param integral The base two center integral.
    void _write_hrr_cpp_header(const I2CIntegral& integral) const;
    
    /// Writes C++ code file for horizontal recursion.
    /// @param integral The base two center integral.
    void _write_hrr_cpp_file(const I2CIntegral& integral) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_defines(      std::ofstream& fstream,
                            const I2CIntegral&   integral,
                            const bool           start) const;
    
    /// Writes definitions of includes for horizontal header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_hrr_hpp_includes(      std::ofstream& fstream,
                                  const I2CIntegral&   integral) const;
    
    /// Writes definitions of includes for horizontal header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_hrr_cpp_includes(      std::ofstream& fstream,
                                  const I2CIntegral&  integral) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const I2CIntegral&   integral,
                          const bool           start) const;
    
public:
    /// Creates horizontal recursion of two-center integrals CPU code generator.
    T2CHRRCPUGenerator() = default;
     
    /// Generates selected two-center integrals up to given angular momentum on A, B centers.
    /// @param max_ang_mom The maximum angular momentum of A, B, C and D centers.
    void generate(const int max_ang_mom) const;
};


#endif /* t2c_hrr_cpu_generators_hpp */
