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

#ifndef t1c_cpu_generator_hpp
#define t1c_cpu_generator_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>

#include "operator.hpp"
#include "tensor_component.hpp"
#include "t2c_defs.hpp"
#include "file_stream.hpp"

// One-center GTOs code generator for CPU.
class T1CCPUGenerator
{
   
    /// Gets two-center inetgral with requested label.
    /// @param label The label of requested two-center integral.
    /// @param ang_a The angular momentum of center A.
    /// @param bra_gdrv The geometrical derivative of bra side.
    /// @return The two-center integral.
    I2CIntegral _get_integral(const std::string& label,
                              const int          ang_a,
                              const int          bra_gdrv) const;
    
    /// Gets file name of file with recursion functions for two center integral.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @return The file name.
    std::string _file_name(const int angmom) const;
    
    /// Writes header file for recursion.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    void _write_cpp_header(const int angmom,
                           const int gdrv) const;
    
    /// Writes C++ code file for recursion.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    void _write_cpp_file(const int angmom,
                         const int gdrv) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_defines(      std::ofstream& fstream,
                            const int            angmom,
                            const bool           start) const;
    
    /// Writes definitions of includes for header file.
    /// @param fstream the file stream.
    void _write_hpp_includes(std::ofstream& fstream) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const bool           start) const;
    
    /// Writes definitions of includes for C++ code file.
    /// @param fstream the file stream.
    /// @param angmom The maximum angular momentum of GTOs.
    void _write_cpp_includes(      std::ofstream& fstream,
                             const int            angmom) const;
    
public:
    /// Creates an one-center GTOs  CPU code generator.
    T1CCPUGenerator() = default;
     
    /// Generates selected GTOs values up to given angular momentum (inclusive).
    /// @param label The label of requested two-center integral.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    void generate(const std::string& label,
                  const int          angmom,
                  const int          gdrv) const;
};

#endif /* t1c_cpu_generator_hpp */
