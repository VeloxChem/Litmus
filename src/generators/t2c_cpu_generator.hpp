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

#ifndef t2c_cpu_generator_hpp
#define t2c_cpu_generator_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>


#include "operator.hpp"
#include "tensor_component.hpp"
#include "t2c_defs.hpp"

// Two-center integrals code generator for CPU.
class T2CCPUGenerator
{
    /// Checks if recursion is available for two-center inetgral with given label.
    /// @param label The label of requested two-center integral.
    bool _is_available(const std::string& label) const;
    
    /// Gets two-center inetgral with requested label.
    /// @param label The label of requested two-center integral.
    /// @param ang_a The angular momentum of center A.
    /// @param ang_b The angular momentum of center B.
    /// @param op_gdrv The geometrical derivative of operator. 
    /// @return The two-center integral.
    I2CIntegral _get_integral(const std::string& label,
                              const int          ang_a,
                              const int          ang_b,
                              const int          op_gdrv = 0) const;
    
    /// Gets file name of file with recursion functions for two center integral.
    /// @param integral The base two center integral.
    /// @return The file name.
    std::string _file_name(const I2CIntegral& integral) const;
    
    /// Writes header file for recursion.
    /// @param integral The base four center integral.
    void _write_cpp_header(const I2CIntegral& integral) const;
    
    /// Writes C++ code file for recursion.
    /// @param integral The base four center integral.
    void _write_cpp_file(const I2CIntegral& integral) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_defines(      std::ofstream& fstream,
                            const I2CIntegral&   integral,
                            const bool           start) const;
    
    /// Writes definitions of includes for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_hpp_includes(      std::ofstream& fstream,
                             const I2CIntegral&   integral) const;
    
    /// Writes definitions of includes for C++ code file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_cpp_includes(      std::ofstream& fstream,
                             const I2CIntegral&   integral) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const I2CIntegral&   integral,
                          const bool           start) const;
        
    /// Writes primitive functions documentation and declaration to header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_prim_funcs_to_cpp_header(      std::ofstream& fstream,
                                         const I2CIntegral&   integral) const;
    
    /// Writes primitive function implementation to C++ code file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_prim_funcs_to_cpp_file(      std::ofstream& fstream,
                                       const I2CIntegral&   integral) const;
    
public:
    /// Creates an electron repulsion integrals CPU code generator.
    T2CCPUGenerator() = default;
     
    /// Generates selected one-electron integrals up to given angular momentum (inclusive) )on A and B centers.
    /// @param label The label of requested two-center integral.
    /// @param angmom The maximum angular momentum of A and B centers.
    /// @param op_gdrv The geometrical derivative of operator. 
    void generate(const std::string& label,
                  const int          angmom,
                  const int          op_gdrv) const;
};

#endif /* t2c_cpu_generator_hpp */
