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

#ifndef t3c_cpu_generators_hpp
#define t3c_cpu_generators_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <utility>

#include "t3c_defs.hpp"

// Three-center integrals code generator for CPU.
class T3CCPUGenerator
{
    /// Checks if recursion is available for three-center inetgral with given label.
    /// @param label The label of requested four-center integral.
    bool _is_available(const std::string& label) const;
    
    /// Gets three-center inetgral with requested label.
    /// @param label The label of requested four-center integral.
    /// @param ang_moms The angular momentum of  A and B centers.
    /// @return The four-center integral.
    I3CIntegral _get_integral(const std::string&        label,
                              const std::array<int, 3>& ang_moms) const;
    
    /// Generates set of integrals required for horizontal Obara-Saika recursion on ket side.
    /// @param integral The base four center integral.
    /// @return The set of integrals.
    SI3CIntegrals _generate_ket_hrr_integral_group(const I3CIntegral& integral) const;
    
    /// Generates set of integrals required for vertical Obara-Saika recursion.
    /// @param integral The base four center integral.
    /// @param integrals The set of four center integrals.
    /// @return The set of integrals.
    SI3CIntegrals _generate_vrr_integral_group(const I3CIntegral&   integral,
                                               const SI3CIntegrals& integrals) const;
    
    /// Writes header file for recursion.
    /// @param hrr_integrals The set of unique integrals for horizontal recursion.
    /// @param vrr_integrals The set of unique integrals for vertical recursion.
    /// @param integral The base two center integral.
    void _write_cpp_header(const SI3CIntegrals& hrr_integrals,
                           const SI3CIntegrals& vrr_integrals,
                           const I3CIntegral& integral) const;
    
    /// Gets file name of file with recursion functions for two center integral.
    /// @param integral The base two center integral.
    /// @return The file name.
    std::string _file_name(const I3CIntegral& integral) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_defines(      std::ofstream& fstream,
                            const I3CIntegral&   integral,
                            const bool           start) const;
    
    /// Writes definitions of includes for header file.
    /// @param fstream the file stream.
    /// @param hrr_integrals The set of unique integrals for horizontal recursion.
    /// @param vrr_integrals The set of unique integrals for vertical recursion.
    /// @param integral The base two center integral.
    void _write_hpp_includes(      std::ofstream& fstream,
                             const SI3CIntegrals& hrr_integrals,
                             const SI3CIntegrals& vrr_integrals,
                             const I3CIntegral&   integral) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const I3CIntegral&   integral,
                          const bool           start) const;
    
    /// Writes primitive header file for recursion.
    /// @param integral The base two center integral.
    void _write_prim_cpp_header(const I3CIntegral& integral) const;
    
    /// Writes C++ code file for primtive recursion.
    /// @param integral The base four center integral.
    void _write_prim_cpp_file(const I3CIntegral& integral) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of define (start or end).
    void _write_prim_hpp_defines(      std::ofstream& fstream,
                                 const I3CIntegral&   integral,
                                 const bool           start) const;
    
    /// Writes definitions of includes for primitive header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_prim_hpp_includes(      std::ofstream& fstream,
                                  const I3CIntegral&   integral) const;
    
    /// Writes definitions of includes for primitive header file.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void _write_prim_cpp_includes(      std::ofstream& fstream,
                                  const I3CIntegral&  integral) const;
    
    /// Writes ket hrr header file for recursion.
    /// @param integral The base two center integral.
    void _write_hrr_cpp_header(const I3CIntegral& integral) const;
    
    /// Writes C++ code file for primtive recursion.
    /// @param integral The base four center integral.
    void _write_hrr_cpp_file(const I3CIntegral& integral) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hrr_hpp_defines(      std::ofstream& fstream,
                                const I3CIntegral&   integral,
                                const bool           start) const;
    
    /// Writes definitions of includes for primitive header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_hrr_hpp_includes(      std::ofstream& fstream,
                                 const I3CIntegral&   integral) const;
    
    /// Writes definitions of includes for primitive header file.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void _write_hrr_cpp_includes(      std::ofstream& fstream,
                                 const I3CIntegral&   integral) const;
    
public:
    /// Creates a three-center integrals CPU code generator.
    T3CCPUGenerator() = default;
     
    /// Generates selected three-center integrals up to given angular momentum (inclusive)  on A, B, C, and D centers.
    /// @param label The label of requested two-center integral.
    /// @param max_ang_mom The maximum angular momentum of C and D centers.
    /// @param max_aux_ang_mom The maximum angular momentum of auxilary A center.
    void generate(const std::string& label,
                  const int          max_ang_mom,
                  const int          max_aux_ang_mom) const;
};

#endif /* t3c_cpu_generators_hpp */
