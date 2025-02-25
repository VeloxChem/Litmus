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

#ifndef t3c_geom_hrr_cpu_generators_hpp
#define t3c_geom_hrr_cpu_generators_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <utility>

#include "t3c_defs.hpp"

// Geometrical derivatives of three-center integrals code generator for CPU.
class T3CGeomHrrCPUGenerator
{
    /// Checks if recursion is available for four-center inetgral with given label.
    /// @param label The label of requested four-center integral.
    bool _is_available(const std::string& label) const;
    
    /// Gets four-center inetgral with requested label.
    /// @param label The label of requested four-center integral.
    /// @param ang_moms The angular momentum of  A, B, C, and D centers.
    /// @param geom_drvs The geometrical derivative of bra side, integrand, and  ket side.
    /// @return The four-center integral.
    I3CIntegral _get_integral(const std::string&        label,
                              const std::array<int, 3>& ang_moms,
                              const std::array<int, 3>& geom_drvs) const;
    
    /// Writes ket hrr header file for recursion.
    /// @param integral The base two center integral.
    void _write_bra_hrr_cpp_header(const I3CIntegral& integral) const;
    
    /// Writes C++ code file for primtive recursion.
    /// @param integral The base four center integral.
    void _write_bra_hrr_cpp_file(const I3CIntegral& integral) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of define (start or end).
    void _write_bra_hrr_hpp_defines(      std::ofstream& fstream,
                                    const I3CIntegral&   integral,
                                    const bool           start) const;
        
    /// Writes definitions of includes for primitive header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_bra_hrr_hpp_includes(      std::ofstream& fstream,
                                     const I3CIntegral&   integral) const;
    
    /// Writes definitions of includes for primitive header file.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void _write_bra_hrr_cpp_includes(      std::ofstream& fstream,
                                     const I3CIntegral&  integral) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const I3CIntegral&   integral,
                          const bool           start) const;

    /// Writes ket hrr header file for recursion.
    /// @param integral The base two center integral.
    void _write_ket_hrr_cpp_header(const I3CIntegral& integral) const;
    
    /// Writes C++ code file for primtive recursion.
    /// @param integral The base four center integral.
    void _write_ket_hrr_cpp_file(const I3CIntegral& integral) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of define (start or end).
    void _write_ket_hrr_hpp_defines(      std::ofstream& fstream,
                                    const I3CIntegral&   integral,
                                    const bool           start) const;
        
    /// Writes definitions of includes for primitive header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void _write_ket_hrr_hpp_includes(      std::ofstream& fstream,
                                     const I3CIntegral&   integral) const;
    
    /// Writes definitions of includes for primitive header file.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void _write_ket_hrr_cpp_includes(      std::ofstream& fstream,
                                     const I3CIntegral&  integral) const;
    
public:
    /// Creates a geometrical derivatives of three-center integrals CPU code generator.
    T3CGeomHrrCPUGenerator() = default;
     
    /// Generates selected three-center integrals up to given angular momentum (inclusive)  on A, B, C, and D centers.
    /// @param label The label of requested two-center integral.
    /// @param max_ang_mom The maximum angular momentum of A, C and D centers.
    /// @param geom_drvs The geometrical derivative of bra side, integrand, and  ket side.
    void generate(const std::string&        label,
                  const int                 max_ang_mom,
                  const std::array<int, 3>& geom_drvs) const;
};

#endif /* t3c_geom_hrr_cpu_generators_hpp */
