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

#ifndef t2c_ecp_cpu_generators_hpp
#define t2c_ecp_cpu_generators_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <utility>

#include "t2c_defs.hpp"

// Two-center integrals code generator for CPU.
class T2CECPCPUGenerator
{
    /// Checks if recursion is available for two-center inetgral with given label.
    /// @param label The label of requested two-center integral.
    bool _is_available(const std::string& label) const;
    
    /// Gets two-center inetgral with requested label.
    /// @param label The label of requested two-center integral.
    /// @param ang_moms The angular momentum of  A and B centers.
    /// @return The two-center integral.
    I2CIntegral _get_integral(const std::string&        label,
                              const std::array<int, 2>& ang_moms) const;
    
    /// Generates set of integrals required for vertical Obara-Saika recursion.
    /// @param integral The base two center integral.
    /// @return The set of integrals.
    SI2CIntegrals _generate_integral_group(const I2CIntegral& integral) const;
    
    /// Filters integrals required for HRR recursion.
    /// @param integrals The set of integrals.
    /// @return The set of integrals required by HRR recursion.
    SI2CIntegrals _filter_hrr_integrals(const SI2CIntegrals& integrals) const;
    
    /// Filters integrals required for VRR recursion.
    /// @param integrals The set of integrals.
    /// @return The set of integrals required by VRR recursion.
    SI2CIntegrals _filter_vrr_integrals(const SI2CIntegrals& integrals) const;
    
    /// Writes header file for recursion.
    /// @param hrr_integrals The set of unique HRR integrals.
    /// @param vrr_integrals The set of unique VRR integrals.
    /// @param integral The base two center integral.
    void _write_cpp_header(const SI2CIntegrals& hrr_integrals,
                           const SI2CIntegrals& vrr_integrals,
                           const I2CIntegral&   integral) const;
    
    /// Gets file name of file with recursion functions for two center integral.
    /// @param integral The base two center integral.
    /// @return The file name.
    std::string _file_name(const I2CIntegral& integral) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param is_prim_rec The flag to indicate primitive recurion.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_defines(      std::ofstream& fstream,
                            const I2CIntegral&   integral,
                            const bool           is_prim_rec,
                            const bool           start) const;
    
    /// Writes definitions of includes for header file.
    /// @param fstream the file stream.
    /// @param hrr_integrals The set of unique HRR integrals.
    /// @param vrr_integrals The set of unique VRR integrals.
    /// @param integral The base two center integral.
    void _write_hpp_includes(      std::ofstream& fstream,
                             const SI2CIntegrals& hrr_integrals,
                             const SI2CIntegrals& vrr_integrals,
                             const I2CIntegral&   integral) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const I2CIntegral&   integral,
                          const bool           start) const;
    
public:
    /// Creates a two-center integrals CPU code generator.
    T2CECPCPUGenerator() = default;
     
    /// Generates selected two-center integrals up to given angular momentum (inclusive)  on A and B centers.
    /// @param label The label of requested two-center integral.
    /// @param max_ang_mom The maximum angular momentum of A and B centers.
    void generate(const std::string& label,
                  const int          max_ang_mom) const;
};

#endif /* t2c_ecp_cpu_generators_hpp */
