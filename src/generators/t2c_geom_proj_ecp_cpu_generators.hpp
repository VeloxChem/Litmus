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

#ifndef t2c_geom_proj_ecp_cpu_generators_hpp
#define t2c_geom_proj_ecp_cpu_generators_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <utility>

#include "t2c_defs.hpp"
#include "t2c_prim_docs.hpp"
#include "file_stream.hpp"

// Two-center integrals code generator for CPU.
class T2CGeomProjECPCPUGenerator
{
    /// Checks if recursion is available for two-center inetgral with given label.
    /// @param label The label of requested two-center integral.
    bool _is_available(const std::string& label) const;
    
    /// Gets two-center inetgral with requested label.
    /// @param label The label of requested two-center integral.
    /// @param ang_moms The angular momentum of  A and B centers.
    /// @param geom_drvs The geometrical derivative of bra and  ket sides.
    /// @return The two-center integral.
    M2Integral _get_integral(const std::string&        label,
                             const std::array<int, 2>& ang_moms,
                             const int                 proj_ang_mom,
                             const std::array<int, 3>& geom_drvs) const;
    
    /// Generates set of integrals required for geometrical derivatives.
    /// @param integral The base four center integral.
    /// @return The set of integrals.
    SM2Integrals _generate_geom_integral_group(const M2Integral& integral) const;
    
    /// Generates set of integrals required for geometrical derivatives.
    /// @param integral The base four center integral.
    /// @param integrals The set of geometrical derivative integrals.
    /// @return The set of integrals.
    SM2Integrals _generate_vrr_integral_group(const M2Integral&   integral,
                                              const SM2Integrals& integrals) const;
    
    /// Writes header file for recursion.
    /// @param geom_integrals The set of unique integrals for geometrical recursion.
    /// @param vrr_integrals The set of unique integrals for vertical recursion.
    /// @param integral The base two center integral.
    /// @param geom_drvs The geometrical derivative of bra side, integrand, and  ket side.
    void _write_cpp_header(const SM2Integrals&       geom_integrals,
                           const SM2Integrals&       vrr_integrals,
                           const M2Integral&         integral,
                           const std::array<int, 3>& geom_drvs) const;
    
    /// Gets file name of file with recursion functions for two center integral.
    /// @param integral The base two center integral.
    /// @return The file name.
    std::string _file_name(const M2Integral& integral) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_defines(      std::ofstream& fstream,
                            const M2Integral&    integral,
                            const bool           start) const;
    
    /// Writes definitions of includes for header file.
    /// @param fstream the file stream.
    /// @param integrals The set of unique integrals.
    /// @param integral The base two center integral.
    void _write_hpp_includes(      std::ofstream&      fstream,
                             const SM2Integrals&       integrals,
                             const M2Integral&         integral,
                             const std::array<int, 3>& geom_drvs) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const M2Integral&    integral,
                          const bool           start) const;
    
public:
    /// Creates a two-center integrals CPU code generator.
    T2CGeomProjECPCPUGenerator() = default;
     
    /// Generates selected two-center integrals up to given angular momentum (inclusive)  on A and B centers.
    /// @param label The label of requested two-center integral.
    /// @param max_ang_mom The maximum angular momentum of A and B centers.
    /// @param proj_ang_mom The maximum angular momentum of projector on center C.
    /// @param geom_drvs The geometrical derivative of bra and  ket sides.
    void generate(const std::string&        label,
                  const int                 max_ang_mom,
                  const int                 proj_ang_mom,
                  const std::array<int, 3>& geom_drvs) const;
};


#endif /* t2c_geom_proj_ecp_cpu_generators_hpp */
