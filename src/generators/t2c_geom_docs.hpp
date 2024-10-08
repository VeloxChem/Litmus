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

#ifndef t2c_geom_docs_hpp
#define t2c_geom_docs_hpp

#include <string>
#include <vector>
#include <fstream>

#include "t2c_defs.hpp"

// Two-center documentation generator for CPU.
class T2CGeomDocuDriver
{
    /// Generates compute string.
    /// @param integral The base two center integral.
    /// @param geom_drvs The geometrical derivative of bra and  ket sides.
    /// @return The compute string.
    std::string _get_compute_str(const I2CIntegral&        integral,
                                 const std::array<int, 3>& geom_drvs) const;
    
    /// Generates vector of buffer strings.
    /// @param geom_integrals The set of unique integrals for geometrical recursion.
    /// @param integral The base two center integral.
    /// @param geom_drvs The geometrical derivative of bra and  ket sides.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_buffers_str(const SI2CIntegrals&      geom_integrals,
                                              const I2CIntegral&        integral,
                                              const std::array<int, 3>& geom_drvs) const;
    
    /// Generates vector of recursion variables strings.
    /// @param integral The base two center integral.
    /// @return The vector of recursion variables strings.
    std::vector<std::string> _get_recursion_variables_str(const I2CIntegral& integral) const;
    
public:
    /// Creates a primtive four-center documentation generator.
    T2CGeomDocuDriver() = default;
    
    /// Writes documentation string for primtive compute function.
    /// @param fstream the file stream.
    /// @param geom_integrals The set of unique integrals for geometrical recursion.
    /// @param integral The base four center integral.
    /// @param geom_drvs The geometrical derivative of bra and  ket sides.
    void write_doc_str(      std::ofstream&      fstream,
                       const SI2CIntegrals&      geom_integrals,
                       const I2CIntegral&        integral,
                       const std::array<int, 3>& geom_drvs) const;
};


#endif /* t2c_geom_docs_hpp */
