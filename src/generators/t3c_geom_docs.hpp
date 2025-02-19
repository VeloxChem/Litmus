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

#ifndef t3c_geom_docs_hpp
#define t3c_geom_docs_hpp

#include <string>
#include <vector>
#include <fstream>

#include "t3c_defs.hpp"

// Three-center documentation generator for CPU.
class T3CGeomDocuDriver
{
    /// Generates compute string.
    /// @param integral The base two center integral.
    /// @return The compute string.
    std::string _get_compute_str(const I3CIntegral& integral) const;
    
    /// Generates vector of buffer strings.
    /// @param geom_integrals The set of unique integrals for geometrical recursion.
    /// @param integral The base two center integral.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_buffers_str(const SI3CIntegrals& geom_integrals,
                                              const I3CIntegral& integral) const;
    
    /// Generates vector of coordinates strings.
    /// @param integral The base two center integral.
    /// @return The vector of coordinates strings.
    std::vector<std::string> _get_coordinates_str(const I3CIntegral& integral) const;
    
    /// Generates vector of recursion variables strings.
    /// @param integral The base two center integral.
    /// @return The vector of recursion variables strings.
    std::vector<std::string> _get_recursion_variables_str(const I3CIntegral& integral) const;
    
    /// Generates vector of matrix strings.
    /// @param integral The base two center integral.
    /// @return The vector of matrix strings.
    std::vector<std::string> _get_matrices_str(const I3CIntegral& integral) const;
    
    /// Generates vector of GTOs block strings.
    /// @param integral The base two center integral.
    /// @return The vector of GTOs block strings,
    std::vector<std::string> _get_gto_pair_blocks_str(const I3CIntegral& integral) const;
  
    /// Generates vector of indices strings.
    /// @return The vector of indices strings.
    std::vector<std::string> _get_indices_str() const;
    
    /// Generates compute string.
    /// @param integral The base two center integral.
    /// @return The compute string.
    std::string _get_bra_geom_compute_str(const I3CIntegral& integral) const;
    
    /// Generates vector of buffer strings.
    /// @param integral The base two center integral.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_bra_geom_buffers_str(const I3CIntegral& integral) const;
    
    /// Generates vector of coordinates strings.
    /// @param integral The base two center integral.
    /// @return The vector of coordinates strings.
    std::vector<std::string> _get_bra_coordinates_str(const I3CIntegral& integral) const;
    
    /// Generates vector of recursion variables strings.
    /// @param integral The base two center integral.
    /// @return The vector of recursion variables strings.
    std::vector<std::string> _get_bra_recursion_variables_str(const I3CIntegral& integral) const;
    
    /// Generates compute string.
    /// @param integral The base two center integral.
    /// @return The compute string.
    std::string _get_ket_geom_compute_str(const I3CIntegral& integral) const;
    
    /// Generates vector of buffer strings.
    /// @param integral The base two center integral.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_ket_geom_buffers_str(const I3CIntegral& integral) const;
    
    /// Generates vector of coordinates strings.
    /// @param integral The base two center integral.
    /// @return The vector of coordinates strings.
    std::vector<std::string> _get_ket_coordinates_str(const I3CIntegral& integral) const;
    
    /// Generates vector of recursion variables strings.
    /// @param integral The base two center integral.
    /// @return The vector of recursion variables strings.
    std::vector<std::string> _get_ket_recursion_variables_str(const I3CIntegral& integral) const;
    
public:
    /// Creates a primtive three-center documentation generator.
    T3CGeomDocuDriver() = default;
    
    /// Writes documentation string for primtive compute function.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void write_doc_str(      std::ofstream& fstream,
                       const I3CIntegral&   integral) const;
    
    /// Writes documentation string for primtive compute function.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void write_bra_geom_doc_str(      std::ofstream& fstream,
                                const I3CIntegral&   integral) const;
    
    /// Writes documentation string for primtive compute function.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void write_ket_geom_doc_str(      std::ofstream& fstream,
                                const I3CIntegral&   integral) const;
    
};

#endif /* t3c_geom_docs_hpp */
