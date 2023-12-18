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

#ifndef t1c_docs_hpp
#define t1c_docs_hpp

#include <string>
#include <vector>
#include <fstream>

#include "t2c_defs.hpp"

// GTOs documentation generator for CPU.
class T1CDocuDriver
{
    /// Generates vector of variable strings.
    /// @return The vector of special variable strings.
    std::vector<std::string> _get_vars_str() const;
    
    /// Generates compute string.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    /// @return The  compute string.
    std::string _get_compute_str(const int angmom,
                                 const int gdrv) const;
    
public:
    /// Creates a GTOs documentation generator.
    T1CDocuDriver() = default;
    
    /// Writes documentation string for compute function.
    /// @param fstream the file stream.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    void write_doc_str(      std::ofstream& fstream,
                       const int            angmom,
                       const int            gdrv) const;
};


#endif /* t1c_docs_hpp */
