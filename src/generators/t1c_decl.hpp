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

#ifndef t1c_decl_hpp
#define t1c_decl_hpp

#include <string>
#include <vector>
#include <fstream>
#include <utility>

#include "t2c_defs.hpp"

// GTOs functions declaration generator for CPU.
class T1CDeclDriver
{
    /// Generates vector of gto strings.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    /// @return The vector of gto strings.
    std::vector<std::string> _get_gto_str(const int angmom,
                                          const int gdrv) const;
    /// Generates vector of variable  strings.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    /// @param terminus The flag to add termination symbol.
    /// @return The vector of index strings.
    std::vector<std::string> _get_vars_str(const int  angmom,
                                           const int  gdrv,
                                           const bool terminus) const;
    
    /// Generates compute function name.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    /// @return The compute function name,
    std::string _func_name(const int angmom,
                           const int gdrv) const;

public:
    /// Creates a GTO functions declaration generator.
    T1CDeclDriver() = default;
    
    /// Writes declaration for compute function.
    /// @param fstream the file stream.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    /// @param terminus The flag to add termination symbol.
    void write_func_decl(      std::ofstream& fstream,
                         const int            angmom,
                         const int            gdrv,
                         const bool           terminus) const;
};


#endif /* t1c_decl_hpp */
