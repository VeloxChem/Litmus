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

#ifndef v3i_geom010_eri_driver_hpp
#define v3i_geom010_eri_driver_hpp

#include <array>

#include "t3c_defs.hpp"

/// Three center electron repulsion integrals driver class.
class V3IGeom010ElectronRepulsionDriver
{
   
public:
    /// Creates a three center electron repulsion integrals driver.
    V3IGeom010ElectronRepulsionDriver() = default;
    
    /// Check if integral is for three-center electron repulsion integral.
    /// @param integral The integral to check.
    /// @return True if reccursion expansion belongs to electron repulsion recursion, False otherwise.
    bool is_electron_repulsion(const I3CIntegral& integral) const;
    
    /// Applies horizontal recursion to bra side of electron repulsion integral.
    /// @param integral The  electron repulsion integral.
    /// @return The set of integrals.
    SI3CIntegrals ket_hrr(const I3CIntegral& integral) const;
    
    /// Applies horizontal recursion to bra side of electron repulsion integral.
    /// @param integral The  electron repulsion integral.
    /// @return The set of integrals.
    SI3CIntegrals ket_aux_hrr(const I3CIntegral& integral) const;
    
    /// Applies bra hrr recursion expansion for set of integral.
    /// @param integral The  integral to apply recursion.
    /// @return The set of integrals.
    SI3CIntegrals apply_ket_hrr_recursion(const I3CIntegral& integral) const;
};

#endif /* v3i_geom010_eri_driver_hpp */
