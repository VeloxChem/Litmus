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

#ifndef v4i_geom20_eri_driver_hpp
#define v4i_geom20_eri_driver_hpp

#include <array>

#include "t4c_defs.hpp"

/// Four center electron repulsion integrals driver class.
class V4IGeom20ElectronRepulsionDriver
{
   
public:
    /// Creates a four center electron repulsion integrals driver.
    V4IGeom20ElectronRepulsionDriver() = default;
    
    /// Check if integral is for four-center electron repulsion integral.
    /// @param integral The integral to check.
    /// @return True if reccursion expansion belongs to electron repulsion recursion, False otherwise.
    bool is_electron_repulsion(const I4CIntegral& integral) const;
    
    /// Applies horizontal recursion to bra side of electron repulsion integral.
    /// @param integral The  electron repulsion integral.
    /// @return The set of integrals.
    SI4CIntegrals bra_hrr(const I4CIntegral& integral) const;
    
    /// Applies bra hrr recursion expansion for set of integral.
    /// @param integral The  integral to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals apply_bra_hrr_recursion(const I4CIntegral& integral) const;
    
    /// Creates bra hrr recursion expansion for set of integral.
    /// @param integrals The  set of integrals to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals create_bra_hrr_recursion(const SI4CIntegrals& integrals) const;
};

#endif /* v4i_geom20_eri_driver_hpp */
