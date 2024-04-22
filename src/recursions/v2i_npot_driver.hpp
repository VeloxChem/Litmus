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

#ifndef v2i_npot_driver_hpp
#define v2i_npot_driver_hpp

#include <array>

#include "t2c_defs.hpp"

/// Two center nuclear potential integrals driver class.
class V2INuclearPotentialDriver
{
   
public:
    /// Creates a two center nuclear potential integrals driver.
    V2INuclearPotentialDriver() = default;
    
    /// Check if integral is for two-center nuclear potential integral.
    /// @param integral The integral to check.
    /// @return True if reccursion expansion belongs to overlap recursion, False otherwise.
    bool is_nuclear_potential(const I2CIntegral& integral) const;
    
    /// Applies vertical recursion to bra side of nuclear potential integral.
    /// @param integral The  overlap integral.
    /// @return The set of integrals.
    SI2CIntegrals bra_vrr(const I2CIntegral& integral) const;
    
    /// Applies vertical recursion to ket side of nuclear potential integral.
    /// @param integral The  overlap integral.
    /// @return The set of integrals.
    SI2CIntegrals ket_vrr(const I2CIntegral& integral) const;
    
    /// Applies vertical recursion to auxilary nuclear potential integral.
    /// @param integral The  overlap integral.
    /// @return The integral.
    I2CIntegral aux_vrr(const I2CIntegral& integral) const;
    
    /// Applies vertical recursion to bra side of nuclear potential integral.
    /// @param integral The  overlap integral.
    /// @return The recursion expansion of integral.
    SI2CIntegrals apply_bra_vrr(const I2CIntegral& integral) const;
    
    /// Applies vertical recursion to ket side of nuclear potential integral.
    /// @param integral The  overlap integral.
    /// @return The recursion expansion of integral.
    SI2CIntegrals apply_ket_vrr(const I2CIntegral& integral) const;
    
    /// Recursively applies Obara-Saika recursion to recursion expansion.
    /// @param integrals The  integral to apply recursion.
    /// @return The set of integrals.
    SI2CIntegrals apply_recursion(const SI2CIntegrals& integrals) const;
    
    /// Creates recursion expansion for set of integral.
    /// @param integrals The  integral to apply recursion.
    /// @return The set of integrals.
    SI2CIntegrals create_recursion(const SI2CIntegrals& integrals) const;
};


#endif /* v2i_npot_driver_hpp */
