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

#ifndef v4i_eri_driver_hpp
#define v4i_eri_driver_hpp

#include <array>

#include "t4c_defs.hpp"

/// Four center electron repulsion integrals driver class.
class V4IElectronRepulsionDriver
{
   
public:
    /// Creates a four center electron repulsion integrals driver.
    V4IElectronRepulsionDriver() = default;
    
    /// Check if integral is for four-center electron repulsion integral.
    /// @param integral The integral to check.
    /// @return True if reccursion expansion belongs to electron repulsion recursion, False otherwise.
    bool is_electron_repulsion(const I4CIntegral& integral) const;
    
    /// Applies horizontal recursion to bra side of electron repulsion integral.
    /// @param integral The  electron repulsion integral.
    /// @return The set of integrals.
    SI4CIntegrals bra_hrr(const I4CIntegral& integral) const;
    
    /// Applies horizontal recursion to ket side of electron repulsion integral.
    /// @param integral The  electron repulsion integral.
    /// @return The set of integrals.
    SI4CIntegrals ket_hrr(const I4CIntegral& integral) const;
    
    /// Applies vertical recursion to bra side center A of given recursion term.
    /// @param integral The  electron repulsion integral.
    /// @return The set of integrals.
    SI4CIntegrals bra_vrr_a(const I4CIntegral& integral) const;
    
    /// Applies vertical recursion to bra side center B of given recursion term.
    /// @param integral The  electron repulsion integral.
    /// @return The set of integrals.
    SI4CIntegrals bra_vrr_b(const I4CIntegral& integral) const;
    
    /// Applies vertical recursion to ket side center C of given recursion term.
    /// @param integral The  electron repulsion integral.
    /// @return The set of integrals.
    SI4CIntegrals ket_vrr_c(const I4CIntegral& integral) const;
    
    /// Applies vertical recursion to ket side center D of given recursion term.
    /// @param integral The  electron repulsion integral.
    /// @return The set of integrals.
    SI4CIntegrals ket_vrr_d(const I4CIntegral& integral) const;
    
    /// Applies vertical recursion to bra side of electron repulsion integral.
    /// @param integral The  electron repulsion integral.
    /// @return The set of integrals.
    SI4CIntegrals bra_vrr(const I4CIntegral& integral) const;
    
    /// Applies vertical recursion to ket side of electron repulsion integral.
    /// @param integral The  electron repulsion integral.
    /// @return The set of integrals.
    SI4CIntegrals ket_vrr(const I4CIntegral& integral) const;
    
    /// Applies bra hrr recursion expansion for set of integral.
    /// @param integral The  integral to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals apply_bra_hrr_recursion(const I4CIntegral& integral) const;
    
    /// Applies ket hrr recursion expansion for set of integral.
    /// @param integral The  integral to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals apply_ket_hrr_recursion(const I4CIntegral& integral) const;
    
    /// Applies bra vrr recursion expansion for set of integral.
    /// @param integral The  integral to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals apply_bra_vrr_recursion(const I4CIntegral& integral) const;
    
    /// Applies ket vrr recursion expansion for set of integral.
    /// @param integral The  integral to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals apply_ket_vrr_recursion(const I4CIntegral& integral) const;
    
    /// Creates bra vrr recursion on center A.
    /// @param integrals The  set of integrals to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals apply_bra_vrr_a(const SI4CIntegrals& integrals) const;
    
    /// Creates bra vrr recursion on center B.
    /// @param integrals The  set of integrals to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals apply_bra_vrr_b(const SI4CIntegrals& integrals) const;
    
    /// Creates ket vrr recursion on center C.
    /// @param integrals The  set of integrals to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals apply_ket_vrr_c(const SI4CIntegrals& integrals) const;
    
    /// Creates bra vrr recursion on center D.
    /// @param integrals The  set of integrals to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals apply_ket_vrr_d(const SI4CIntegrals& integrals) const;
    
    /// Creates bra hrr recursion expansion for set of integral.
    /// @param integrals The  set of integrals to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals create_bra_hrr_recursion(const SI4CIntegrals& integrals) const;
    
    /// Creates ket hrr recursion expansion for set of integral.
    /// @param integrals The  set of integrals to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals create_ket_hrr_recursion(const SI4CIntegrals& integrals) const;
    
    /// Creates vrr recursion expansion for set of integral.
    /// @param integrals The  set of integrals to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals create_vrr_recursion(const SI4CIntegrals& integrals) const;
    
    /// Creates full vrr recursion expansion for set of integral.
    /// @param integrals The  set of integrals to apply recursion.
    /// @return The set of integrals.
    SI4CIntegrals create_full_vrr_recursion(const SI4CIntegrals& integrals) const;
};


#endif /* v4i_eri_driver_hpp */
