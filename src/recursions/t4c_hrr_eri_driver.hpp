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

#ifndef t4c_hrr_eri_driver_hpp
#define t4c_hrr_eri_driver_hpp

#include <optional>
#include <array>

#include "tensor_component.hpp"
#include "t4c_defs.hpp"

/// Four center horizontal recursion electron repulsion integrals driver class.
class T4CHrrElectronRepulsionDriver
{
    /// Cartesian coordinate tensor components.
    std::array<TensorComponent, 3> _rxyz;
    
public:
    /// Creates a four center electron repulsion integrals driver.
    T4CHrrElectronRepulsionDriver();
    
    /// Check if recursion term is for four-center electron repulsion integral.
    /// @param rterm The recursion term.
    /// @return True if reccursion expansion belongs to electron repulsion recursion, False otherwise.
    bool is_electron_repulsion(const R4CTerm& rterm) const;
    
    /// Applies horizontal recursion to bra side  of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R4CDist> bra_hrr(const R4CTerm& rterm,
                                   const char     axis) const;
    
    /// Applies horizontal recursion to ket side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R4CDist> ket_hrr(const R4CTerm& rterm,
                                   const char     axis) const;
    
    /// Applies horizontal recursion to bra side recursion term containing
    /// electron repulsion integral.
    /// @param rterm The recursion term with electron repulsion integral.
    /// @return The recursion expansion of given recursion term.
    R4CDist apply_bra_hrr(const R4CTerm& rterm) const;
    
    /// Applies horizontal recursion to ket side  recursion term containing
    /// electron repulsion integral.
    /// @param rterm The recursion term with electron repulsion integral.
    /// @return The recursion expansion of given recursion term.
    R4CDist apply_ket_hrr(const R4CTerm& rterm) const;
    
    /// Recursively applies Obara-Saika recursion to recursion expansion.
    /// @param rdist The recursion expansion.
    void apply_recursion(R4CDist& rdist) const;
    
    /// Recursively applies horizontal recursion to bra side center A of given recursion expansion.
    /// @param rdist The recursion expansion.
    void apply_bra_hrr(R4CDist& rdist) const;
    
    /// Recursively applies horizontal recursion to ket side center C of given recursion expansion.
    /// @param rdist The recursion expansion.
    void apply_ket_hrr(R4CDist& rdist) const;
    
    /// Creates recursion group from vector of electron repulsion integral components.
    /// @param vints The  vector of electron repulsion integral components.
    /// @return The recursion group.
    R4Group create_recursion(const VT4CIntegrals& vints) const;
    
    /// Recursively applies Obara-Saika recursion to recursion group.
    /// @param rgroup The recursion group.
    void apply_recursion(R4Group& rgroup) const;
};

#endif /* t4c_hrr_eri_driver_hpp */
