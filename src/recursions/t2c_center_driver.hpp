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

#ifndef t2c_center_driver_hpp
#define t2c_center_driver_hpp

#include <optional>
#include <array>

#include "tensor_component.hpp"
#include "t2c_defs.hpp"

/// Two center geometrical prefix operator driver class.
class T2CCenterDriver
{
    /// Cartesian coordinate tensor components.
    std::array<TensorComponent, 3> _rxyz;
    
public:
    /// Creates a two center geometrical prefix operator driver.
    T2CCenterDriver();
    
    /// Check if recursion term is auxilary with respect of geometrical prefix operator.
    /// @param rterm The recursion term.
    /// @param index The index of prefix operator.
    /// @return True if reccursion term is auxilary, False otherwise.
    bool is_auxilary(const R2CTerm& rterm,
                     const int      index) const;
    
    /// Applies vertical recursion to geometrical prefix operator acting bra or ket side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @param index The index of prefix operator.
    /// @return The recursion expansion of given recursion term.
    std::optional<R2CDist> bra_ket_vrr(const R2CTerm& rterm,
                                       const char     axis,
                                       const int      index) const;
    
    /// Applies vertical recursion to geometrical prefix operator acting bra or ket side of given recursion term.
    /// @param rterm The recursion term with overlap integral.
    /// @param index The index of prefix operator.
    /// @return The recursion expansion of given recursion term.
    R2CDist apply_bra_ket_vrr(const R2CTerm& rterm,
                              const int      index) const;
    
    
    /// Recursively applies Obara-Saika recursion to recursion expansion.
    /// @param rdist The recursion expansion.
    void apply_recursion(R2CDist& rdist) const;
    
    /// Recursively applies vertical recursion to geometrical prefix operator acting bra or ket side of given recursion expansion.
    /// @param rdist The recursion expansion.
    /// @param index The index of prefix operator.
    void apply_bra_ket_vrr(      R2CDist& rdist,
                           const int      index) const;
    
    /// Creates recursion group from vector of integral components.
    /// @param vints The  vector of integral components.
    /// @return The recursion group.
    R2Group create_recursion(const VT2CIntegrals& vints) const;
};


#endif /* t2c_center_driver_hpp */
