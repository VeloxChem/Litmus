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

#ifndef geom_center_driver_hpp
#define geom_center_driver_hpp

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
};


#endif /* geom_center_driver_hpp */
