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

#ifndef t2c_linmom_driver_hpp
#define t2c_linmom_driver_hpp

#include <optional>
#include <array>

#include "tensor_component.hpp"
#include "t2c_defs.hpp"

/// Two center multipole integrals driver class.
class T2CLinearMomentumDriver
{

public:
    /// Creates a two center overlap integrals driver.
    T2CLinearMomentumDriver() = default;

    /// Check if recursion term is for two-center multipole integral.
    /// @param rterm The recursion term.
    /// @return True if reccursion expansion belongs to multipole recursion, False otherwise.
    bool is_linear_momentum(const R2CTerm& rterm) const;

    /// Applies vertical recursion to bra side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R2CDist> op_vrr(const R2CTerm& rterm) const;

    R2CDist apply_op_vrr(const R2CTerm& rterm) const;

};

#endif /* t2c_linmom_driver_hpp */

