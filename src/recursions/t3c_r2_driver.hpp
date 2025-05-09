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

#ifndef t3c_r2_driver_hpp
#define t3c_r2_driver_hpp

#include <optional>
#include <array>

#include "tensor_component.hpp"
#include "t2c_defs.hpp"

/// Three center overlap gradient integrals driver class.
class T3CR2Driver
{
    /// Cartesian coordinate tensor components.
    std::array<TensorComponent, 3> _rxyz;
    
public:
    /// Creates a three center overlap gradient integrals driver.
    T3CR2Driver();
    
    /// Check if recursion term is for three-center overlap integral.
    /// @param rterm The recursion term.
    /// @return True if reccursion expansion belongs to overlap recursion, False otherwise.
    bool is_r2(const R2CTerm& rterm) const;
    
    /// Applies vertical recursion to auxilary side of given recursion term.
    /// @param rterm The recursion term.
    /// @return The recursion expansion of given recursion term.
    R2CDist aux_vrr(const R2CTerm& rterm) const;
};

#endif /* t3c_r2_driver_hpp */
