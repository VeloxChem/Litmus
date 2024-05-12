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

#ifndef v4i_center_driver_hpp
#define v4i_center_driver_hpp

#include <optional>
#include <array>

#include "t4c_defs.hpp"

/// Geometrical derivatives four center integrals driver class.
class V4ICenterDriver
{

public:
    /// Creates a geometrical derivatives four center integrals driver.
    V4ICenterDriver() = default;

    /// Check if integral is  auxilary geometrical derivatives four center integral.
    /// @param integral The integral to check.
    /// @param index The index of center to check for auxilary value.
    /// @return True if integral is auxilary, False otherwise.
    bool is_auxiliary(const I4CIntegral& integral, const int index) const;

    /// Applies vertical recursion to bra side of geometrical derivatives four center.
    /// @param integral The  geometrical derivatives four center.
    /// @param index The index of center targeted by recursion.
    /// @return The set of integrals.
    SI4CIntegrals bra_ket_vrr(const I4CIntegral& integral, const int index) const;
    
    /// Applies vertical recursion to bra side of overlap integral.
    /// @param integral The  overlap integral.
    /// @return The recursion expansion of integral.
    SI4CIntegrals apply_bra_ket_vrr(const I4CIntegral& integral) const;
};

#endif /* v4i_center_driver_hpp */
