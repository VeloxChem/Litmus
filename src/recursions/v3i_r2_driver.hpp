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

#ifndef v3i_r2_driver_hpp
#define v3i_r2_driver_hpp

#include <array>

#include "t2c_defs.hpp"

/// Three center overlap integrals driver class.
class V3IR2Driver
{
   
public:
    /// Creates a three center r2 integrals driver.
    V3IR2Driver() = default;
    
    /// Check if integral is for three-center r2 integral.
    /// @param integral The integral to check.
    /// @return True if reccursion expansion belongs to r2 recursion, False otherwise.
    bool is_r2(const I2CIntegral& integral) const;
    
    /// Applies vertical recursion to auxilary three center overlap integral.
    /// @param integral The  overlap integral.
    /// @return The integral.
    SI2CIntegrals aux_vrr(const I2CIntegral& integral) const;
};

#endif /* v3i_r2_driver_hpp */
