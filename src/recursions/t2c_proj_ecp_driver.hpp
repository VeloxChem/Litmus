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

#ifndef t2c_proj_ecp_driver_hpp
#define t2c_proj_ecp_driver_hpp

#include <optional>
#include <array>

#include "tensor_component.hpp"
#include "t2c_defs.hpp"

/// Two center projected ECP integrals driver class.
class T2CProjectedECPDriver
{
    /// Cartesian coordinate tensor components.
    std::array<TensorComponent, 3> _rxyz;
    
public:
    /// Creates a two center local ECP integrals driver.
    T2CProjectedECPDriver();
    
    /// Check if recursion term is for two-center projected ECP integral.
    /// @param rterm The recursion term.
    /// @return True if reccursion expansion belongs to overlap recursion, False otherwise.
    bool is_projected_ecp(const R2CTerm& rterm) const;
    
    /// Applies vertical recursion to bra side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R2CDist> bra_vrr(const R2CTerm& rterm,
                                   const char     axis) const;
    
    /// Applies vertical recursion to ket side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R2CDist> ket_vrr(const R2CTerm& rterm,
                                   const char     axis) const;
    
    /// Applies vertical recursion to bra side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R2CDist> common_bra_vrr(const R2CTerm& rterm,
                                          const char     axis) const;
    
    /// Applies vertical recursion to ket side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R2CDist> common_ket_vrr(const R2CTerm& rterm,
                                          const char     axis) const;
    
    /// Applies vertical recursion to bra side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R2CDist> red_bra_vrr(const R2CTerm& rterm,
                                       const char     axis) const;
    
    /// Applies vertical recursion to ket side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R2CDist> red_ket_vrr(const R2CTerm& rterm,
                                       const char     axis) const;
        
    /// Applies vertical recursion to bra side recursion term containing overlap
    /// integral.
    /// @param rterm The recursion term with overlap integral.
    /// @return The recursion expansion of given recursion term.
    R2CDist apply_bra_vrr(const R2CTerm& rterm) const;
    
    /// Applies vertical recursion to ket side recursion term containing overlap
    /// integral.
    /// @param rterm The recursion term with overlap integral.
    /// @return The recursion expansion of given recursion term.
    R2CDist apply_ket_vrr(const R2CTerm& rterm) const;
    
    /// Applies vertical recursion to bra side recursion term containing overlap
    /// integral.
    /// @param rterm The recursion term with overlap integral.
    /// @return The recursion expansion of given recursion term.
    R2CDist apply_red_bra_vrr(const R2CTerm& rterm) const;
    
    /// Applies vertical recursion to ket side recursion term containing overlap
    /// integral.
    /// @param rterm The recursion term with overlap integral.
    /// @return The recursion expansion of given recursion term.
    R2CDist apply_red_ket_vrr(const R2CTerm& rterm) const;
};

#endif /* t2c_proj_ecp_driver_hpp */
