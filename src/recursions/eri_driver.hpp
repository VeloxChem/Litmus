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

#ifndef eri_driver_hpp
#define eri_driver_hpp

#include <optional>
#include <set>

#include "recursion_expansion.hpp"
#include "recursion_term.hpp"
#include "integral_component.hpp"
#include "two_center_pair_component.hpp"

using T2CPair = TwoCenterPairComponent;

using T4CIntegral = IntegralComponent<T2CPair, T2CPair>;

using R4CTerm = RecursionTerm<T4CIntegral>;

using R4CDist = RecursionExpansion<T4CIntegral>;

/// Electron repulsion integrals driver class.
class EriDriver
{
    /// Cartesian coordinate tensor components.
    std::array<TensorComponent, 3> _rxyz;
    
public:
    /// Creates an electron repulsion integrals driver.
    EriDriver();
    
    /// Applies horizontal recursion to bra side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of horizontal recursion. 
    /// @return The recursion expansion of given recursion term.
    std::optional<R4CDist> bra_hrr(const R4CTerm& rterm,
                                   const char     axis) const;
    
    /// Applies horizontal recursion to bra side recursion term containing electron repulsion
    /// integral.
    /// @param rterm The recursion term with electron repulsion integral.
    /// @param vints The set of electron repulsion integrals (updated by selecting recursion).
    /// @return The recursion expansion of given recursion term.
    R4CDist apply_bra_hrr(const R4CTerm&         rterm,
                          std::set<T4CIntegral>& vints) const;
};

#endif /* eri_driver_hpp */
