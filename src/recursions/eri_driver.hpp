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

#include "graph.hpp"
#include "recursion_group.hpp"
#include "recursion_expansion.hpp"
#include "recursion_term.hpp"
#include "integral_component.hpp"
#include "two_center_pair_component.hpp"

using T2CPair = TwoCenterPairComponent;

using T4CIntegral = IntegralComponent<T2CPair, T2CPair>;

using ST4CIntegrals = SIntegralComponents<T2CPair, T2CPair>;

using R4CTerm = RecursionTerm<T4CIntegral>;

using V4CTerms = VRecursionTerms<T4CIntegral>;

using R4CDist = RecursionExpansion<T4CIntegral>;

using R4Group = RecursionGroup<T4CIntegral>;

using V4Groups = VRecursionGroups<T4CIntegral>;

using R4Graph = Graph<R4Group>;

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
    
    /// Applies horizontal recursion to ket side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of horizontal recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R4CDist> ket_hrr(const R4CTerm& rterm,
                                   const char     axis) const;
    
    /// Applies vertical recursion to bra side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R4CDist> bra_vrr(const R4CTerm& rterm,
                                   const char     axis) const;
    
    /// Applies vertical recursion to ket side of given recursion term.
    /// @param rterm The recursion term.
    /// @param axis The axis of vertical recursion.
    /// @return The recursion expansion of given recursion term.
    std::optional<R4CDist> ket_vrr(const R4CTerm& rterm,
                                   const char     axis) const;
    
    /// Applies horizontal recursion to bra side recursion term containing electron repulsion
    /// integral.
    /// @param rterm The recursion term with electron repulsion integral.
    /// @param sints The set of electron repulsion integrals (updated by selecting recursion).
    /// @return The recursion expansion of given recursion term.
    R4CDist apply_bra_hrr(const R4CTerm&       rterm,
                                ST4CIntegrals& sints) const;
    
    /// Applies horizontal recursion to ket side recursion term containing electron repulsion
    /// integral.
    /// @param rterm The recursion term with electron repulsion integral.
    /// @param sints The set of electron repulsion integrals (updated by selecting recursion).
    /// @return The recursion expansion of given recursion term.
    R4CDist apply_ket_hrr(const R4CTerm&       rterm,
                                ST4CIntegrals& sints) const;
    
    /// Applies vertical recursion to bra side recursion term containing electron repulsion
    /// integral.
    /// @param rterm The recursion term with electron repulsion integral.
    /// @param sints The set of electron repulsion integrals (updated by selecting recursion).
    /// @return The recursion expansion of given recursion term.
    R4CDist apply_bra_vrr(const R4CTerm&       rterm,
                                ST4CIntegrals& sints) const;
    
    /// Applies vertical recursion to ket side recursion term containing electron repulsion
    /// integral.
    /// @param rterm The recursion term with electron repulsion integral.
    /// @param sints The set of electron repulsion integrals (updated by selecting recursion).
    /// @return The recursion expansion of given recursion term.
    R4CDist apply_ket_vrr(const R4CTerm&       rterm,
                                ST4CIntegrals& sints) const;
    
    /// Applies horizontal recursion to bra side of given vector of recursion terms containing
    /// electron repulsion integral.
    /// @param rterms The vector of recursion terms with electron repulsion integral.
    /// @param sints The set of electron repulsion integrals (updated by selecting recursion).
    /// @return The recursion group.
    R4Group apply_bra_hrr(const V4CTerms&      rterms,
                                ST4CIntegrals& sints) const;
    
    /// Applies horizontal recursion to ket side of given vector of recursion terms containing
    /// electron repulsion integral.
    /// @param rterms The vector of recursion terms with electron repulsion integral.
    /// @param sints The set of electron repulsion integrals (updated by selecting recursion).
    /// @return The recursion group.
    R4Group apply_ket_hrr(const V4CTerms&      rterms,
                                ST4CIntegrals& sints) const;
    
    /// Applies vertical recursion to bra side of given vector of recursion terms containing
    /// electron repulsion integral.
    /// @param rterms The vector of recursion terms with electron repulsion integral.
    /// @param sints The set of electron repulsion integrals (updated by selecting recursion).
    /// @return The recursion group.
    R4Group apply_bra_vrr(const V4CTerms&      rterms,
                                ST4CIntegrals& sints) const;
    
    /// Applies vertical recursion to ket side of given vector of recursion terms containing
    /// electron repulsion integral.
    /// @param rterms The vector of recursion terms with electron repulsion integral.
    /// @param sints The set of electron repulsion integrals (updated by selecting recursion).
    /// @return The recursion group.
    R4Group apply_ket_vrr(const V4CTerms&      rterms,
                                ST4CIntegrals& sints) const;
    
    /// Recursively applies horizontal recursion to bra side of given graph.
    /// @param rgraph The recursion graph.
    /// @param sints The set of electron repulsion integrals (updated by selecting recursion).
    void apply_bra_hrr(R4Graph&       rgraph,
                       ST4CIntegrals& sints) const;
    
    /// Recursively applies horizontal recursion to ket side of given graph.
    /// @param rgraph The recursion graph.
    /// @param sints The set of electron repulsion integrals (updated by selecting recursion).
    void apply_ket_hrr(R4Graph&       rgraph,
                       ST4CIntegrals& sints) const;
    
    /// Recursively applies vertical recursion to bra side of given graph.
    /// @param rgraph The recursion graph.
    /// @param sints The set of electron repulsion integrals (updated by selecting recursion).
    void apply_bra_vrr(R4Graph&       rgraph,
                       ST4CIntegrals& sints) const;
    
    /// Recursively applies vertical recursion to ket side of given graph.
    /// @param rgraph The recursion graph.
    /// @param sints The set of electron repulsion integrals (updated by selecting recursion).
    void apply_ket_vrr(R4Graph&       rgraph,
                       ST4CIntegrals& sints) const;
};

#endif /* eri_driver_hpp */