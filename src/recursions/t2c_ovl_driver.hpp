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

#ifndef t2c_ovl_driver_hpp
#define t2c_ovl_driver_hpp

#include <optional>
#include <set>
#include <array>

#include "tensor_component.hpp"
#include "t2c_defs.hpp"

/// Two center overlap integrals driver class.
class T2COverlapDriver
{
    /// Cartesian coordinate tensor components.
    std::array<TensorComponent, 3> _rxyz;
    
public:
    /// Creates a two center overlap integrals driver.
    T2COverlapDriver();
    
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
        
    /// Applies vertical recursion to bra side recursion term containing overlap
    /// integral.
    /// @param rterm The recursion term with overlap integral.
    /// @param sints The set of overlap integrals (updated by selecting recursion).
    /// @return The recursion expansion of given recursion term.
    R2CDist apply_bra_vrr(const R2CTerm&       rterm,
                                ST2CIntegrals& sints) const;
    
    /// Applies vertical recursion to ket side recursion term containing overlap
    /// integral.
    /// @param rterm The recursion term with overlap integral.
    /// @param sints The set of overlap integrals (updated by selecting recursion).
    /// @return The recursion expansion of given recursion term.
    R2CDist apply_ket_vrr(const R2CTerm&       rterm,
                                ST2CIntegrals& sints) const;
        
    /// Applies vertical recursion to bra side of given vector of recursion terms containing
    /// overlap integral.
    /// @param rterms The vector of recursion terms with overlap integral.
    /// @param sints The set of overlap integrals (updated by selecting recursion).
    /// @return The recursion group.
    R2Group apply_bra_vrr(const V2CTerms&      rterms,
                                ST2CIntegrals& sints) const;
    
    /// Applies vertical recursion to ket side of given vector of recursion terms containing
    /// overlap integral.
    /// @param rterms The vector of recursion terms with overlap integral.
    /// @param sints The set of overlap integrals (updated by selecting recursion).
    /// @return The recursion group.
    R2Group apply_ket_vrr(const V2CTerms&      rterms,
                                ST2CIntegrals& sints) const;
    
    /// Recursively applies vertical recursion to bra side of given recursion groups container.
    /// @param rgroups The recursion groups container.
    /// @param sints The set of overlap integrals (updated by selecting recursion).
    void apply_bra_vrr(R2GroupContainer& rgroups,
                       ST2CIntegrals&    sints) const;
    
    /// Recursively applies vertical recursion to ket side of given recursion groups container.
    /// @param rgroups The recursion groups container.
    /// @param sints The set of overlap integrals (updated by selecting recursion).
    void apply_ket_vrr(R2GroupContainer& rgroups,
                       ST2CIntegrals&    sints) const;
    
    /// Recursively applies Obara-Saika recursion to recursion groups container.
    /// @param rgroups The recursion groups container.
    /// @param sints The set of overlap integrals (updated by selecting recursion).
    void apply_recursion(R2GroupContainer& rgroups,
                         ST2CIntegrals&    sints) const;
    
    /// Creates recursion graph for given angular momentum values.
    /// @param anga The angular momentum of center A.
    /// @param angb The angular momentum of center B.
    /// @param angc The angular momentum of center C.
    /// @param angd The angular momentum of center D.
    /// @param diag The flag to generate only diagonal terms.
    /// @return The recursion graph.
    //R4Graph create_graph(const int  anga,
    //                     const int  angb,
    //                     const int  angc,
    //                     const int  angd,
    //                     const bool diag) const;
    
    /// Creates vector of recursion graphs with upper limit of angular momentum values.
    /// @param mang The maximum angular momentum.
    /// @param diag The flag to generate only diagonal terms
    /// @return The vector of recursion graphs.
    //V4Graphs create_graphs(const int  mang,
    //                       const bool diag) const;
};

#endif /* t2c_ovl_driver_hpp */
