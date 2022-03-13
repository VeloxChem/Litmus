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

#ifndef four_center_integral_hpp
#define four_center_integral_hpp

#include <vector>

#include "operator.hpp"
#include "two_center_pair.hpp"
#include "four_center_integral_component.hpp"

/// Four center integral.
class FourCenterIntegral
{
    /// Two center pair on bra side of integral.
    TwoCenterPair _bra_pair;
    
    /// Two center pair on ket side of integral.
    TwoCenterPair _ket_pair;
    
    /// Integrand operator of integral.
    Operator _integrand;
    
    /// Order of integral.
    int _order;
    
    ///  The vector of prefix operators acting on integral.
    VOperators _prefixes;
    
public:
    /// Creates an empty four center integral.
    FourCenterIntegral();
    
    /// Creates a integral from the given operator and two center pairs.
    /// @param bra_pair The two center pair on bra side of integral.
    /// @param ket_pair The two center pair on ket side of integral.
    /// @param integrand The integrand operator of integral.
    /// @param order The order of integral.
    /// @param prefixes The vector of prefix operators acting on integral.
    FourCenterIntegral(const TwoCenterPair& bra_pair,
                       const TwoCenterPair& ket_pair,
                       const Operator&      integrand,
                       const int            order,
                       const VOperators&    prefixes);
    
    /// Creates a four center Gaussian integral with given operator and angular momentum values.
    /// @param a_angmom The angular momentum of center GA in integral.
    /// @param b_angmom The angular momentum of center GB in integral.
    /// @param c_angmom The angular momentum of center GC in integral.
    /// @param d_angmom The angular momentum of center GD in integral.
    /// @param integrand The integrand operator of integral.
    /// @param order The order of integral.
    /// @param prefixes The vector of prefix operators acting on integral.
    FourCenterIntegral(const int          a_angmom,
                       const int          b_angmom,
                       const int          c_angmom,
                       const int          d_angmom,
                       const Operator&    integrand,
                       const int          order = 0,
                       const VOperators&  prefixes = VOperators({}));
    
    /// Creates a integral from the given integral component.
    /// @param t4ccomp The integral component to create integral.
    FourCenterIntegral(const FourCenterIntegralComponent& t4ccomp);
    
    /// Compares this integral with other integral.
    /// @param other The other integral to compare.
    /// @return true if integrals, false otherwise.
    bool operator==(const FourCenterIntegral& other) const;
    
    /// Compares this integral with other integral.
    /// @param other The other integral to compare.
    /// @return true if integrals are not equal, false otherwise.
    bool operator!=(const FourCenterIntegral& other) const;
    
    /// Compares this integral with other integral.
    /// @param other The other integral to compare.
    /// @return true if this integral is less than other integral, false otherwise.
    bool operator<(const FourCenterIntegral& other) const;
    
    /// Creates primitive textual label of this integral.
    /// @return The string with primitive textual label of integral.
    std::string label(const bool use_order = false) const;
    
    /// Creates a vector with integral components of this integral.
    /// @param only_diag The flag to indicate generation of only diagonal integral components.
    /// @return The vector of integral components.
    VFourCenterIntegralComponents components(const bool only_diag = false) const;
};

#endif /* four_center_integral_hpp */
