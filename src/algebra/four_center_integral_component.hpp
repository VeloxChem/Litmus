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

#ifndef four_center_integral_component_hpp
#define four_center_integral_component_hpp

#include "operator_component.hpp"
#include "two_center_pair_component.hpp"

/// Four center integral component.
class FourCenterIntegralComponent
{
    /// Two center pair on bra side of integral component.
    TwoCenterPairComponent _bra_pair;
    
    /// Two center pair on ket side of integral component.
    TwoCenterPairComponent _ket_pair;
    
    /// Integrand operator of integral component.
    OperatorComponent _integrand;
    
    /// Order of integral component.
    int _order;
    
    ///  The vector of prefix operators acting on integral component.
    VOperatorComponents _prefixes;
    
public:
    /// Creates an empty four center integral component.
    FourCenterIntegralComponent();
    
    /// Creates a integral component from the given operator component and two center pair components.
    /// @param bra_pair The two center pair component on bra side of integral component.
    /// @param ket_pair The two center pair component on ket side of integral component.
    /// @param integrand The integrand operator component of integral component.
    /// @param order The order of integral component.
    /// @param prefixes The vector of prefix operator components acting on integral component.
    FourCenterIntegralComponent(const TwoCenterPairComponent& bra_pair,
                                const TwoCenterPairComponent& ket_pair,
                                const OperatorComponent&      integrand,
                                const int                     order = 0,
                                const VOperatorComponents&    prefixes = VOperatorComponents({}));
    
    /// Compares this integral component with other integral component.
    /// @param other The other integral component to compare.
    /// @return true if integral components, false otherwise.
    bool operator==(const FourCenterIntegralComponent& other) const;
    
    /// Compares this integral component with other integral component.
    /// @param other The other integral  component to compare.
    /// @return true if integral components are not equal, false otherwise.
    bool operator!=(const FourCenterIntegralComponent& other) const;
    
    /// Compares this integral component with other integral component.
    /// @param other The other integral component to compare.
    /// @return true if this integral  component is less than other integral component, false otherwise.
    bool operator<(const FourCenterIntegralComponent& other) const;
    
    /// Gets bra pair component of integral component.
    /// @return The bra pair component of integral component.
    TwoCenterPairComponent bra_pair() const;
    
    /// Gets ket pair component of integral component.
    /// @return The ket pair component of integral component.
    TwoCenterPairComponent ket_pair() const;
    
    /// Gets integrand of integral component.
    /// @return The integrand of integral component.
    OperatorComponent integrand() const;

    /// Gets order of integral component.
    /// @return The order of integral component.
    int order() const;
    
    /// Gets vector of prefix operator components  of integral component.
    /// @return The vector of prefix operator components of integral component.
    VOperatorComponents prefixes() const;
    
    /// Creates primitive textual label of this integral component.
    /// @param use_order The flag to include order of integral into its label.
    /// @return The string with primitive textual label of integral component.
    std::string label(const bool use_order = false) const;
    
    /// Creates an optional integral component from this integral component by shifting axial value
    /// along the selected axis on targeted center.
    /// @param axis The axis to shift axial value.
    /// @param value The value to shift axial value.
    /// @param center The targeted shift center.
    /// @return The optional integral component.
    std::optional<FourCenterIntegralComponent> shift(const char axis,
                                                     const int  value,
                                                     const int  center) const;
};

#endif /* four_center_integral_component_hpp */
