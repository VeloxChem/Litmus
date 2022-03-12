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

#ifndef two_center_pair_component_hpp
#define two_center_pair_component_hpp

#include <string>
#include <array>

#include "tensor_component.hpp"

/// Two centers  pair component.
class TwoCenterPairComponent
{
    /// Names of expansion centers.
    std::array<std::string, 2> _names;
    
    /// Tensorial shapes of expansion centers.
    std::array<TensorComponent, 2> _shapes;
    
public:
    /// Creates an empty two centers pair component.
    TwoCenterPairComponent();
    
    /// Creates a two center pair component from the given names and tensorial shapes.
    /// @param names The name of expansion centers.
    /// @param shapes The tensorial shapes of expansion centers.
    TwoCenterPairComponent(const std::array<std::string, 2>&     names,
                           const std::array<TensorComponent, 2>& shapes);
    
    /// Retrieves axial value along requested center.
    /// @param index The index of center to retrieve axial value.
    /// @return The  tensorial shape of requested center in two center pair component.
    const TensorComponent& operator[](const int index) const;
    
    /// Compares this two center pair component with other two center pair component.
    /// @param other The other two center pair component to compare.
    /// @return true if two center pair components are equal, false otherwise.
    bool operator==(const TwoCenterPairComponent& other) const;
    
    /// Compares this two center pair component  with other two center pair component.
    /// @param other The other two center pair component to compare.
    /// @return true if two center pair components are not equal, false otherwise.
    bool operator!=(const TwoCenterPairComponent& other) const;
    
    /// Compares this two center pair component with other two center pair component.
    /// @param other The other two center pair component to compare.
    /// @return true if this two center pair component  is less than other two center pair component, false otherwise.
    bool operator<(const TwoCenterPairComponent& other) const;
    
    /// Gets names of centers in two center pair component.
    /// @return The names of centers.
    std::array<std::string, 2> names() const;
    
    /// Gets shapes of centers in two center pair component.
    /// @return The shapes of centers.
    std::array<TensorComponent, 2> shapes() const;
    
    /// Creates primitive textual representation of this two center pair component.
    /// @return The string with primitive textual representation of two center pair component.
    std::string to_string() const;
    
    /// Creates primitive textual label of this two center pair component.
    /// @return The string with primitive textual label of two center pair component.
    std::string label() const;
    
    /// Creates an optional two center pair component from this two center pair component by shifting axial value
    /// along the selected axis on targeted center.
    /// @param axis The axis to shift axial value.
    /// @param value The value to shift axial value.
    /// @param center The targeted shift center.
    /// @return The optional two center pair component.
    std::optional<TwoCenterPairComponent> shift(const char axis,
                                                const int  value,
                                                const int  center) const;
};

using VTwoCenterPairComponents = std::vector<TwoCenterPairComponent>;

#endif /* two_center_pair_component_hpp */
