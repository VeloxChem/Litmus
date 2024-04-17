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

#ifndef one_center_component_hpp
#define one_center_component_hpp

#include <string>

#include "tensor_component.hpp"

/// One center  expansion component.
class OneCenterComponent
{
    /// Name of expansion center.
    std::string _name;
    
    /// Tensorial shape of expansion center.
    TensorComponent _shape;
    
public:
    /// Creates an empty one center expansion component.
    OneCenterComponent();
    
    /// Creates a one center expansion component from the given name and tensorial shape.
    /// @param name The name of expansion center.
    /// @param shape The tensorial shape of expansion center.
    OneCenterComponent(const std::string&     name,
                       const TensorComponent& shape);
    
    /// Retrieves axial value along requested center.
    /// @param index The index of center to retrieve axial value.
    /// @return The  tensorial shape of requested center in one center expansion component.
    const TensorComponent& operator[](const int index) const;
    
    /// Compares this one center expansion component with other one center expansion component.
    /// @param other The other one center expansion component to compare.
    /// @return true if one center expansion components are equal, false otherwise.
    bool operator==(const OneCenterComponent& other) const;
    
    /// Compares this one center expansion component  with other one center expansion component.
    /// @param other The other one center expansion component to compare.
    /// @return true if one center expansion components are not equal, false otherwise.
    bool operator!=(const OneCenterComponent& other) const;
    
    /// Compares this one center expansion component with other one center expansion component.
    /// @param other The other one center expansion component to compare.
    /// @return true if this one center expansion component  is less than other one center expansion component, false otherwise.
    bool operator<(const OneCenterComponent& other) const;
    
    /// Checks if this one center expansion component is similar to other one center expansion component.
    /// @param other The other one center expansion component to compare.
    /// @return True if one center expansion components  are similar, false otherwise.
    bool similar(const OneCenterComponent& other) const;
    
    /// Gets name of one center expansion component.
    /// @return The name of one center expansion component.
    std::string name() const {return _name;};
    
    /// Gets shape of one center expansion component.
    /// @return The shape of one center expansion component.
    TensorComponent shape() const {return _shape;};
    
    /// Gets number of centers  in one center expansion component.
    /// @return The number of centers in one center expansion component.
    int centers() const {return 1;};
    
    /// Creates primitive textual representation of this one center expansion component.
    /// @return The string with primitive textual representation of one center expansion component.
    std::string to_string() const;
    
    /// Creates primitive textual label of this one center expansion component.
    /// @return The string with primitive textual label of one center expansion component.
    std::string label() const;
    
    /// Creates an optional one center expansion component from this one center expansion component by shifting axial value
    /// along the selected axis on targeted center.
    /// @param axis The axis to shift axial value.
    /// @param value The value to shift axial value.
    /// @param center The targeted shift center.
    /// @return The optional one center expansion component.
    std::optional<OneCenterComponent> shift(const char axis,
                                            const int  value,
                                            const int  center) const;
};

using VOneCenterComponents = std::vector<OneCenterComponent>;

#endif /* one_center_component_hpp */
