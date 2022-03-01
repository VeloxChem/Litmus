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

#ifndef operator_component_hpp
#define operator_component_hpp

#include <string>
#include <optional>
#include <vector>

#include "tensor_component.hpp"

/// Operator component class.
class OperatorComponent
{
    /// Name of operator component.
    std::string _name;
    
    /// Tensorial shape of operator component.
    TensorComponent _shape;
    
    /// The target of operator component.
    std::string _target;
    
    /// The targeted center of operator component.
    int _center;
    
public:
    /// Creates an empty operator component.
    OperatorComponent();
    
    /// Creates an operator component  from the given name and tensorial shapes.
    /// @param name The name to create operator component.
    /// @param shape The tensorial shape to create operator component.
    /// @param target The target of operator component action.
    /// @param center The targeted center of operator component action.
    OperatorComponent(const std::string&     name,
                      const TensorComponent& shape  = TensorComponent(0, 0, 0),
                      const std::string&     target = "none",
                      const int              center = -1);
    
    /// Compares this operator component with other operator component.
    /// @param other The other operator component to compare.
    /// @return true if operator components are equal, false otherwise.
    bool operator==(const OperatorComponent& other) const;
    
    /// Compares this operator component with other operator component.
    /// @param other The other operator component to compare.
    /// @return true if operator components are not equal, false otherwise.
    bool operator!=(const OperatorComponent& other) const;
    
    /// Compares this operator component  with other operator component.
    /// @param other The other operator component to compare.
    /// @return true if this operator component is less than other operator component, false otherwise.
    bool operator<(const OperatorComponent& other) const;
    
    /// Gets name of operator component.
    /// @return The name of operator component.
    std::string name() const;
    
    /// Gets tensorial shape of operator component.
    /// @return The tensorial shape of operator component.
    TensorComponent shape() const;
    
    /// Gets target of operator component action.
    /// @return The target of operator component action.
    std::string target() const;
    
    /// Gets center of target of operator component action.
    /// @return The center of target of operator component action.
    int center() const;
    
    /// Creates primitive textual representation of this operator component.
    /// @return The string with primitive textual representation of operator component.
    std::string to_string() const;
    
    /// Creates primitive textual label of this operator component.
    /// @return The string with primitive textual label of operator component.
    std::string label() const;
    
    /// Creates an optional operator component from this operator component by shifting axial value
    /// along the selected axis.
    /// @param center The tensorial center to shift.
    /// @param axis The axis to shift axial value.
    /// @param value The value to shift axial value.
    /// @param noscalar The flag for scalar component: false  to keep,  true otherwise.
    /// @return The optional operator component.
    std::optional<OperatorComponent> shift(const int  center,
                                           const char axis,
                                           const int  value,
                                           const bool noscalar = false) const;
};

using VOperatorComponents = std::vector<OperatorComponent>;

#endif /* operator_component_hpp */
