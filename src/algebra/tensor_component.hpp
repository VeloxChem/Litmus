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

#ifndef tensor_component_hpp
#define tensor_component_hpp

#include <string>
#include <optional>
#include <vector>

/// Tensor component class.
class TensorComponent
{
    /// Axial values of tensor component along X axis.
    int _ax;
    
    /// Axial values of tensor component along X axis.
    int _ay;
    
    /// Axial values of tensor component along X axis.
    int _az;
    
public:
    
    /// Creates a tensor component of zero order.
    TensorComponent();
    
    /// Creates a tensor component from the given axial values.
    /// @param ax The axial value along axis X to create tensor component.
    /// @param ay The axial value along axis Y to create tensor component.
    /// @param az The axial value along axis Z to create tensor component.
    TensorComponent(const int ax,
                    const int ay,
                    const int az);
    
    /// Retrieves axial value along requested axis.
    /// @param axis The axis to retrieve axial value.
    /// @return The axial value of tensor component along axis.
    int operator[](const char axis) const;
    
    /// Compares this tensor component with other tensor component.
    /// @param other The other tensor component to compare.
    /// @return true if tensor components are equal, false otherwise.
    bool operator==(const TensorComponent& other) const;
    
    /// Compares this tensor component with other tensor component.
    /// @param other The other tensor component to compare.
    /// @return true if tensor components are not equal, false otherwise.
    bool operator!=(const TensorComponent& other) const;
    
    /// Compares this tensor component with other tensor component.
    /// @param other The other tensor component to compare.
    /// @return true if this tensor component is less than other tensor component, false otherwise.
    bool operator<(const TensorComponent& other) const;
    
    /// Creates primitive textual representation of this tensor component.
    /// @return The string with primitive textual representation of tensor component.
    std::string to_string() const;
    
    /// Creates primitive textual label of this tensor component.
    /// @return The string with primitive textual label of tensor component.
    std::string label() const;
    
    /// Computes order of this tensor component.
    /// @return The order of tensor component.
    int order() const;
    
    /// Determines maximum axial value among all axes of this tensor component.
    /// @return The maximum axial value of tensor component.
    int maximum() const;
    
    /// Determines primary axis of this tensor component.
    /// @return: A primary axis of tensor component.
    char primary() const;
    
    /// Creates an optional tensor component from this tensor component by shifting axial value
    /// along the selected axis.
    /// @param axis The axis to shift axial value.
    /// @param value The value to shift axial value.
    /// @return The optional tensor component.
    std::optional<TensorComponent> shift(const char axis,
                                         const int  value) const;
};

using VTensorComponents = std::vector<TensorComponent>; 

#endif /* tensor_component_hpp */
