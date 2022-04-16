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

#ifndef tensor_hpp
#define tensor_hpp

#include <string>

#include "tensor_component.hpp"

/// Tensor class.
class Tensor
{
    /// Order of tensor.
    int _order;
    
public:
    /// Creates a tensor of zero order.
    Tensor();
    
    /// Creates a tensor of the given order.
    /// @param order The order to create tensor.
    Tensor(const int order);
    
    /// Creates a tensor from the given tensor component.
    /// @param tcomp The tensor component to create tensor.
    Tensor(const TensorComponent& tcomp);
    
    /// Compares this tensor with other tensor.
    /// @param other The other tensor to compare.
    /// @return true if tensors are equal, false otherwise.
    bool operator==(const Tensor& other) const;
    
    /// Compares this tensor with other tensor.
    /// @param other The other tensor to compare.
    /// @return true if tensors are not equal, false otherwise.
    bool operator!=(const Tensor& other) const;
    
    /// Compares this tensor with other tensor.
    /// @param other The other tensor to compare.
    /// @return true if this tensor is less than other tensor, false otherwise.
    bool operator<(const Tensor& other) const;
    
    /// Gets order of this tensor.
    /// @return The order of this tensor.
    int order() const;
        
    /// Creates primitive textual representation of this tensor.
    /// @return The string with primitive textual representation of tensor.
    std::string to_string() const;
    
    /// Creates label of this tensor.
    /// @return The string with label of tensor.
    std::string label() const;
    
    /// Creates a vector with tensor components of this tensor.
    /// @return The vector of tensor components.
    VTensorComponents components() const;
};

using VTensors = std::vector<Tensor>;

#endif /* tensor_hpp */
