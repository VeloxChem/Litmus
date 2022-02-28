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

#ifndef operator_hpp
#define operator_hpp

#include <string>
#include <vector>

#include "tensor.hpp"
#include "operator_component.hpp"

/// Operator class.
class Operator
{
    /// Name of operator component.
    std::string _name;
    
    /// Tensorial shape of operator.
    Tensor _shape;
    
    /// The target of operator.
    std::string _target;
    
    /// The targeted center of operator.
    int _center;
    
public:
    /// Creates an empty operator.
    Operator();
    
    /// Creates an operator from the given name and tensorial shapes.
    /// @param name The name to create operator.
    /// @param shape The tensorial shape to create operator.
    /// @param target The target of operator action.
    /// @param center The targeted center of operator action.
    Operator(const std::string& name,
             const Tensor&      shape  = Tensor(0),
             const std::string& target = std::string("self"),
             const int          center = 0);
    
    /// Creates an operator from the given operator component.
    /// @param opcomp The operator component to create operator.
    Operator(const OperatorComponent& opcomp);
    
    /// Compares this operator  with other operator.
    /// @param other The other operator to compare.
    /// @return true if operators are equal, false otherwise.
    bool operator==(const Operator& other) const;
    
    /// Compares this operator with other operator.
    /// @param other The other operator to compare.
    /// @return true if operators are not equal, false otherwise.
    bool operator!=(const Operator& other) const;
    
    /// Compares this operator with other operator.
    /// @param other The other operator to compare.
    /// @return true if this operator is less than other operator, false otherwise.
    bool operator<(const Operator& other) const;
    
    /// Creates primitive textual representation of this operator.
    /// @return The string with primitive textual representation of operator.
    std::string to_string() const;
    
    /// Creates primitive textual label of this operator.
    /// @return The string with primitive textual label of operator.
    std::string label() const;
    
    /// Creates a vector with operator components of this operator.
    /// @return The vector of operator components.
    VOperatorComponents components() const;
};

using VOperators = std::vector<Operator>;

#endif /* operator_hpp */
