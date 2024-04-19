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

#ifndef factor_hpp
#define factor_hpp

#include <string>

#include "tensor_component.hpp"

/// Factor class.
class Factor
{
    /// Name of factor.
    std::string _name;
    
    /// Label of factor.
    std::string _label;
    
    /// Tensorial shape of factor.
    TensorComponent _shape;
        
public:
    /// Creates an empy factor.
    Factor();
    
    /// Creates a factror from given name, label and tensorial shape.
    /// @param name The name of factor.
    /// @param label The label of factor.
    /// @param shape The tensorial shape of factor.
    Factor(const std::string&     name,
           const std::string&     label,
           const TensorComponent& shape = TensorComponent(0, 0, 0));
    
    /// Compares this factor with other factor.
    /// @param other The other factor to compare.
    /// @return true if factors are equal, false otherwise.
    bool operator==(const Factor& other) const;
    
    /// Compares this factor with other factor.
    /// @param other The other factor to compare.
    /// @return true if factors are not equal, false otherwise.
    bool operator!=(const Factor& other) const;
    
    /// Compares this factor with other factor.
    /// @param other The other factor to compare.
    /// @return true if this factor is less than other factor, false otherwise.
    bool operator<(const Factor& other) const;
    
    /// Gets tensorial order of factor.
    /// @return The tensorial order of factor.
    int order() const;
    
    /// Creates primitive textual representation of this factor.
    /// @return The string with primitive textual representation of factor.
    std::string to_string() const;
    
    /// Creates label of this factor.
    /// @param nocomp The flag to include tensor component. 
    /// @return The string with label of factor.
    std::string label(const bool nocomp = false) const;
    
    /// Gets name of this factor.
    /// @return The string with name of factor.
    std::string name() const;
};

#endif /* factor_hpp */
