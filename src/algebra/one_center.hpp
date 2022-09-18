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

#ifndef one_center_hpp
#define one_center_hpp

#include "tensor.hpp"
#include "two_center_pair_component.hpp"

#include <string>

#include "tensor.hpp"

/// One center expansion.
class OneCenter
{
    /// Names of expansion center.
    std::string _name;
    
    /// Tensorial shapes of expansion center.
    Tensor _shape;
    
public:
    /// Creates an empty one centers expansion.
    OneCenter();
    
    /// Creates a one center expansion from the given name and tensorial shape.
    /// @param name The name of expansion center.
    /// @param shape The tensorial shape of expansion center.
    OneCenter(const std::string& name,
              const Tensor&      shape);
    
    /// Creates a one center expansion from the given name and angular momentum.
    /// @param name The name of expansion center.
    /// @param angmom The angular momentum of expansion center.
    OneCenter(const std::string& name,
              const int          angmom);
    
    /// Creates a two center pair from the given two center pair component.
    /// @param t2pcomp The two center pair component to create two center pair.
    //TwoCenterPair(const TwoCenterPairComponent& t2pcomp);
    
    /// Compares this one center expansion with other one center expansion.
    /// @param other The other one center expansion to compare.
    /// @return true if one center expansions are equal, false otherwise.
    bool operator==(const OneCenter& other) const;
    
    /// Compares this one center expansion with other one center expansion.
    /// @param other The other one center expansion to compare.
    /// @return true if one center expansions are not equal, false otherwise.
    bool operator!=(const OneCenter& other) const;
    
    /// Compares this one center expansion with other one center expansion.
    /// @param other The other one center expansion to compare.
    /// @return true if this one center expansion is less than other one center expansion, false otherwise.
    bool operator<(const OneCenter& other) const;
    
    /// Gets tensorial shape of this one center expansion.
    /// @return The tensorial shape of one center expansion.
    int shape() const {return _shape.order();};
    
    /// Creates primitive textual representation of this one center expansion.
    /// @return The string with primitive textual representation of one center expansion.
    std::string to_string() const;
    
    /// Creates primitive textual label of this one center expansion.
    /// @return The string with primitive textual label of one center expansion.
    std::string label() const;
    
    /// Creates a vector with two center pair components of this two center pair.
    /// @return The vector of two center pair components.
    //VTwoCenterPairComponents components() const;
};

#endif /* one_center_hpp */
