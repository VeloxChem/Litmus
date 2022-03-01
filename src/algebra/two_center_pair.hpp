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

#ifndef two_center_pair_hpp
#define two_center_pair_hpp

#include <string>
#include <array>
#include <vector>

#include "tensor.hpp"

/// Two centers  pair.
class TwoCenterPair
{
    /// Names of expansion centers.
    std::array<std::string, 2> _names;
    
    /// Tensorial shapes of expansion centers.
    std::array<Tensor, 2> _shapes;
    
public:
    /// Creates an empty two centers pair.
    TwoCenterPair();
    
    /// Creates a two center pair from the given names and tensorial shapes.
    /// @param names The name of expansion centers.
    /// @param shapes The tensorial shapes of expansion centers.
    TwoCenterPair(const std::array<std::string, 2>& names,
                  const std::array<Tensor, 2>&      shapes);
    
    /// Creates a two center pair from the given two center pair component.
    /// @param t2pcomp The operator component to create operator.
    //TwoCenterPair(const TwoCenterPairComponent& t2pcomp);
    
    /// Compares this two center pair  with other two center pair.
    /// @param other The other two center pair to compare.
    /// @return true if two center pairs are equal, false otherwise.
    bool operator==(const TwoCenterPair& other) const;
    
    /// Compares this two center pair with other two center pair.
    /// @param other The other two center pair to compare.
    /// @return true if two center pairs are not equal, false otherwise.
    bool operator!=(const TwoCenterPair& other) const;
    
    /// Compares this two center pair with other two center pair.
    /// @param other The other two center pair to compare.
    /// @return true if this two center pair is less than other two center pair, false otherwise.
    bool operator<(const TwoCenterPair& other) const;
    
    /// Creates primitive textual representation of this two center pair.
    /// @return The string with primitive textual representation of two center pair.
    std::string to_string() const;
    
    /// Creates primitive textual label of this two center pair.
    /// @return The string with primitive textual label of two center pair.
    std::string label() const;
    
    /// Creates a vector with two center pair components of this two center pair.
    /// @return The vector of two center pair components.
    //VTwoCenterPairComponent components() const;
};

#endif /* two_center_pair_hpp */
