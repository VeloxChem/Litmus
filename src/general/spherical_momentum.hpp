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

#ifndef spherical_momentum_hpp
#define spherical_momentum_hpp

#include <string>
#include <vector>
#include <utility>

/// SphericalMomentum class.
class SphericalMomentum
{
    /// Vector  of scaling factors.
    std::vector<std::string> _factors;
    
    /// Vector of Cartesian indexes.
    std::vector<int> _indexes;
    
    /// Vector of dimensions.
    std::vector<int> _dimensions;
        
public:
    /// Creates an angular momentum.
    /// @param angmom The angular momentum.
    SphericalMomentum(const int angmom);
    
    /// Gets pairs of factors and spherical angular component  involving requested Cartesian angular component.
    /// @param cartcomp The Cartesian component to select.
    /// @return The vector of pairs (factor, spherical angular component)
    std::vector<std::pair<int, std::string>> select_pairs(const int cartcomp) const;
    
    /// Gets  vector of factors.
    /// @param angmom The angular momentum.
    /// @return The vector of factors.
    std::vector<std::string> get_factors(const int angmom) const;
};

#endif /* spherical_momentum_hpp */
