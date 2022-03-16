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

#ifndef recursion_expansion_hpp
#define recursion_expansion_hpp

#include <string>
#include <vector>

#include <recursion_term.hpp>

/// Recursion expansion class.
template <class T>
class RecursionExpansion
{
    /// Root recursion term.
    T _root;
    
    /// Vector of recursion terms in expansion of root recursion term..
    VRecursionTerms<T> _expansion;
  
public:
    /// Creates an empty recursion expansion.
    RecursionExpansion();
    
    /// Creates a  recursion expansion from the given recursion term and vector of recursion terms.
    /// @param root The root recursion term of recursion expansion.
    /// @param expansion The vector of recursion terms to expand root recursion term.
    RecursionExpansion(const T&                  root,
                       const VRecursionTerms<T>& expansion = VRecursionTerms<T>({}));
    
    /// Compares this recursion expansion with other recursion expansion.
    /// @param other The other recursion expansion to compare.
    /// @return true if recursion expansions are equal, false otherwise.
    bool operator==(const RecursionExpansion& other) const;
    
    /// Compares this recursion expansion with other recursion expansion.
    /// @param other The other expansion term to compare.
    /// @return true if recursion expansionss are not equal, false otherwise.
    bool operator!=(const RecursionExpansion& other) const;
    
    /// Compares this recursion expansion with other recursion expansion.
    /// @param other The other recursion expansion to compare.
    /// @return true if this recursion expansion is less than other recursion expansion, false otherwise.
    bool operator<(const RecursionTerm& other) const;
};

template <class T>
RecursionExpansion<T>::RecursionExpansion()
    
    : _root()

    , _expansion(std::vector<RecursionTerm<T>>({}))
{
    
}

template <class T>
RecursionExpansion<T>::RecursionExpansion(const T&                  root,
                                          const VRecursionTerms<T>& expansion)

    : _root()

    , _expansion(expansion)
{
    
}

template <class T>
bool
RecursionExpansion<T>::operator==(const RecursionExpansion<T>& other) const
{
    if (this == &other) return true;

    if (_root != other._root)
    {
        return false;
    }
    else
    {
        return _expansion == other._expansion;
    }
}

template <class T>
bool
RecursionExpansion<T>::operator!=(const RecursionExpansion<T>& other) const
{
    return !((*this) == other);
}

template <class T>
bool
RecursionExpansion<T>::operator<(const RecursionExpansion<T>& other) const
{
    if (_root != other._root)
    {
        return _root < other._root;
    }
    else
    {
        return _expasion < other._expansion;
    }
}

#endif /* recursion_expansion_hpp */
