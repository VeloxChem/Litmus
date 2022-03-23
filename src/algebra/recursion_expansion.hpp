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

#include <vector>
#include <set>

#include <recursion_term.hpp>

/// Recursion expansion class.
template <class T>
class RecursionExpansion
{
    /// Root recursion term.
    RecursionTerm<T> _root;
    
    /// Vector of recursion terms in expansion of root recursion term..
    VRecursionTerms<T> _expansion;
  
public:
    /// Creates an empty recursion expansion.
    RecursionExpansion();
    
    /// Creates a  recursion expansion from the given recursion term and vector of recursion terms.
    /// @param root The root recursion term of recursion expansion.
    /// @param expansion The vector of recursion terms to expand root recursion term.
    RecursionExpansion(const RecursionTerm<T>&   root,
                       const VRecursionTerms<T>& expansion = VRecursionTerms<T>({}));
    
    /// Retrieves requested recursion term from recursion terms expansion of root recursion term.
    /// @param index The index of recursion term.
    /// @return The  requested recursion term.
    const RecursionTerm<T>& operator[](const int index) const;
    
    /// Compares this recursion expansion with other recursion expansion.
    /// @param other The other recursion expansion to compare.
    /// @return true if recursion expansions are equal, false otherwise.
    bool operator==(const RecursionExpansion<T>& other) const;
    
    /// Compares this recursion expansion with other recursion expansion.
    /// @param other The other recursion expansion to compare.
    /// @return true if recursion expansions are not equal, false otherwise.
    bool operator!=(const RecursionExpansion<T>& other) const;
    
    /// Compares this recursion expansion with other recursion expansion.
    /// @param other The other recursion expansion to compare.
    /// @return true if this recursion expansion is less than other recursion expansion, false otherwise.
    bool operator<(const RecursionExpansion<T>& other) const;
    
    /// Adds recursion term to expansion of root recursion expansion.
    /// @param rterm The recursion term to add.
    void add(const RecursionTerm<T>& rterm);
    
    /// Gets root recursion term of this recursion expansion.
    /// @return The root recursion term.
    RecursionTerm<T> root() const;
    
    /// Gets number of recursion terms in expansion of root recursion term.
    /// @return The number of recursion terms in expansion of root recursion term.
    int terms() const;
    
    /// Gets unique integrals in recursion expansion.
    /// @return The set of unique integrals in recursion expansion.
    std::set<T> unique_integrals() const;
    
    /// Counts number of new integrals in recursion expansion with respect to given integral set.
    /// @param integrals The set of integrals used to count new integrals.
    /// @return The number of new integrals in recursion expansion.
    int count_new_integrals(const std::set<T>& integrals) const;
    
    /// Checks if recursion expansion contains auxilary root  at specific
    /// angular center.
    /// @param center The angular center to check.
    /// @return True if recursion expansion contains auxilary root, false otherwise.
    bool auxilary(const int center) const;
};

template <class T>
RecursionExpansion<T>::RecursionExpansion()
    
    : _root()

    , _expansion(VRecursionTerms<T>({}))
{
    
}

template <class T>
RecursionExpansion<T>::RecursionExpansion(const RecursionTerm<T>&   root,
                                          const VRecursionTerms<T>& expansion)

    : _root(root)

    , _expansion(expansion)
{
    
}

template <class T>
const RecursionTerm<T>&
RecursionExpansion<T>::operator[](const int index) const
{
    return _expansion[index]; 
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
        return _expansion < other._expansion;
    }
}

template <class T>
void
RecursionExpansion<T>::add(const RecursionTerm<T>& rterm)
{
    _expansion.push_back(rterm); 
}

template <class T>
RecursionTerm<T>
RecursionExpansion<T>::root() const
{
    return _root; 
}

template <class T>
int
RecursionExpansion<T>::terms() const
{
    return static_cast<int>(_expansion.size());
}

template <class T>
std::set<T>
RecursionExpansion<T>::unique_integrals() const
{
    std::set<T> sints;
    
    for (const auto& rterm : _expansion)
    {
        sints.insert(rterm.integral()); 
    }
    
    return sints;
}

template <class T>
int
RecursionExpansion<T>::count_new_integrals(const std::set<T>& integrals) const
{
    int nints = 0;
    
    for (const auto& tint : unique_integrals())
    {
        if (integrals.find(tint) == integrals.cend()) nints++;
    }
    
    return nints;
}

template <class T>
bool
RecursionExpansion<T>::auxilary(const int center) const
{
    return _root.auxilary(center);
}

template <class T>
using VRecursionExpansions = std::vector<RecursionExpansion<T>>; 

#endif /* recursion_expansion_hpp */
