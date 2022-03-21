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

#ifndef recursion_group_hpp
#define recursion_group_hpp

#include <recursion_expansion.hpp>

/// Recursion group class.
template <class T>
class RecursionGroup
{
    /// Group of recursion expansions.
    VRecursionExpansions<T> _expansions;
  
public:
    /// Creates an empty recursion group.
    RecursionGroup();
    
    /// Creates a  recursion group from the vector of recursion expansions.
    /// @param expansions The vector of recursion expansions.
    RecursionGroup(const VRecursionExpansions<T>& expansions);
    
    /// Retrieves requested recursion expansion from recursion group.
    /// @param index The index of recursion expansion.
    /// @return The  requested recursion expansion.
    const RecursionExpansion<T>& operator[](const int index) const;
    
    /// Compares this recursion group with other recursion group.
    /// @param other The other recursion group to compare.
    /// @return true if recursion groups are equal, false otherwise.
    bool operator==(const RecursionGroup<T>& other) const;
    
    /// Compares this recursion group with other recursion group.
    /// @param other The other recursion group to compare.
    /// @return true if recursion groups are not equal, false otherwise.
    bool operator!=(const RecursionGroup<T>& other) const;
    
    /// Compares this recursion group with other recursion group.
    /// @param other The other recursion group to compare.
    /// @return true if this recursion group is less than other recursion group, false otherwise.
    bool operator<(const RecursionGroup<T>& other) const;
    
    /// Adds recursion expansion to recursion expansion group.
    /// @param expansion The recursion expansion to add.
    void add(const RecursionExpansion<T>& expansion);
    
    /// Gets number of recursion expansions in recursion group.
    /// @return The number of recursion expansions in recursion group.
    int expansions() const;
    
    /// Splits expansions in recursion group into 2D vector of unique recursion
    /// terms according to integral type.
    /// @return The vector of integral groups.
    template <class U>
    MRecursionTerms<T> split_terms() const;
};

template <class T>
RecursionGroup<T>::RecursionGroup()
    
    : _expansions(VRecursionExpansions<T>({}))
{
    
}

template <class T>
RecursionGroup<T>::RecursionGroup(const VRecursionExpansions<T>& expansions)

    : _expansions(expansions)
{
    
}

template <class T>
const RecursionExpansion<T>&
RecursionGroup<T>::operator[](const int index) const
{
    return _expansions[index];
}

template <class T>
bool
RecursionGroup<T>::operator==(const RecursionGroup<T>& other) const
{
    if (this == &other) return true;

    return _expansions == other._expansions;
}

template <class T>
bool
RecursionGroup<T>::operator!=(const RecursionGroup<T>& other) const
{
    return !((*this) == other);
}

template <class T>
bool
RecursionGroup<T>::operator<(const RecursionGroup<T>& other) const
{
    return _expansions < other._expansions;
}

template <class T>
void
RecursionGroup<T>::add(const RecursionExpansion<T>& expansion)
{
    _expansions.push_back(expansion);
}

template <class T>
int
RecursionGroup<T>::expansions() const
{
    return static_cast<int>(_expansions.size());
}

template <class T>
template <class U>
MRecursionTerms<T>
RecursionGroup<T>::split_terms() const
{
    // collect unique integral components
    
    std::set<T> sints;
    
    for (const auto& tval : _expansions)
    {
        const auto rints = tval.unique_integrals();
        
        sints.insert(rints.cbegin(), rints.cend());
    }
    
    // collect unique integrals
    
    std::set<U> tints;
   
    for(const auto& tval: sints)
    {
        tints.insert(U(tval));
    }
   
    // create matrix of recursion terms
    
    MRecursionTerms<T> mrterms(tints.size(), VRecursionTerms<T>());
    
    int idx = 0;
    
    for (const auto& tint : tints)
    {
        for (const auto tval : sints)
        {
            if (tint == U(tval))
            {
                mrterms[idx].push_back(RecursionTerm<T>(tval));
            }
        }
        
        idx++;
    }
    
    return mrterms;
}

#endif /* recursion_group_hpp */
