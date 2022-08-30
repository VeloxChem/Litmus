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

#include <vector>
#include <algorithm>
#include <optional>
#include <climits>

#include "recursion_expansion.hpp"
#include "signature.hpp"

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
    const RecursionExpansion<T>& operator[](const size_t index) const;
    
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
    
    /// Checks if this recursion group is similar to other recursion group.
    /// @param other The other recursion group to compare.
    /// @return True if recursion groups  are similar, false otherwise.
    bool similar(const RecursionGroup<T>& other) const;
    
    /// Adds recursion expansion to recursion expansion group.
    /// @param expansion The recursion expansion to add.
    void add(const RecursionExpansion<T>& expansion);
    
    /// Merges given recursion expansion into this recursion  group.
    /// @param other The recursion group to merge.
    void merge(const RecursionGroup<T>& other);
    
    /// Reduces order of this recursion term to minimal possible.
    void reduce();
    
    /// Cheks if recursion group contains the given recursion expansion.
    /// @param rexp The recursion expansion to check.
    /// @return Ture if recursion group contains given recursion expansion, false otherwise.
    bool contains(const RecursionExpansion<T>& rexp) const;
    
    /// Gets number of recursion expansions in recursion group.
    /// @return The number of recursion expansions in recursion group.
    size_t expansions() const;
    
    /// Splits expansions in recursion group into 2D vector of unique recursion
    /// terms according to integral type.
    /// @return The 2D vector of unique recursion terms.
    template <class U>
    MRecursionTerms<T> split_terms() const;
    
    /// Generates vector of recursions terms from roots of expansions.
    /// @return The vector of recursion terms.
    VRecursionTerms<T> roots() const;
    
    /// Checks if recursion group contains only empty recursion expansions.
    /// @return True if recursion group contains only empty recursion expansions,
    /// false otherwise.
    bool empty() const;
    
    /// Checks if recursion group contains only auxilary recursion expansions on specific
    /// angular center.
    /// @param center The angular center to check.
    /// @return True if recursion group contains only auxilary recursion expansions,
    /// false otherwise.
    bool auxilary(const int center) const;
    
    /// Gets optional base integral of unniform recursion group.
    /// @return The base integral of uniform recursion group.
    template <class U>
    std::optional<U> base() const;
    
    /// Determines minimum order of integrals in this recursion group.
    /// @return The minimum order of integrlas.
    int min_order() const;
    
    /// Creates signature of recursion group.
    /// @return The signature of recursion group.
    Signature<T> signature() const;
    
    /// Gets unique factors in this recursion group.
    /// @return The unique recursion factors.
    std::set<Factor> factors() const;
    
    /// Gets unique prefactors in this recursion group.
    /// @return The unique prefactors.
    std::set<Fraction> prefactors() const;
    
    /// Gets map of factors in recursion group.
    /// @return The map of factors in recursion group.
    std::map<Factor, int> map_of_factors() const;
    
    /// Gets set of unique integral components in recursion group.
    /// @return The set of unique integral components in recursion group.
    std::set<T> components() const;
    
    /// Gets set of unique integral components of specific integral in recursion group.
    /// @param integral The integral to match. 
    /// @return The set of unique integral components of  specific integral  in recursion group.
    template <class U>
    std::set<T> components(const U& integral) const;
    
    /// Gets set of unique integrals  in recursion group.
    /// @return The set of unique integrals in recursion group.
    template <class U>
    std::set<U> integrals() const;
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
        std::sort(_expansions.begin(), _expansions.end());
}

template <class T>
const RecursionExpansion<T>&
RecursionGroup<T>::operator[](const size_t index) const
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
bool
RecursionGroup<T>::similar(const RecursionGroup<T>& other) const
{
    for (const auto& rhsrt : _expansions)
    {
        for (const auto lhsrt : other._expansions)
        {
            if (!lhsrt.similar(rhsrt))
            {
                return false;
            }
        }
    }
    
    return true; 
}

template <class T>
void
RecursionGroup<T>::add(const RecursionExpansion<T>& expansion)
{
    _expansions.push_back(expansion);
    
    std::sort(_expansions.begin(), _expansions.end());
}

template <class T>
void
RecursionGroup<T>::reduce()
{
    const auto morder = min_order();
    
    if (morder > 0)
    {
        for (size_t i = 0; i < _expansions.size(); i++)
        {
            _expansions[i].reduce(morder); 
        }
    }
}

template <class T>
bool
RecursionGroup<T>::contains(const RecursionExpansion<T>& rexp) const
{
    for (const auto& tval : _expansions)
    {
        if (tval.root() == rexp.root())
        {
            return true;
        }
    }
    
    return false; 
}

template <class T>
void
RecursionGroup<T>::merge(const RecursionGroup<T>& other)
{
    for (const auto& tval : other._expansions)
    {
        if (!contains(tval)) _expansions.push_back(tval);
    }
    
    std::sort(_expansions.begin(), _expansions.end());
}

template <class T>
size_t
RecursionGroup<T>::expansions() const
{
    return _expansions.size();
}

template <class T>
VRecursionTerms<T>
RecursionGroup<T>::roots() const
{
    VRecursionTerms<T> vrterms;
    
    for (const auto& tval : _expansions)
    {
        const auto rterm = tval.root();
        
        vrterms.push_back(RecursionTerm<T>(rterm.integral()));
    }
    
    return vrterms;
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
    
    size_t idx = 0;
    
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

template <class T>
bool
RecursionGroup<T>::empty() const
{
    for (const auto& tval : _expansions)
    {
        if (tval.terms() >  0)
        {
            return false;
        }
    }
    
    return true; 
}

template <class T>
bool
RecursionGroup<T>::auxilary(const int center) const
{
    for (const auto& tval : _expansions)
    {
        if (!tval.auxilary(center))
        {
            return false;
        }
    }
    
    return true;
}

template <class T>
template <class U>
std::optional<U>
RecursionGroup<T>::base() const
{
    std::set<U> tints;
    
    for (const auto& tval : _expansions)
    {
        const auto root = tval.root();
        
        tints.insert(U(root.integral()));
    }
    
    if (tints.size()  == 1)
    {
        return *(tints.begin());
    }
    else
    {
        return std::nullopt;
    }
}

template <class T>
int
RecursionGroup<T>::min_order() const
{
    int morder = INT_MAX;
    
    for (const auto& tval : _expansions)
    {
        if (const auto tord = tval.min_order(); tord < morder)
        {
            morder = tord;
        }
    }
    
    return morder;
}

template <class T>
Signature<T>
RecursionGroup<T>::signature() const
{
    Signature<T> tsign;
    
    const auto morder = min_order();
    
    for (const auto& tval : _expansions)
    {
        auto tint = (tval.root()).integral();
        
        if (const auto rint = tint.shift_order(-morder))
        {
            tsign.add(*rint, "out");
        }
        
        for (const auto& tint : tval.unique_integrals())
        {
            if (const auto rint = tint.shift_order(-morder))
            {
                tsign.add(*rint, "inp");
            }
        }
        
        for (const auto& fact : tval.unique_factors())
        {
            tsign.add(fact);
        }
    }
    
    return tsign;
}

template <class T>
std::set<Factor>
RecursionGroup<T>::factors() const
{
    std::set<Factor> sfacts;
    
    for (const auto& tval : _expansions)
    {
        if (const auto facts = tval.factors(); !facts.empty())
        {
            sfacts.insert(facts.cbegin(), facts.cend());
        }
    }
    
    return sfacts;
}

template <class T>
std::set<Fraction>
RecursionGroup<T>::prefactors() const
{
    std::set<Fraction> facts;
    
    for (const auto& tval : _expansions)
    {
        if (const auto tfacts = tval.prefactors(); !tfacts.empty())
        {
            facts.insert(tfacts.cbegin(), tfacts.cend());
        }
    }
    
    return facts; 
}

template <class T>
std::map<Factor, int>
RecursionGroup<T>::map_of_factors() const
{
    std::map<Factor, int> mfacts;
    
    for (const auto& tval : _expansions)
    {
        for (const auto& fact : tval.map_of_factors())
        {
            if (const auto tkval = mfacts.find(fact.first); tkval != mfacts.end())
            {
                mfacts[tkval->first] = tkval->second + fact.second;
            }
            else
            {
                mfacts[fact.first] = fact.second;
            }
        }
    }
    
    return mfacts;
}

template <class T>
std::set<T>
RecursionGroup<T>::components() const
{
    std::set<T> tcomps;
    
    for (const auto& tval : _expansions)
    {
        tcomps.insert((tval.root()).integral());
        
        for (size_t i = 0; i < tval.terms(); i++)
        {
            tcomps.insert(tval[i].integral());
        }
    }
    
    return tcomps;
}

template <class T>
template <class U>
std::set<T>
RecursionGroup<T>::components(const U& integral) const
{
    std::set<T> tcomps;
    
    for (const auto& tcomp : components())
    {
        if (U(tcomp) == integral)
        {
            tcomps.insert(tcomp);
        }
    }
    
    return tcomps;
}

template <class T>
template <class U>
std::set<U>
RecursionGroup<T>::integrals() const
{
    std::set<U> tints;
    
    for (const auto& tcomp : components())
    {
        tints.insert(U(tcomp));
    }
    
    return tints;
}

template <class T>
using VRecursionGroups = std::vector<RecursionGroup<T>>; 

#endif /* recursion_group_hpp */
