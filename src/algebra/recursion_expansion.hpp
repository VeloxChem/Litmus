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

#include "recursion_term.hpp"

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
    const RecursionTerm<T>& operator[](const size_t index) const;
    
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
    
    /// Creates new recursion expansion by splitting off recursion terms containing given integral.
    /// @param integral The other recursion expansion to compare.
    /// @return The new recursion expansion.
    RecursionExpansion<T> split(const T& integral) const;
    
    /// Simplifies recursion expansion.
    void simplify();
    
    /// Reduces order of recursion expansion by given order.
    /// @param order The order to be substracted from recursion expansion terms.
    void reduce(const int32_t order);
    
    /// Checks if this recursion expansion is similar to other recursion expansion.
    /// @param other The other recursion expansion to compare.
    /// @return True if recursion expansions  are similar, false otherwise.
    bool similar(const RecursionExpansion<T>& other) const;
    
    /// Adds recursion term to expansion of root recursion expansion.
    /// @param rterm The recursion term to add.
    void add(const RecursionTerm<T>& rterm);
    
    /// Gets root recursion term of this recursion expansion.
    /// @return The root recursion term.
    RecursionTerm<T> root() const;
    
    /// Gets number of recursion terms in expansion of root recursion term.
    /// @return The number of recursion terms in expansion of root recursion term.
    size_t terms() const;
    
    /// Gets unique integrals in recursion expansion.
    /// @return The set of unique integrals in recursion expansion.
    std::set<T> unique_integrals() const;
    
    /// Gets unique factors in recursion expansion.
    /// @return The set of unique factors in recursion expansion.
    std::set<Factor> unique_factors() const;
    
    /// Counts number of new integrals in recursion expansion with respect to given integral set.
    /// @param integrals The set of integrals used to count new integrals.
    /// @return The number of new integrals in recursion expansion.
    size_t count_new_integrals(const std::set<T>& integrals) const;
    
    /// Checks if recursion expansion contains auxilary root  at specific
    /// angular center.
    /// @param center The angular center to check.
    /// @return True if recursion expansion contains auxilary root, false otherwise.
    bool auxilary(const int center) const;
    
    /// Determines minimum order of integrals in this recursion expansion.
    /// @return The minimum order of integrlas.
    int min_order() const;
    
    /// Gets unique factors in this recursion expansion.
    /// @return The unique recursion factors.
    std::set<Factor> factors() const;
    
    /// Gets unique prefactors in this recursion expansion.
    /// @return The unique recursion prefactors.
    std::set<Fraction> prefactors() const;
    
    /// Gets map of factors in recursion expansion.
    /// @return The map of factors in recursion expansion.
    std::map<Factor, int> map_of_factors() const;
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
RecursionExpansion<T>::operator[](const size_t index) const
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
RecursionExpansion<T>
RecursionExpansion<T>::split(const T& integral) const
{
    VRecursionTerms<T> rterms;
    
    for (const auto& rterm : _expansion)
    {
        if (rterm.integral() == integral)
        {
            rterms.push_back(rterm);
        }
    }
    
    return RecursionExpansion<T>(_root, rterms);
}

template <class T>
void
RecursionExpansion<T>::simplify()
{
    SRecursionTerms<T> rterms;
    
    for (size_t i = 0; i < _expansion.size(); i++)
    {
        auto rterm = _expansion[i];
        
        rterm.prefactor(Fraction(0));
        
        rterms.insert(rterm);
    }
    
    VRecursionTerms<T> new_expansion;
    
    for (const auto& rterm : rterms)
    {
        auto frac = Fraction(0);
        
        for (size_t i = 0; i < _expansion.size(); i++)
        {
            if (rterm.same_base(_expansion[i]))
            {
                frac = frac + _expansion[i].prefactor();
            }
        }
        
        auto mterm = rterm;
        
        mterm.prefactor(frac);
        
        new_expansion.push_back(mterm);
    }
    
    _expansion = new_expansion; 
}

template <class T>
void
RecursionExpansion<T>::reduce(const int32_t order)
{
    _root = *(_root.shift_order(-order));
    
    for (size_t i = 0; i < _expansion.size(); i++)
    {
        _expansion[i] = *(_expansion[i].shift_order(-order)); 
    }
}

template <class T>
bool
RecursionExpansion<T>::similar(const RecursionExpansion<T>& other) const
{
    return _root.similar(other._root); 
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
size_t
RecursionExpansion<T>::terms() const
{
    return _expansion.size();
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
std::set<Factor>
RecursionExpansion<T>::unique_factors() const
{
    std::set<Factor> facts;
    
    for (const auto& rterm : _expansion)
    {
        for (const auto& fact : rterm.factors())
        {
            facts.insert(fact);
        }
    }

    return facts;
}

template <class T>
size_t
RecursionExpansion<T>::count_new_integrals(const std::set<T>& integrals) const
{
    size_t nints = 0;
    
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
int
RecursionExpansion<T>::min_order() const
{
    int morder = _root.order();
    
    for (const auto& rterm : _expansion)
    {
        if (const auto tord = rterm.order(); tord < morder)
        {
            morder = tord;
        }
    }
    
    return morder;
}

template <class T>
std::set<Factor>
RecursionExpansion<T>::factors() const
{
    std::set<Factor> sfacts;
    
    if (const auto facts = _root.factors(); !facts.empty())
    {
        sfacts.insert(facts.cbegin(), facts.cend());
    }
    
    for (const auto& rterm : _expansion)
    {
        if (const auto facts = rterm.factors(); !facts.empty())
        {
            sfacts.insert(facts.cbegin(), facts.cend());
        }
    }
    
    return sfacts;
}

template <class T>
std::set<Fraction>
RecursionExpansion<T>::prefactors() const
{
    std::set<Fraction> facts;
    
    for (const auto& rterm : _expansion)
    {
        facts.insert((rterm.prefactor()).abs());
    }
    
    return facts;
}

template <class T>
std::map<Factor, int>
RecursionExpansion<T>::map_of_factors() const
{
    std::map<Factor, int> mfacts(_root.map_of_factors());
    
    for (const auto& rterm : _expansion)
    {
        for (const auto& fact : rterm.map_of_factors())
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
using VRecursionExpansions = std::vector<RecursionExpansion<T>>; 

#endif /* recursion_expansion_hpp */
