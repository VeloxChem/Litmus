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

#ifndef recursion_term_hpp
#define recursion_term_hpp

#include <map>
#include <string>
#include <set>
#include <vector>

#include "factor.hpp"
#include "fraction.hpp"
#include "operator_component.hpp"

/// Recursion term class.
template <class T>
class RecursionTerm
{
    /// Integral component of  recursion term.
    T _integral;
    
    /// Map of factors of four center recursion term.
    std::map<Factor, int> _factors;
  
    /// Scalar fractional prefactor of four center recursion term.
    Fraction _prefactor;
    
public:
    /// Creates an empty recursion term.
    RecursionTerm();
    
    /// Creates a  recursion term from the given integral component, map of factors, and scalar prefactor.
    /// @param integral The integral component of recursion term.
    /// @param factors The map of factors of recursion term.
    /// @param prefactor The scalar fractional prefactor of recursion term.
    RecursionTerm(const T&                     integral,
                  const std::map<Factor, int>& factors = std::map<Factor, int>({}),
                  const Fraction&              prefactor = Fraction(1));
    
    /// Retrieves axial value along requested center of integral in recursion term .
    /// @param center The index of center to retrieve axial value.
    /// @return The  tensorial shape of requested center of integral component in recursion term .
    const TensorComponent& operator[](const int center) const;
    
    /// Compares this recursion term with other recursion term.
    /// @param other The other recursion term to compare.
    /// @return true if recursion terms are equal, false otherwise.
    bool operator==(const RecursionTerm<T>& other) const;
    
    /// Compares this recursion term with other recursion term.
    /// @param other The other recursion term to compare.
    /// @return true if recursion terms are not equal, false otherwise.
    bool operator!=(const RecursionTerm<T>& other) const;
    
    /// Compares this recursion term with other recursion term.
    /// @param other The other recursion term to compare.
    /// @return true if this recursion term is less than other recursion term, false otherwise.
    bool operator<(const RecursionTerm<T>& other) const;
    
    /// Sets fractional prefactor of recursion term.
    /// @param factor The fractional prefactor of recursion term.
    void prefactor(const Fraction& factor);
    
    /// Checks if this recursion term is similar to other recursion term.
    /// @param other The other recursion term to compare.
    /// @return True if recursion terms  are similar, false otherwise.
    bool similar(const RecursionTerm<T>& other) const;
    
    /// Checks if this recursion term has same base as the other recursion term.
    /// @param other The other recursion term to compare.
    /// @return True if recursion terms  have same base, false otherwise.
    bool same_base(const RecursionTerm<T>& other) const;
    
    /// Gets bra side of recursion term.
    /// @return The bra side of recursion term.
    template <class U>
    U bra() const;
    
    /// Gets ket side of recursion term.
    /// @return The ket side of recursion term.
    template <class U>
    U ket() const;
    
    /// Gets integrand of recursion term.
    /// @return The integrand of recursion term.
    OperatorComponent integrand() const;

    /// Gets order of recursion term.
    /// @return The order of recursion term.
    int order() const;
    
    /// Gets vector of prefix operator components of recursion term.
    /// @return The vector of prefix operator components of recursion term.
    VOperatorComponents prefixes() const;
    
    /// Gets vector of prefix operator orders.
    /// @return The vector of prefix operator orders.
    std::vector<int> prefixes_order() const;
    
    /// Gets integral of recursion term.
    /// @return The integral of recursion term.
    T integral() const;

    /// Gets fractional prefactor of recursion term.
    /// @return The fractional prefactor of recursion term.
    Fraction prefactor() const;
    
    /// Gets set of factors in recursion term.
    /// @return The set of factors in recursion term.
    std::set<Factor> factors() const;
    
    /// Gets map of factors in recursion term.
    /// @return The map of factors in recursion term.
    std::map<Factor, int> map_of_factors() const;
    
    /// Gets order of specific factors in recursion term.
    /// @return The order of specific factors in recursion term.
    int factor_order(const Factor& factor) const;
    
    /// Creates primitive textual label of this recursion term.
    /// @param use_order The flag to include order of integral into its label.
    /// @return The string with primitive textual label of recursion term.
    std::string label(const bool use_order = false) const;
    
    /// Creates an term from this recursion term by replacing integrand with given
    /// operator component.
    /// @param integrand The integrand operator component of recursion term.
    /// @return The created recursion term.
    RecursionTerm replace(const OperatorComponent& integrand) const;
    
    /// Removes recursion factor with given name for recursion term.
    /// @param name The name of factor to be removed.
    /// @return The created recursion term.
    RecursionTerm remove(const std::string& name) const;
    
    /// Creates an optional recursion term from this recursion term by shifting axial value
    /// along the selected axis on targeted center.
    /// @param axis The axis to shift axial value.
    /// @param value The value to shift axial value.
    /// @param center The targeted shift center.
    /// @return The optional recursion term.
    std::optional<RecursionTerm> shift(const char axis,
                                       const int  value,
                                       const int  center) const;
    
    /// Creates an optional recursion term from this recursion term by shifting axial value
    /// along the selected axis in selected prefix operator.
    /// @param axis The axis to shift axial value.
    /// @param value The value to shift axial value.
    /// @param index The index of targeted operator.
    /// @param noscalar The flag for scalar component: false  to keep,  true otherwise.
    /// @return The optional recursion term.
    std::optional<RecursionTerm> shift_prefix(const char axis,
                                              const int  value,
                                              const int  index,
                                              const bool noscalar = false) const;
    
    /// Creates an optional recursion term from this recursion term by shifting order of
    /// recursion term..
    /// @param value The value to shift order of recursion term.
    /// @return The optional recursion term.
    std::optional<RecursionTerm> shift_order(const int value) const;
    
    /// Creates an optional recursion term from this recursion term by shifting axial value
    /// along the selected axis on integrand.
    /// @param axis The axis to shift axial value.
    /// @param value The value to shift axial value.
    /// @return The optional recursion term.
    std::optional<RecursionTerm> shift_operator(const char axis,
                                                const int  value) const;
    
    /// Adds new factor or updates order of existing factor in this recursion term.
    /// @param factor The factor to scale recurion term.
    /// @param multiplier The fractional multiplier accompanying factor to add.
    void add(const Factor&   factor,
             const Fraction& multiplier = Fraction(1));

    void add(const Fraction& multiplier = Fraction(1));
    
    /// Scales prefactor of this recursion term with the given factor.
    /// @param factor The fractional factor to scale recurion term.
    void scale(const Fraction& factor);
    
    /// Checks if recursion term is auxilary at specific
    /// angular center.
    /// @param center The angular center to check.
    /// @return True if recursion term is auxilary, false otherwise.
    bool auxilary(const int center) const;
    
    /// Removes prefixes from recursion term..
    void clear_prefixes();
};

template <class T>
RecursionTerm<T>::RecursionTerm()
    
    : _integral(T())

    , _factors(std::map<Factor, int>({}))

    , _prefactor(Fraction(1))
{
    
}

template <class T>
RecursionTerm<T>::RecursionTerm(const T&                     integral,
                                const std::map<Factor, int>& factors,
                                const Fraction&              prefactor)

    : _integral(integral)

    , _factors(factors)

    , _prefactor(prefactor)
{
    
}

template <class T>
const TensorComponent&
RecursionTerm<T>::operator[](const int center) const
{
    return _integral[center];
}

template <class T>
bool
RecursionTerm<T>::operator==(const RecursionTerm<T>& other) const
{
    if (this == &other) return true;

    if (_integral != other._integral)
    {
        return false;
    }
    else if (_factors != other._factors)
    {
        return false;
    }
    else
    {
        return _prefactor == other._prefactor;
    }
}

template <class T>
bool
RecursionTerm<T>::operator!=(const RecursionTerm<T>& other) const
{
    return !((*this) == other);
}

template <class T>
bool
RecursionTerm<T>::operator<(const RecursionTerm<T>& other) const
{
    if (_integral != other._integral)
    {
        return _integral < other._integral;
    }
    else if (_factors != other._factors)
    {
        return _factors < other._factors;
    }
    else
    {
        return _prefactor < other._prefactor;
    }
}

template <class T>
void
RecursionTerm<T>::prefactor(const Fraction& factor)
{
    _prefactor = factor; 
}

template <class T>
bool
RecursionTerm<T>::similar(const RecursionTerm<T>& other) const
{
    return _integral.similar(other._integral); 
}

template <class T>
bool
RecursionTerm<T>::same_base(const RecursionTerm<T>& other) const
{
    if (_integral != other._integral)
    {
        return false;
    }
    
    if (_factors != other._factors)
    {
        return false;
    }
    
    return true;
}

template <class T>
template <class U>
U
RecursionTerm<T>::bra() const
{
    return _integral.bra();
}

template <class T>
template <class U>
U
RecursionTerm<T>::ket() const
{
    return _integral.ket();
}

template <class T>
OperatorComponent
RecursionTerm<T>::integrand() const
{
    return _integral.integrand();
}

template <class T>
int
RecursionTerm<T>::order() const
{
    return _integral.order();
}

template <class T>
VOperatorComponents
RecursionTerm<T>::prefixes() const
{
    return _integral.prefixes();
}

template <class T>
std::vector<int>
RecursionTerm<T>::prefixes_order() const
{
    return _integral.prefixes_order();
}

template <class T>
T
RecursionTerm<T>::integral() const
{
    return _integral;
}

template <class T>
Fraction
RecursionTerm<T>::prefactor() const
{
    return _prefactor;
}

template <class T>
std::set<Factor>
RecursionTerm<T>::factors() const
{
    std::set<Factor> facts;
    
    for (const auto& tkval : _factors)
    {
        facts.insert(tkval.first);
    }
    
    return facts; 
}

template <class T>
std::map<Factor, int>
RecursionTerm<T>::map_of_factors() const
{
    return _factors;
}

template <class T>
int
RecursionTerm<T>::factor_order(const Factor& factor) const
{
    if (const auto tkval = _factors.find(factor); tkval != _factors.cend())
    {
        return tkval->second;
    }
    else
    {
        return 0;
    }
}

template <class T>
std::string
RecursionTerm<T>::label(const bool use_order) const
{
    return _integral.label(use_order);
}

template <class T>
RecursionTerm<T>
RecursionTerm<T>::replace(const OperatorComponent& integrand) const
{
    const auto tint = _integral.replace(integrand);
    
    return RecursionTerm<T>(tint, _factors, _prefactor);
}

template <class T>
RecursionTerm<T>
RecursionTerm<T>::remove(const std::string& name) const
{
    std::map<Factor, int> facts;
    
    for (const auto& [fact, repval] : _factors)
    {
        if (fact.name() != name)
        {
            facts.insert({fact, repval});
        }
    }
    
    return RecursionTerm<T>(_integral, facts, _prefactor);
}

template <class T>
std::optional<RecursionTerm<T>>
RecursionTerm<T>::shift(const char axis,
                        const int  value,
                        const int  center) const
{
    if (const auto tint = _integral.shift(axis, value, center))
    {
        return RecursionTerm<T>(*tint, _factors, _prefactor);
    }
    else
    {
        return std::nullopt;
    }
}

template <class T>
std::optional<RecursionTerm<T>>
RecursionTerm<T>::shift_prefix(const char axis,
                               const int  value,
                               const int  index,
                               const bool noscalar) const
{
    if (const auto tint = _integral.shift_prefix(axis, value, index, noscalar))
    {
        return RecursionTerm<T>(*tint, _factors, _prefactor);
    }
    else
    {
        return std::nullopt;
    }
}

template <class T>
std::optional<RecursionTerm<T>>
RecursionTerm<T>::shift_order(const int value) const
{
    if (const auto tint = _integral.shift_order(value))
    {
        return RecursionTerm<T>(*tint, _factors, _prefactor);
    }
    else
    {
        return std::nullopt;
    }
}

template <class T>
std::optional<RecursionTerm<T>>
RecursionTerm<T>::shift_operator(const char axis,
                                 const int  value) const
{
    if (const auto tint = _integral.shift_operator(axis, value))
    {
        return RecursionTerm<T>(*tint, _factors, _prefactor);
    }
    else
    {
        return std::nullopt;
    }
}

template <class T>
void
RecursionTerm<T>::add(const Factor&   factor,
                      const Fraction& multiplier)
{
    if (const auto tkval = _factors.find(factor); tkval != _factors.end())
    {
        _factors[tkval->first] = tkval->second + 1;
    }
    else
    {
        _factors[factor] = 1;
    }
    
    _prefactor = _prefactor * multiplier;
}

template <class T>
void
RecursionTerm<T>::add(const Fraction& multiplier)
{
    _prefactor = _prefactor * multiplier;
}

template <class T>
void
RecursionTerm<T>::scale(const Fraction& factor)
{
    _prefactor = _prefactor * factor;
}

template <class T>
bool
RecursionTerm<T>::auxilary(const int center) const
{
    if (!_integral.prefixes().empty())
    {
        return false; 
    }
    
    if (_integral[center].order() == 0)
    {
        return true;
    }
    else
    {
        return false; 
    }
}

template <class T>
void
RecursionTerm<T>::clear_prefixes()
{
    _integral = _integral.base(); 
}

template <class T>
using VRecursionTerms = std::vector<RecursionTerm<T>>;

template <class T>
using SRecursionTerms = std::set<RecursionTerm<T>>;

template <class T>
using MRecursionTerms = std::vector<std::vector<RecursionTerm<T>>>;

#endif /* recursion_term_hpp */
