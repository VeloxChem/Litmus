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

#ifndef four_center_integral_hpp
#define four_center_integral_hpp

#include <vector>

#include "operator.hpp"
#include "integral_component.hpp"
#include "components.hpp"

/// Integral class.
template <class T, class U>
class Integral
{
    /// Bra side of integral.
    T _bra;
    
    /// Ket side of integral.
    U _ket;
    
    /// Integrand operator of integral.
    Operator _integrand;
    
    /// Order of integral.
    int _order;
    
    ///  The vector of prefix operators acting on integral.
    VOperators _prefixes;
    
public:
    /// Creates an empty integral.
    Integral();
    
    /// Creates a integral from the given operator, bra and ket sides.
    /// @param bra The bra side of integral.
    /// @param ket The ket side of integral.
    /// @param integrand The integrand operator of integral.
    /// @param order The order of integral.
    /// @param prefixes The vector of prefix operators acting on integral.
    Integral(const T&          bra,
             const U&          ket,
             const Operator&   integrand,
             const int         order = 0,
             const VOperators& prefixes = VOperators({}));
    
    /// Creates a integral from the given integral component.
    /// @param tcomp The integral component to create integral.
    template <class V, class W>
    Integral(const IntegralComponent<V, W>& tcomp);
    
    /// Retrieves angular momentum of requested center in this integral.
    /// @param center The index of center to retrieve angular momentum.
    /// @return The   angular momentum of requested center in.
    int operator[](const int center) const;
    
    /// Compares this integral with other integral.
    /// @param other The other integral to compare.
    /// @return true if integrals, false otherwise.
    bool operator==(const Integral<T, U>& other) const;
    
    /// Compares this integral with other integral.
    /// @param other The other integral to compare.
    /// @return true if integrals are not equal, false otherwise.
    bool operator!=(const Integral<T, U>& other) const;
    
    /// Compares this integral with other integral.
    /// @param other The other integral to compare.
    /// @return true if this integral is less than other integral, false otherwise.
    bool operator<(const Integral<T, U>& other) const;
    
    /// Sets order of integral.
    /// @param order The order of integral.
    void set_order(const int order);
    
    /// Gets operator of this integral.
    /// @return The operator of this integral.
    Operator integrand() const;
    
    /// Gets order of integral.
    /// @return The order of integral.
    int32_t order() const;
    
    /// Checks if integrals are without prefix operators.
    /// @return True if integral is simple, False otherwise.
    bool is_simple() const;
    
    /// Checks if integral has simple integrand.
    /// @return True if integrand is simple, False otherwise.
    bool is_simple_integrand() const;
    
    /// Creates primitive textual label of this integral.
    /// @return The string with primitive textual label of integral.
    std::string label(const bool use_order = false) const;
    
    /// Creates a vector with integral components of this integral.
    /// @return The vector of integral components.
    template <class V, class W>
    VIntegralComponents<V, W> components() const;
    
    /// Creates a vector with diagonal integral components of this integral.
    /// @return The vector of integral components.
    template <class V, class W>
    VIntegralComponents<V, W> diag_components() const;
};

template <class T, class U>
Integral<T, U>::Integral()

    : _bra(T())

    , _ket(U())

    , _integrand(Operator())

    , _order(0)

    , _prefixes({})
{
    
}

template <class T, class U>
Integral<T, U>::Integral(const T&          bra,
                         const U&          ket,
                         const Operator&   integrand,
                         const int         order,
                         const VOperators& prefixes)

    : _bra(bra)

    , _ket(ket)

    , _integrand(integrand)

    , _order(order)

    , _prefixes(prefixes)
{
    
}

template <class T, class U>
template <class V, class W>
Integral<T, U>::Integral(const IntegralComponent<V, W>& tcomp)

    : _bra(tcomp.bra())

    , _ket(tcomp.ket())

    , _integrand(tcomp.integrand())

    , _order(tcomp.order())

    , _prefixes({})
{
    for (const auto& opcomp : tcomp.prefixes())
    {
        _prefixes.push_back(Operator(opcomp));
    }
}

template <class T, class U>
int
Integral<T, U>::operator[](const int center) const
{
    if (const auto bcenters = _bra.centers(); center < bcenters)
    {
        return _bra[center];
    }
    else
    {
        return _ket[center - bcenters];
    }
}

template <class T, class U>
bool
Integral<T, U>::operator==(const Integral<T, U>& other) const
{
    if (this == &other) return true;

    if (_bra != other._bra)
    {
        return false;
    }
    else if (_ket != other._ket)
    {
        return false;
    }
    else if (_integrand != other._integrand)
    {
        return false;
    }
    if (_order != other._order)
    {
        return false;
    }
    else
    {
        return _prefixes == other._prefixes;
    }
}

template <class T, class U>
bool
Integral<T, U>::operator!=(const Integral<T, U>& other) const
{
    return !((*this) == other);
}

template <class T, class U>
bool
Integral<T, U>::operator<(const Integral<T, U>& other) const
{
    if (_bra != other._bra)
    {
        return _bra < other._bra;
    }
    else if (_ket != other._ket)
    {
        return _ket < other._ket;
    }
    else if (_integrand != other._integrand)
    {
        return _integrand < other._integrand;
    }
    if (_order != other._order)
    {
        return _order < other._order;
    }
    else
    {
        return _prefixes < other._prefixes;
    }
}

template <class T, class U>
void
Integral<T, U>::set_order(const int order)
{
    _order = order; 
}

template <class T, class U>
Operator
Integral<T, U>::integrand() const
{
    return _integrand;
}

template <class T, class U>
int32_t
Integral<T, U>::order() const
{
    return _order;
}

template <class T, class U>
bool
Integral<T, U>::is_simple() const
{
    return _prefixes.empty();
}

template <class T, class U>
bool
Integral<T, U>::is_simple_integrand() const
{
    return _integrand.components().size() == 1;
}

template <class T, class U>
std::string
Integral<T, U>::label(const bool use_order) const
{
    std::string intstr;
    
    intstr.append(_bra.label());
    
    intstr.append(_ket.label());
    
    if (use_order)
    {
        intstr.append("_" + std::to_string(_order));
    }

    return intstr;
}

template <class T, class U>
template <class V, class W>
VIntegralComponents<V, W>
Integral<T, U>::components() const
{
    VIntegralComponents<V, W> vcomps;
    
    if (!_prefixes.empty())
    {
        for (const auto& prefix : make_components<OperatorComponent>(_prefixes))
        {
            for (const auto& opcomp : _integrand.components())
            {
                for (const auto& bra : _bra.components())
                {
                    for (const auto& ket : _ket.components())
                    {
                        vcomps.push_back(IntegralComponent<V, W>(bra, ket, opcomp, _order, prefix));
                    }
                }
            }
        }
    }
    else
    {
        for (const auto& opcomp : _integrand.components())
        {
            for (const auto& bra : _bra.components())
            {
                for (const auto& ket : _ket.components())
                {
                    vcomps.push_back(IntegralComponent<V, W>(bra, ket, opcomp, _order));
                }
            }
        }
    }
    
    return vcomps;
}

template <class T, class U>
template <class V, class W>
VIntegralComponents<V, W>
Integral<T, U>::diag_components() const
{
    VIntegralComponents<V, W> vcomps;
    
    if (!_prefixes.empty())
    {
        for (const auto& prefix : make_components<OperatorComponent>(_prefixes))
        {
            for (const auto& opcomp : _integrand.components())
            {
                const auto bcomps = _bra.components();
                
                const auto kcomps = _ket.components();
                
                for (size_t i = 0; i < bcomps.size(); i++)
                {
                    vcomps.push_back(IntegralComponent<V, W>(bcomps[i], kcomps[i], opcomp, _order, prefix));
                }
            }
        }
    }
    else
    {
        for (const auto& opcomp : _integrand.components())
        {
            const auto bcomps = _bra.components();
            
            const auto kcomps = _ket.components();
            
            for (size_t i = 0; i < bcomps.size(); i++)
            {
                vcomps.push_back(IntegralComponent<V, W>(bcomps[i], kcomps[i], opcomp, _order));
            }
        }
    }
    
    return vcomps;
}



#endif /* four_center_integral_hpp */
