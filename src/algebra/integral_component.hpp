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

#ifndef four_center_integral_component_hpp
#define four_center_integral_component_hpp

#include "operator_component.hpp"

/// Integral component.
template <class T, class U>
class IntegralComponent
{
    /// Bra side of integral component.
    T _bra;
    
    /// Ket side of integral component.
    U _ket;
    
    /// Integrand operator of integral component.
    OperatorComponent _integrand;
    
    /// Order of integral component.
    int _order;
    
    ///  The vector of prefix operators acting on integral component.
    VOperatorComponents _prefixes;
    
public:
    /// Creates an empty  component.
    IntegralComponent();
    
    /// Creates a integral component from the given operator component, bra and ket sides.
    /// @param bra The bra side of integral component.
    /// @param ket The ket side of integral component.
    /// @param integrand The integrand operator component of integral component.
    /// @param order The order of integral component.
    /// @param prefixes The vector of prefix operator components acting on integral component.
    IntegralComponent(const T&                   bra,
                      const U&                   ket,
                      const OperatorComponent&   integrand,
                      const int                  order = 0,
                      const VOperatorComponents& prefixes = VOperatorComponents({}));
    
    /// Retrieves axial value along requested center of integral.
    /// @param center The index of center to retrieve axial value.
    /// @return The  tensorial shape of requested center in integral component.
    const TensorComponent& operator[](const int center) const;
    
    /// Compares this integral component with other integral component.
    /// @param other The other integral component to compare.
    /// @return true if integral components, false otherwise.
    bool operator==(const IntegralComponent& other) const;
    
    /// Compares this integral component with other integral component.
    /// @param other The other integral  component to compare.
    /// @return true if integral components are not equal, false otherwise.
    bool operator!=(const IntegralComponent& other) const;
    
    /// Compares this integral component with other integral component.
    /// @param other The other integral component to compare.
    /// @return true if this integral  component is less than other integral component, false otherwise.
    bool operator<(constIntegralComponent& other) const;
    
    /// Gets bra side of integral component.
    /// @return The bra side of integral component.
    T bra() const;
    
    /// Gets ket side of integral component.
    /// @return The ket side of integral component.
    U ket() const;
    
    /// Gets integrand of integral component.
    /// @return The integrand of integral component.
    OperatorComponent integrand() const;

    /// Gets order of integral component.
    /// @return The order of integral component.
    int order() const;
    
    /// Gets vector of prefix operator components  of integral component.
    /// @return The vector of prefix operator components of integral component.
    VOperatorComponents prefixes() const;
    
    /// Creates primitive textual label of this integral component.
    /// @param use_order The flag to include order of integral into its label.
    /// @return The string with primitive textual label of integral component.
    std::string label(const bool use_order = false) const;
    
    /// Creates an  integral component from this integral component by replacing
    /// integrand with given operator component.
    /// @param integrand The integrand operator component of integral component.
    /// @return The created integral component.
    IntegralComponent replace(const OperatorComponent& integrand) const;
    
    /// Creates an optional integral component from this integral component by shifting axial value
    /// along the selected axis on targeted center.
    /// @param axis The axis to shift axial value.
    /// @param value The value to shift axial value.
    /// @param center The targeted shift center.
    /// @return The optional integral component.
    std::optional<IntegralComponent> shift(const char axis,
                                           const int  value,
                                           const int  center) const;
    
    /// Creates an optional integral component from this integral component by shifting axial value
    /// along the selected axis in selected prefix operator.
    /// @param axis The axis to shift axial value.
    /// @param value The value to shift axial value.
    /// @param index The index of targeted operator.
    /// @param noscalar The flag for scalar component: false  to keep,  true otherwise.
    /// @return The optional operator component.
    std::optional<IntegralComponent> shift_prefix(const char axis,
                                                  const int  value,
                                                  const int  index,
                                                  const bool noscalar = false) const;
};

template <class T, class U>
IntegralComponent::IntegralComponent()

    : _bra(T())

    , _ket(U())

    , _integrand(OperatorComponent())

    , _order(0)

    , _prefixes({})
{
    
}

template <class T, class U>
IntegralComponent::IntegralComponent(const T&                   bra,
                                     const U&                   ket,
                                     const OperatorComponent&   integrand,
                                     const int                  order,
                                     const VOperatorComponents& prefixes)

    : _bra(bra)

    , _ket(ket)

    , _integrand(integrand)

    , _order(order)

    , _prefixes(prefixes)
{
    
}

template <class T, class U>
const TensorComponent&
IntegralComponent::operator[](const int center) const
{
    const auto bcenters = _bra.centers();
    
    if (center < bcenters)
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
IntegralComponent::operator==(const IntegralComponent& other) const
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
IntegralComponent::operator!=(const IntegralComponent& other) const
{
    return !((*this) == other);
}

template <class T, class U>
bool
IntegralComponent::operator<(const IntegralComponent& other) const
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
T
IntegralComponent::bra() const
{
    return _bra;
}

template <class T, class U>
U
IntegralComponent::ket() const
{
    return _ket;
}

template <class T, class U>
OperatorComponent
IntegralComponent::integrand() const
{
    return _integrand;
}

template <class T, class U>
int
IntegralComponent::order() const
{
    return _order;
}

template <class T, class U>
VOperatorComponents
IntegralComponent::prefixes() const
{
    return _prefixes;
}

template <class T, class U>
std::string
IntegralComponent::label(const bool use_order) const
{
    std::string intstr;
    
    if (!_prefixes.empty())
    {
        for (const auto& prefix : _prefixes)
        {
            intstr.append(prefix.label() + "_");
        }
    }
    
    if (const auto lblstr = _integrand.label(); lblstr != "0")
    {
        intstr.append(lblstr + "_");
    }
    
    intstr.append(_bra_pair.label() + "_");
    
    intstr.append(_ket_pair.label());
    
    if (use_order)
    {
        intstr.append("_" + std::to_string(_order));
    }

    return intstr;
}

template <class T, class U>
IntegralComponent<T, U>
IntegralComponent::replace(const OperatorComponent& integrand) const
{
    return IntegralComponent<T, U>(_bra, _ket, integrand, _order, _prefixes);
}


template <class T, class U>
std::optional<IntegralComponent<T, U>>
IntegralComponent::shift(const char axis,
                         const int  value,
                         const int  center) const
{
    const auto bcenters = _bra.centers();
    
    if (center < bcenters)
    {
        if (const auto tbra = _bra.shift(axis, value, center))
        {
            return IntegralComponent<T, U>(*tbra, _ket, _integrand, _order, _prefixes);
        }
        else
        {
            return std::nullopt;
        }
    }
    else
    {
        if (const auto tket = _ket.shift(axis, value, center - bcenters))
        {
            return IntegralComponent<T, U>(_bra, *tket, _integrand, _order, _prefixes);
        }
        else
        {
            return std::nullopt;
        }
    }
}

template <class T, class U>
std::optional<IntegralComponent<T, U>>
IntegralComponent::shift_prefix(const char axis,
                                const int  value,
                                const int  index,
                                const bool noscalar) const
{
    if (const auto opcomp = _prefixes[index].shift(axis, value, noscalar))
    {
        auto new_prefixes = _prefixes;
        
        new_prefixes[index] = *opcomp;
        
        return IntegralComponent<T, U>(_bra, _ket, _integrand, _order, new_prefixes);
    }
    else
    {
        if (const auto opcomp = _prefixes[index].shift(axis, value))
        {
            auto new_prefixes = _prefixes;
            
            new_prefixes.erase(new_prefixes.begin() + index);
            
            return IntegralComponent<T, U>(_bra, _ket, _integrand, _order, new_prefixes);
        }
        
        return std::nullopt;
    }
}

template <class T, class U>
using VIntegralComponents = std::vector<IntegralComponent<T, U>>;

#endif /* four_center_integral_component_hpp */
