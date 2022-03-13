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

#include "four_center_integral_component.hpp"

FourCenterIntegralComponent::FourCenterIntegralComponent()

    : _bra_pair(TwoCenterPairComponent())

    , _ket_pair(TwoCenterPairComponent())

    , _integrand(OperatorComponent())

    , _order(0)

    , _prefixes({})
{
    
}

FourCenterIntegralComponent::FourCenterIntegralComponent(const TwoCenterPairComponent& bra_pair,
                                                         const TwoCenterPairComponent& ket_pair,
                                                         const OperatorComponent&      integrand,
                                                         const int                     order,
                                                         const VOperatorComponents&    prefixes)

    : _bra_pair(bra_pair)

    , _ket_pair(ket_pair)

    , _integrand(integrand)

    , _order(order)

    , _prefixes(prefixes)
{
    
}

bool
FourCenterIntegralComponent::operator==(const FourCenterIntegralComponent& other) const
{
    if (this == &other) return true;

    if (_bra_pair != other._bra_pair)
    {
        return false;
    }
    else if (_ket_pair != other._ket_pair)
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

bool
FourCenterIntegralComponent::operator!=(const FourCenterIntegralComponent& other) const
{
    return !((*this) == other);
}

bool
FourCenterIntegralComponent::operator<(const FourCenterIntegralComponent& other) const
{
    if (_bra_pair != other._bra_pair)
    {
        return _bra_pair < other._bra_pair;
    }
    else if (_ket_pair != other._ket_pair)
    {
        return _ket_pair < other._ket_pair;
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

std::string
FourCenterIntegralComponent::to_string() const
{
    std::string intstr;
    
    if (!_prefixes.empty())
    {
        intstr.append("[");
        
        for (const auto& prefix : _prefixes)
        {
            intstr.append(prefix.to_string() + ";");
        }
        
        intstr.append("]");
    }
    
    intstr.append(_bra_pair.to_string());
    
    intstr.append(_integrand.to_string());
    
    intstr.append(_ket_pair.to_string());
    
    intstr.append("^(" + std::to_string(_order) + ")");

    return intstr;
}

std::string
FourCenterIntegralComponent::label(const bool use_order) const
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
