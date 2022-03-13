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

#include "four_center_integral.hpp"

FourCenterIntegral::FourCenterIntegral()

    : _bra_pair(TwoCenterPair())

    , _ket_pair(TwoCenterPair())

    , _integrand(Operator())

    , _order(0)

    , _prefixes({})
{
    
}

FourCenterIntegral::FourCenterIntegral(const TwoCenterPair& bra_pair,
                                       const TwoCenterPair& ket_pair,
                                       const Operator&      integrand,
                                       const int            order,
                                       const VOperators&    prefixes)

    : _bra_pair(bra_pair)

    , _ket_pair(ket_pair)

    , _integrand(integrand)

    , _order(order)

    , _prefixes(prefixes)
{
    
}

FourCenterIntegral::FourCenterIntegral(const int          a_angmom,
                                       const int          b_angmom,
                                       const int          c_angmom,
                                       const int          d_angmom,
                                       const Operator&    integrand,
                                       const int          order,
                                       const VOperators&  prefixes)

    : _bra_pair(TwoCenterPair("GA", a_angmom, "GB", b_angmom))

    , _ket_pair(TwoCenterPair("GC", c_angmom, "GD", d_angmom))

    , _integrand(integrand)

    , _order(order)

    , _prefixes(prefixes)
{
    
}

bool
FourCenterIntegral::operator==(const FourCenterIntegral& other) const
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
FourCenterIntegral::operator!=(const FourCenterIntegral& other) const
{
    return !((*this) == other);
}

bool
FourCenterIntegral::operator<(const FourCenterIntegral& other) const
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
FourCenterIntegral::label(const bool use_order) const
{
    std::string intstr;
    
    intstr.append(_bra_pair.label());
    
    intstr.append(_ket_pair.label());
    
    if (use_order)
    {
        intstr.append("_" + std::to_string(_order));
    }

    return intstr;
}

