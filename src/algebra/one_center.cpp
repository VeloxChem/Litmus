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

#include "one_center.hpp"

OneCenter::OneCenter()

    : _name(std::string())

    , _shape(Tensor(0))
{
    
}

OneCenter::OneCenter(const std::string& name,
                     const Tensor&      shape)
    
    : _name(name)

    , _shape(shape)
{
    
}

OneCenter::OneCenter(const std::string& name,
                     const int          angmom)
    : _name(name)

    , _shape(Tensor(angmom))
{
    
}

//TwoCenterPair::TwoCenterPair(const TwoCenterPairComponent& t2pcomp)
//
//    : _names(t2pcomp.names())
//
//    , _shapes(std::array<Tensor, 2>({Tensor(0), Tensor(0)}))
//{
//    const auto tcomps = t2pcomp.shapes();
//
//    _shapes[0] = Tensor(tcomps[0].order());
//
//    _shapes[1] = Tensor(tcomps[1].order());
//}

bool
OneCenter::operator==(const OneCenter& other) const
{
    if (this == &other) return true;

    if (_name != other._name)
    {
        return false;
    }
    else
    {
        return _shape == other._shape;
    }
}

bool
OneCenter::operator!=(const OneCenter& other) const
{
    return !((*this) == other);
}

bool
OneCenter::operator<(const OneCenter& other) const
{
    if (_name != other._name)
    {
        return _name < other._name;
    }
    else
    {
        return _shape < other._shape;
    }
}

std::string
OneCenter::to_string() const
{
    return "{" + _name + ":" + _shape.to_string() + "}";
}

std::string
OneCenter::label() const
{
    return _shape.label();
}

//VTwoCenterPairComponents
//TwoCenterPair::components() const
//{
//    VTwoCenterPairComponents t2pcomps;
//    
//    for (const auto& fcomp : _shapes[0].components())
//    {
//        for (const auto& scomp : _shapes[1].components())
//        {
//            t2pcomps.push_back(TwoCenterPairComponent(_names, {fcomp, scomp}));
//        }
//    }
//    
//    return t2pcomps;
//}
