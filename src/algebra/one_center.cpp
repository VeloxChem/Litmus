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

OneCenter::OneCenter(const OneCenterComponent& tcomp)

    : _name(tcomp.name())

    , _shape(Tensor((tcomp.shape()).order()))
{

}

int
OneCenter::operator[](const int index) const
{
    return _shape.order();
}

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

std::optional<OneCenter>
OneCenter::shift(const int value,
                 const int center) const
{
    if (center == 0)
    {
        if (const int torder =_shape.order() + value; torder >= 0)
        {
            return OneCenter(_name, Tensor(torder));
        }
        else
        {
            return std::nullopt;
        }
    }
    else
    {
        return std::nullopt;
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

VOneCenterComponents
OneCenter::components() const
{
    VOneCenterComponents tcomps;
    
    for (const auto& tcomp : _shape.components())
    {
        tcomps.push_back(OneCenterComponent(_name, tcomp));
    }
    
    return tcomps;
}
