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

#include "two_center_pair_component.hpp"

TwoCenterPairComponent::TwoCenterPairComponent()

    : _names(std::array<std::string, 2>({"", ""}))

    , _shapes(std::array<TensorComponent, 2>({TensorComponent(0, 0, 0), TensorComponent(0, 0, 0)}))
{
    
}

TwoCenterPairComponent::TwoCenterPairComponent(const std::array<std::string, 2>&     names,
                                               const std::array<TensorComponent, 2>& shapes)
    
    : _names(names)

    , _shapes(shapes)
{
    
}

const TensorComponent&
TwoCenterPairComponent::operator[](const int index) const
{
    return _shapes[index];
}

bool
TwoCenterPairComponent::operator==(const TwoCenterPairComponent& other) const
{
    if (this == &other) return true;

    if (_names != other._names)
    {
        return false;
    }
    else
    {
        return _shapes == other._shapes;
    }
}

bool
TwoCenterPairComponent::operator!=(const TwoCenterPairComponent& other) const
{
    return !((*this) == other);
}

bool
TwoCenterPairComponent::operator<(const TwoCenterPairComponent& other) const
{
    if (_names != other._names)
    {
        return _names < other._names;
    }
    else 
    {
        return _shapes < other._shapes;
    }
}

bool
TwoCenterPairComponent::similar(const TwoCenterPairComponent& other) const
{
    if (_names != other._names)
    {
        return false;
    }
    else if (!_shapes[0].similar(other._shapes[0]))
    {
        return false;
    }
    else
    {
        return _shapes[1].similar(other._shapes[1]);
    }
}

std::array<std::string, 2>
TwoCenterPairComponent::names() const
{
    return _names;
}

std::array<TensorComponent, 2>
TwoCenterPairComponent::shapes() const
{
    return _shapes;
}

std::string
TwoCenterPairComponent::to_string() const
{
    return "{" + _names[0] + ":" + _shapes[0].to_string() + ";" +
    
           _names[1] + ":"  + _shapes[1].to_string() + "}";
}

std::string
TwoCenterPairComponent::label() const
{
    return _shapes[0].label() + "_" + _shapes[1].label();
}

std::optional<TwoCenterPairComponent>
TwoCenterPairComponent::shift(const char axis,
                              const int  value,
                              const int  center) const
{
    if (const auto tcomp = _shapes[center].shift(axis, value))
    {
        auto new_shapes = shapes();
        
        new_shapes[center] = *tcomp;
        
        return TwoCenterPairComponent(_names, new_shapes);
        
    }
    else
    {
        return std::nullopt;
    }
}

