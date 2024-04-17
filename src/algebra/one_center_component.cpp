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

#include "one_center_component.hpp"

OneCenterComponent::OneCenterComponent()

    : _name(std::string())

    , _shape(TensorComponent(0, 0, 0))
{
    
}

OneCenterComponent::OneCenterComponent(const std::string&     name,
                                       const TensorComponent& shape)
    
    : _name(name)

    , _shape(shape)
{
    
}

const TensorComponent&
OneCenterComponent::operator[](const int index) const
{
    return _shape;
}

bool
OneCenterComponent::operator==(const OneCenterComponent& other) const
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
OneCenterComponent::operator!=(const OneCenterComponent& other) const
{
    return !((*this) == other);
}

bool
OneCenterComponent::operator<(const OneCenterComponent& other) const
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

bool
OneCenterComponent::similar(const OneCenterComponent& other) const
{
    if (_name != other._name)
    {
        return false;
    }
    else
    {
        return _shape.similar(other._shape);
    }
}

std::string
OneCenterComponent::to_string() const
{
    return "{" + _name + ":" + _shape.to_string() + "}";
}

std::string
OneCenterComponent::label() const
{
    return _shape.label();
}

std::optional<OneCenterComponent>
OneCenterComponent::shift(const char axis,
                          const int  value,
                          const int  center) const
{
    if (const auto tcomp = _shape.shift(axis, value))
    {
        return OneCenterComponent(_name, *tcomp);
    }
    else
    {
        return std::nullopt;
    }
}
