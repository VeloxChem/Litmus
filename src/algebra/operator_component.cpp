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

#include "operator_component.hpp"

OperatorComponent::OperatorComponent()

    : _name("")

    , _shape(TensorComponent(0, 0, 0))

    , _target("none")

    , _center(-1)
{
}

OperatorComponent::OperatorComponent(const std::string&     name,
                                     const TensorComponent& shape,
                                     const std::string&     target,
                                     const int              center)

    : _name(name)
    
    , _shape(shape)

    , _target(target)

    , _center(center)
{
    
}

int
OperatorComponent::operator[](const char axis) const
{
    return _shape[axis];
}

bool
OperatorComponent::operator==(const OperatorComponent& other) const
{
    if (this == &other) return true;

    if (_name != other._name)
    {
        return false;
    }
    else if (_shape != other._shape)
    {
        return false;
    }
    else if (_target != other._target)
    {
        return false;
    }
    else
    {
        return _center == other._center;
    }
}

bool
OperatorComponent::operator!=(const OperatorComponent& other) const
{
    return !((*this) == other);
}

bool
OperatorComponent::operator<(const OperatorComponent& other) const
{
    if (_name != other._name)
    {
        return _name < other._name;
    }
    else if (_shape != other._shape)
    {
        return _shape < other._shape;
    }
    else if (_target != other._target)
    {
        return _target < other._target;
    }
    else
    {
        return _center < other._center;
    }
}

bool
OperatorComponent::similar(const OperatorComponent& other) const
{
    if (this == &other) return true;

    if (_name != other._name)
    {
        return false;
    }
    else if (_target != other._target)
    {
        return false;
    }
    else if (_center != other._center)
    {
        return false;
    }
    else
    {
        return _shape.similar(other._shape);
    }
}

std::string
OperatorComponent::name() const
{
    return _name;
}

TensorComponent
OperatorComponent::shape() const
{
    return _shape;
}

std::string
OperatorComponent::target() const
{
    return _target;
}

int
OperatorComponent::center() const
{
    return _center;
}

std::string
OperatorComponent::to_string() const
{
    return "{" + _name + ":" + _shape.to_string() + "}" +
    
           "[" + _target + ":" + std::to_string(_center) + "]";
}

std::string
OperatorComponent::label() const
{
    return _shape.label();
}

std::optional<OperatorComponent>
OperatorComponent::shift(const char axis,
                         const int  value,
                         const bool noscalar) const
{
    if (const auto tcomp = _shape.shift(axis, value))
    {
        if (noscalar && (tcomp->order() == 0))
        {
            return std::nullopt;
        }
        else
        {
            return OperatorComponent(_name, *tcomp, _target, _center);
        }
    }
    else
    {
        return std::nullopt;
    }
}
