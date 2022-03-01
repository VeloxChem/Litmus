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

bool
OperatorComponent::operator==(const OperatorComponent& other) const
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
    else
    {
        return _shape < other._shape;
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
    std::string str = "{" + _name + ":";
    
    //for (const auto& tcomp : _shapes) str += tcomp.to_string();

    return str + "}";
}

std::string
OperatorComponent::label() const
{
    std::string str;
    
    //for (const auto& tcomp : _shapes) str += tcomp.label();
    
    return str;
}

std::optional<OperatorComponent>
OperatorComponent::shift(const int  center,
                         const char axis,
                         const int  value,
                         const bool noscalar) const
{
    if (const auto tcomp = _shape.shift(axis, value))
    {
//        auto new_shapes = _shapes;
//
//        if (noscalar && (tcomp->order() == 0))
//        {
//            new_shapes.erase(new_shapes.begin() + center);
//
//            if (new_shapes.empty()) return std::nullopt;
//        }
//        else
//        {
//            new_shapes[center] = *tcomp;
//        }
//
//        return OperatorComponent(_name, new_shapes);
    }
    else
    {
        return std::nullopt;
    }
    
    return std::nullopt;
}
