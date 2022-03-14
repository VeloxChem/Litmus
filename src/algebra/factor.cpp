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

#include "factor.hpp"

Factor::Factor()

    : _name("")

    , _label("")

    , _shape(TensorComponent(0, 0, 0))
{
}

Factor::Factor(const std::string&     name,
               const std::string&     label,
               const TensorComponent& shape)

    : _name(name)
    
    , _label(label)

    , _shape(shape)
{
    
}

bool
Factor::operator==(const Factor& other) const
{
    if (this == &other) return true;

    if (_name != other._name)
    {
        return false;
    }
    else if (_label != other._label)
    {
        return false;
    }
    else
    {
        return _shape == other._shape;
    }
}

bool
Factor::operator!=(const Factor& other) const
{
    return !((*this) == other);
}

bool
Factor::operator<(const Factor& other) const
{
    if (_name != other._name)
    {
        return _name < other._name;
    }
    else if (_label != other._label)
    {
        return _shape < other._shape;
    }
    else
    {
        return _shape < other._shape;
    }
}

std::string
Factor::to_string() const
{
    return "{" + _name + "(" + _label + "):" + _shape.to_string() + "}";
}

std::string
Factor::label() const
{
    if (_shape.order() > 0)
    {
        return _label + "_" + _shape.label();
    }
    else
    {
        return _label;
    }
}

