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

#include "operator.hpp"


Operator::Operator()

    : _name("")

    , _shape(Tensor(0))

    , _target("")

    , _center(-1)
{
    
}

Operator::Operator(const std::string& name,
                   const Tensor&      shape,
                   const std::string& target,
                   const int          center)
    
    : _name(name)

    , _shape(shape)

    , _target(target)

    , _center(center)
{
    
}

//Operator::Operator(const OperatorComponent& opcomp)
//
//    : _name(opcomp.name())
//
//    , _shape(Tensor(opcomp.shape()))
//
//    , _target(opcomp.target())
//{
//    for (const auto& tcomp : opcomp.shapes())
//    {
//        _shapes.push_back(tcomp.order());
//    }
//}

bool
Operator::operator==(const Operator& other) const
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
Operator::operator!=(const Operator& other) const
{
    return !((*this) == other);
}

bool
Operator::operator<(const Operator& other) const
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

std::string
Operator::to_string() const
{
    std::string str = "{" + _name + ":";
    
    //for (const auto& tcomp : _shapes) str += tcomp.to_string();

    return str + "}";
}

std::string
Operator::label() const
{
    std::string str;
    
    //for (const auto& tcomp : _shapes) str += tcomp.label();
    
    return str;
}

VOperatorComponents
Operator::components() const
{
    VOperatorComponents opcomps;
    
//    for (const auto& tcomps : make_components<TensorComponent>(_shapes))
//    {
//        opcomps.push_back(OperatorComponent(_name, tcomps));
//    }
    
    return opcomps;
}
