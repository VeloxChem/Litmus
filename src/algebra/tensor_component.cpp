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

#include "tensor_component.hpp"

TensorComponent::TensorComponent()

    : _ax(0), _ay(0), _az(0)
{
    
}

TensorComponent::TensorComponent(const int ax,
                                 const int ay,
                                 const int az)
    : _ax(ax), _ay(ay), _az(az)
{
    
}

int
TensorComponent::operator[](const char axis) const
{
    if (axis == 'x') return _ax;
    
    if (axis == 'y') return _ay;
    
    if (axis == 'z') return _az;
    
    return -1;
}

bool
TensorComponent::operator==(const TensorComponent& other) const
{
    if (this == &other) return true;

    if (_ax != other._ax)
    {
        return false;
    }
    else if (_ay != other._ay)
    {
        return false;
    }
    else
    {
        return _az == other._az;
    }
}

bool
TensorComponent::operator!=(const TensorComponent& other) const
{
    return !((*this) == other);
}

bool
TensorComponent::operator<(const TensorComponent& other) const
{
    if (_ax != other._ax)
    {
        return _ax < other._ax;
    }
    else if (_ay != other._ay)
    {
        return _ay < other._ay;
    }
    else
    {
        return _az < other._az;
    }
}

std::string
TensorComponent::to_string() const
{
    return "(" + std::to_string(_ax) +
    
           "," + std::to_string(_ay) +
    
           "," + std::to_string(_az) + ")";
}

std::string
TensorComponent::label() const
{
    if (order() == 0)
    {
        return std::string("0");
    }
    else
    {
        return std::string(_ax, 'x') +
        
               std::string(_ay, 'y') +
    
               std::string(_az, 'z');
    }
}

int
TensorComponent::order() const
{
    return _ax + _ay + _az;
}

int
TensorComponent::maximum() const
{
    return std::max(_ax, std::max(_ay, _az));
}

char
TensorComponent::primary() const
{
    if (_ax > 0) return 'x';
    
    if (_ay > 0) return 'y';
    
    if (_az > 0) return 'z';
       
    return 'x';
}

std::optional<TensorComponent>
TensorComponent::shift(const char axis,
                       const int  value) const
{
    if (const auto tax = _ax + value; (axis == 'x') && (tax >= 0))
    {
        return TensorComponent(tax, _ay, _az);
    }
    
    if (const auto tay = _ay + value; (axis == 'y') && (tay >= 0))
    {
        return TensorComponent(_ax, tay, _az);
    }
    
    if (const auto taz = _az + value; (axis == 'z') && (taz >= 0))
    {
        return TensorComponent(_ax, _ay, taz);
    }
        
    return std::nullopt;
}
