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

#include "tensor.hpp"

Tensor::Tensor()

    : _order(0)
{
    
}

Tensor::Tensor(const int order)

    : _order(order)
{
    
}

Tensor::Tensor(const TensorComponent& tcomp)

    : _order(tcomp.order())
{
    
}

bool
Tensor::operator==(const Tensor& other) const
{
    if (this == &other) return true;

    return _order == other._order;
}

bool
Tensor::operator!=(const Tensor& other) const
{
    return !((*this) == other);
}

bool
Tensor::operator<(const Tensor& other) const
{
    return _order < other._order;
}

int
Tensor::order() const
{
    return _order;
}

std::string
Tensor::to_string() const
{
    return "(" + std::to_string(_order) + ")";
}

std::string
Tensor::label() const
{
    if (_order > 16)
        return "l" + std::to_string(_order);
    else
    {
        const std::string names("SPDFGHIKLMNOQRTUV");
        
        return std::string(1, names[_order]);
    }
}

VTensorComponents
Tensor::components() const
{
    VTensorComponents vtcomps({TensorComponent()});
    
    for (int i = 1; i <= _order; i++)
    {
        const auto ctcomps = vtcomps;
            
        vtcomps.clear();
            
        for (const auto axis : std::string("xyz"))
        {
            for (const auto& ctcomp : ctcomps)
            {
                if (const auto tcomp = ctcomp.shift(axis, 1); tcomp->primary() == axis)
                {
                    vtcomps.push_back(tcomp.value());
                }
            }
        }
    }
    
    return vtcomps;
}
