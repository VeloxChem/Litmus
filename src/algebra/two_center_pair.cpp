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

#include "two_center_pair.hpp"

TwoCenterPair::TwoCenterPair()

    : _names(std::array<std::string, 2>({"", ""}))

    , _shapes(std::array<Tensor, 2>({Tensor(0), Tensor(0)}))
{
    
}

TwoCenterPair::TwoCenterPair(const std::array<std::string, 2>& names,
                             const std::array<Tensor, 2>&      shapes)
    
    : _names(names)

    , _shapes(shapes)
{
    
}

TwoCenterPair::TwoCenterPair(const std::string& f_name,
                             const int          f_angmom,
                             const std::string& s_name,
                             const int          s_angmom)

    : _names({f_name, s_name})

    , _shapes({Tensor(f_angmom), Tensor(s_angmom)})
{
    
}

TwoCenterPair::TwoCenterPair(const TwoCenterPairComponent& t2pcomp)

    : _names(t2pcomp.names())

    , _shapes(std::array<Tensor, 2>({Tensor(0), Tensor(0)}))
{
    const auto tcomps = t2pcomp.shapes();
        
    _shapes[0] = Tensor(tcomps[0].order());
    
    _shapes[1] = Tensor(tcomps[1].order());
}

bool
TwoCenterPair::operator==(const TwoCenterPair& other) const
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
TwoCenterPair::operator!=(const TwoCenterPair& other) const
{
    return !((*this) == other);
}

bool
TwoCenterPair::operator<(const TwoCenterPair& other) const
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

std::string
TwoCenterPair::to_string() const
{
    return "{" + _names[0] + ":" + _shapes[0].to_string() + ";" +
    
           _names[1] + ":"  + _shapes[1].to_string() + "}";
}

std::string
TwoCenterPair::label() const
{
    return _shapes[0].label() + _shapes[1].label();
}

VTwoCenterPairComponents
TwoCenterPair::components() const
{
    VTwoCenterPairComponents t2pcomps;
    
    for (const auto& fcomp : _shapes[0].components())
    {
        for (const auto& scomp : _shapes[1].components())
        {
            t2pcomps.push_back(TwoCenterPairComponent(_names, {fcomp, scomp})); 
        }
    }
    
    return t2pcomps;
}
