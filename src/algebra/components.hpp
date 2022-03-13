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

#ifndef components_hpp
#define components_hpp

#include <vector>

/// Template function for converting vector of tensor like objects
/// into direct product of tensor component like objects.
/// @param tvalues the vector of tensor like objects.
/// @return the vector of tensor component like objects
template <class T, class U>
std::vector<std::vector<T>>
make_components(const std::vector<U>& tvalues)
{
    std::vector<std::vector<T>> mtcomps;
    
    if (tvalues.size() > 0)
    {
        for (const auto& tcomp : tvalues[0].components())
        {
            mtcomps.push_back(std::vector({tcomp}));
        }
        
        for (size_t i = 1; i < tvalues.size(); i++)
        {
            auto mccomps = mtcomps;
            
            mtcomps.clear();
            
            for (const auto& vcomps : mccomps)
            {
                for (const auto& tcomp : tvalues[i].components())
                {
                    mtcomps.push_back(vcomps);
                    
                    (mtcomps.rbegin())->push_back(tcomp);
                }
            }
        }
    }
    
    return mtcomps;
}

#endif /* components_hpp */
