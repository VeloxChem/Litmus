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

#ifndef unique_map_hpp
#define unique_map_hpp

#include <set>
#include <map>

/// Map of unique components class.
template <class T, class U>
class UniqueMap
{
    /// Map of unique components of tensorial values.
    std::map<T, std::set<U>> _components;
    
public:
    /// Creates a map of unique components.
    UniqueMap();
    
    /// Adds unique component.
    /// @param component The component to be added.
    void add(const U& component);
    
    /// Adds set of components.
    /// @param component The set of components to be added.
    void add(const std::set<U>& components);
    
    /// Checks if unique components map contains this component.
    /// @param component The component to be found.
    bool find(const U& component) const
    
    /// Gets number of unique components in unique components map.
    /// @return The number of unique components.
    size_t size() const;
};

template <class T, class U>
UniqueMap<T, U>::UniqueMap()
    
    : _components(std::map<T, std::set<U>>())
{
    
}

template <class T, class U>
void
UniqueMap<T, U>::add(const U& component)
{
    const auto key = T(component);
    
    if (auto tval = _components.find(key); tval != _components.end())
    {
        _components[key].insert(component);
    }
    else
    {
        _components[key] = std::set<U>({component, });
    }
}

template <class T, class U>
void
UniqueMap<T, U>::add(const std::set<U>& components)
{
    for (const auto& tcomp : components) add(tcomp);
}

template <class T, class U>
bool
UniqueMap<T, U>::find(const U& component) const
{
    const auto key = T(component);
    
    if (auto tval = _components.find(key); tval != _components.end())
    {
        return _components[key].find(component) != _components[key].end();
    }
    else
    {
        return false;
    }
}

template <class T, class U>
size_t
UniqueMap<T, U>::size() const
{
    size_t ncomps = 0;
    
    for (const auto& [key, tval] : _components)
    {
        ncomps += tval.size();
    }
    
    return ncomps;
}

#endif /* unique_map_hpp */
