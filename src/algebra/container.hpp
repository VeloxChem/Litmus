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

#ifndef container_hpp
#define container_hpp

#include <vector>
#include <optional>
#include <algorithm>

/// Container class.
template <class T>
class Container
{
    /// Vector of recursion groups.
    std::vector<T> _rec_groups;
    
public:
    /// Creates an empty container.
    Container();
    
    /// Creates a  container from the vector of recursion groups.
    /// @param rec_groups The vector of recursion groups.
    Container(const std::vector<T>& rec_groups);
    
    /// Creates a  container from the recursion group.
    /// @param rec_group The recursion group.
    Container(const T& rec_group);
    
    /// Retrieves requested recursion group from group.
    /// @param index The index of vertice.
    /// @return The  requested recursion group.
    const T& operator[](const size_t index) const;
    
    /// Compares this container with other container.
    /// @param other The other container to compare.
    /// @return true if containers are equal, false otherwise.
    bool operator==(const Container<T>& other) const;
    
    /// Compares this container with other container.
    /// @param other The other container to compare.
    /// @return true if containers are not equal, false otherwise.
    bool operator!=(const Container<T>& other) const;
    
    /// Adds recursion group to container.
    /// @param rec_group The recursion group to be added.
    void add(const T& rec_group);
    
    /// Replaces selected recursion group with given recursion group.
    /// @param rec_group The replacement recursion group.
    /// @param index The index of recursion group to be replaced.
    void replace(const T&     rec_group,
                 const size_t index);
    
    /// Reduces this graph to minimal representation by merging similar vertices.
    template <class U>
    void reduce();
    
    /// Number of vertices in this graph.
    /// @return The numbet of vertices.
    size_t recursion_groups() const;
    
    /// Gets base integral of this graph.
    /// @return The base integral of this graph.
    template <class U>
    std::optional<U> base() const;
};

template <class T>
Container<T>::Container()
    
    : _rec_groups(std::vector<T>({}))
{
    
}

template <class T>
Container<T>::Container(const std::vector<T>& rec_groups)

    : _rec_groups(rec_groups)
{
    
}

template <class T>
Container<T>::Container(const T& rec_group)

    : _rec_groups({rec_group,})
{
    
}

template <class T>
const T&
Container<T>::operator[](const size_t index) const
{
    return _rec_groups[index];
}

template <class T>
bool
Container<T>::operator==(const Container<T>& other) const
{
    if (this == &other) return true;

    return _rec_groups == other._rec_groups;
}

template <class T>
bool
Container<T>::operator!=(const Container<T>& other) const
{
    return !((*this) == other);
}

template <class T>
void
Container<T>::add(const T& rec_group)
{
    _rec_groups.push_back(rec_group);
}

template <class T>
void
Container<T>::replace(const T&     rec_group,
                      const size_t index)
{
    _rec_groups[index] = rec_group;
}

template <class T>
template <class U>
void
Container<T>::reduce()
{
    std::vector<T> new_rec_groups;
    
    std::vector<size_t> indexes;
    
    const auto ngroups = _rec_groups.size();
    
    for (size_t i = 0; i < ngroups; i++)
    {
        if (std::find(indexes.begin(), indexes.end(), i) == indexes.end())
        {
            auto tgroup = _rec_groups[i];
            
            auto tbase = *tgroup. template base<U>();
            
            for (size_t j = i + 1; j < ngroups; j++)
            {
                auto rgroup = _rec_groups[j];
                
                auto rbase = *rgroup. template base<U>();
                
                if (tbase == rbase)
                {
                    tgroup.merge(rbase);
                    
                    indexes.push_back(j);
                }
            }
            
            new_rec_groups.push_back(tgroup);
        }
    }
    
    _rec_groups = new_rec_groups; 
}

template <class T>
size_t
Container<T>::recursion_groups() const
{
    return _rec_groups.size();
}

template <class T>
template <class U>
std::optional<U>
Container<T>::base() const
{
    if (!_rec_groups.empty())
    {
        return _rec_groups[0]. template base<U>();
    }
    else
    {
        return std::nullopt;
    }
}

#endif /* container_hpp */
