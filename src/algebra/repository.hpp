// LITMUS: An Automated Molecular params Generator
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

#ifndef repository_hpp
#define repository_hpp

#include <vector>
#include <map>

#include "signature.hpp"
#include "recursion_group.hpp"
#include "graph.hpp"

// Repository class.
template <class T, class U>
class Repository
{
    /// Vector of graphs.
    VGraphs<T> _graphs;
    
    /// Map of unique  signature : data  pairs.
    std::map<Signature<U>, T> _rgmap;

public:
    /// Creates an empty repository.
    Repository();
    
    /// Creates a repository for vector of  graphs and signatures map.
    /// @param graphs The vector of graphs.
    /// @param rgmap The signatures map of unique signature: data pairs .
    Repository(const VGraphs<T>&                graphs,
               const std::map<Signature<U>, T>& rgmap);
    
    /// Compares this repository with other repository.
    /// @param other The other repository to compare.
    /// @return true if repositories, false otherwise.
    bool operator==(const Repository<T, U>& other) const;
    
    /// Compares this repository with other repository.
    /// @param other The other repository to compare.
    /// @return true if repositories are not equal, false otherwise.
    bool operator!=(const Repository<T, U>& other) const;
    
    /// Adds vector of  graphs into this repository.
    /// @param graphs The vector of graphs.
    void add(const VGraphs<T>& graphs);
};

template <class T, class U>
Repository<T, U>::Repository()
    
    : _graphs()

    , _rgmap()
{
    
}

template <class T, class U>
Repository<T, U>::Repository(const VGraphs<T>&                graphs,
                             const std::map<Signature<U>, T>& rgmap)

    : _graphs(graphs)

    , _rgmap(rgmap)
{
    
}

template <class T, class U>
bool
Repository<T, U>::operator==(const Repository<T, U>& other) const
{
    if (this == &other) return true;

    if (_graphs != other._graphs)
    {
        return false;
    }
    else
    {
        return _rgmap == other._rgmap;
    }
}

template <class T, class U>
bool
Repository<T, U>::operator!=(const Repository<T, U>& other) const
{
    return !((*this) == other);
}


template <class T, class U>
void
Repository<T, U>::add(const VGraphs<T>& graphs)
{
    for (const auto& tval : graphs)
    {
        _graphs.push_back(tval);
        
        _rgmap.merge(tval.template signatures<U>());
    }
}

#endif /* repository_hpp */
