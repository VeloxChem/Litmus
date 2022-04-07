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
#include <set>
#include <iostream>

#include "signature.hpp"
#include "recursion_group.hpp"
#include "graph.hpp"

// Repository class.
template <class T, class U>
class Repository
{
    /// Vector of graphs.
    VDynGraphs<T> _graphs;
    
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
    
    /// Destroys a repository.
    ~Repository();
    
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
    
    /// Gets  set  of base integrals in this repository.
    /// @return The sef of base integrals in this repository.
    template <class V>
    std::set<V> base() const;
    
    /// Gets  map of unique signature : data pairs with same base integral.
    /// @param bdata The base integral of signature.
    /// @return The sef of base integrals in this repository.
    template <class V>
    std::map<Signature<U>, T> base_map(const V& bdata) const;
    
    /// Computes number of recursion groups in repository.
    /// @return The number of recursion groups.
    size_t rec_groups() const;
    
    /// Prints summary of this repository.
    void summary() const;
    
    /// Prints detailes summary of signatures of this repository.
    template <class V>
    void details() const;
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

    : _graphs()

    , _rgmap(rgmap)
{
        for (const auto&  tval : graphs)
        {
            _graphs.push_back(new Graph<T>(tval));
        }
}

template <class T, class U>
Repository<T, U>::~Repository()
{
    for (auto& tval : _graphs)
    {
        delete tval;
    }
}

template <class T, class U>
bool
Repository<T, U>::operator==(const Repository<T, U>& other) const
{
    if (this == &other) return true;

    if (_rgmap != other._rgmap)
    {
        return false;
    }
    else
    {
        if (_graphs.size() == other._graphs.size())
        {
            for (size_t i = 0; i < _graphs.size(); i++)
            {
                if ((*_graphs[i]) != (*other._graphs[i]))
                {
                    return false;
                }
            }
            
            return true; 
        }
        else
        {
            return false;
        }
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
        _graphs.push_back(new Graph<T>(tval));
        
        _rgmap.merge(tval.template signatures<U>());
    }
}

template <class T, class U>
template <class V>
std::set<V>
Repository<T, U>::base() const
{
    std::set<V> sints;
    
    for (const auto& tval : _rgmap)
    {
        if (const auto tint = tval.first.template base<V>())
        {
            sints.insert(*tint);
        }
    }
    
    return sints;
}

template <class T, class U>
template <class V>
std::map<Signature<U>, T>
Repository<T, U>::base_map(const V& bdata) const
{
    std::map<Signature<U>, T> tmap;
    
    for (const auto& tval : _rgmap)
    {
        if (const auto tint = tval.first.template base<V>())
        {
            if (bdata == *tint)
            {
                tmap[tval.first] = tval.second;
            }
        }
    }
    
    return tmap;
}

template <class T, class U>
size_t
Repository<T, U>::rec_groups() const
{
    size_t ngroups = 0;
    
    for (const auto& tval : _graphs)
    {
        ngroups += tval->vertices();
    }
    
    return ngroups;
}

template <class T, class U>
void
Repository<T, U>::summary() const
{
    std::cout << "Number of Recursion Graphs  : " << _graphs.size() << std::endl;
    
    std::cout << "Number of Recursion Groups  : " << rec_groups() << std::endl;
    
    std::cout << "Number of Unique Signatures : " << _rgmap.size() << std::endl;
}

template <class T, class U>
template <class V>
void
Repository<T, U>::details() const
{
    for (const auto& tint : base<V>())
    {
        const auto tmap = base_map<V>(tint);
        
        std::cout << "Integral angular form [" << tint.label() << "]" << std::endl;
        
        std::cout << "Unique recursion patterns: " << tmap.size() << std::endl;
        
        for (const auto& tval : tmap)
        {
            std::cout << " o-params: " << tval.first.nparams("out");
            
            std::cout << " i-params: " << tval.first.nparams("inp");;
         
            std::cout << " i-facts: " << tval.first.nfactors() << std::endl;
        }
    }
}

#endif /* repository_hpp */
