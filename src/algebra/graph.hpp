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

#ifndef graph_hpp
#define graph_hpp

#include <vector>
#include <set>
#include <algorithm>
#include <tuple>
#include <map>

#include "generics.hpp"
#include "signature.hpp"

/// Recursion group class.
template <class T>
class Graph
{
    /// Vertices of graph.
    std::vector<T> _vertices;
    
    /// Adjacency list (edges) of graph vertices.
    std::vector<std::set<size_t>> _edges;
  
public:
    /// Creates an empty graph.
    Graph();
    
    /// Creates a  graph from the vertices and adjacency list.
    /// @param vertices The vector of vertices.
    /// @param edges The two dimensional adjacency list (edges).
    Graph(const std::vector<T>&                vertices,
          const std::vector<std::set<size_t>>& edges);
    
    /// Creates a  graph from the single vertice.
    /// @param vertice The vertice to create graph.
    Graph(const T& vertice);
    
    /// Retrieves requested vertice from graph.
    /// @param index The index of vertice.
    /// @return The  requested vertice.
    const T& operator[](const size_t index) const;
    
    /// Compares this graph with other graph.
    /// @param other The other graph to compare.
    /// @return true if graphs are equal, false otherwise.
    bool operator==(const Graph<T>& other) const;
    
    /// Compares this graph with other graph.
    /// @param other The other graph to compare.
    /// @return true if graphs are not equal, false otherwise.
    bool operator!=(const Graph<T>& other) const;
    
    /// Compares this graph with other graph.
    /// @param other The other graph to compare.
    /// @return true if this graph is less than other graph, false otherwise.
    bool operator<(const Graph<T>& other) const;
    
    /// Adds vertice to graph.
    /// @param vertice The vertice to be added.
    /// @param root The index of root vertice to which added vertice is connected.
    void add(const T&     vertice,
             const size_t root);
    
    /// Adds vertice to graph.
    /// @param vertice The vertice to be added.
    /// @param root The root vertice to which added vertice is connected.
    void add(const T& vertice,
             const T& root);
    
    /// Replaces selected vertice with given vertice.
    /// @param vertice The replacement vertice.
    /// @param index The index of vertice to be replaced.
    void replace(const T&     vertice,
                 const size_t index);
    
    /// Replaces selected vertice with given vertice.
    /// @param ivertice The index of vertice to merge into.
    /// @param jvertice The index of vertice to be merged.
    void merge(const size_t ivertice,
               const size_t jvertice);
    
    /// Creates inverted graph from this graph.
    /// @return The inverted graph.
    Graph<T> invert() const;
    
    /// Sorts this graph according to normal or reverse order.
    /// @param rorder The flag to indicate the reverse ordering in solrting.
    template <class U>
    void sort(const bool rorder);
    
    /// Reduces this graph to minimal representation by merging similar vertices.
    void reduce();
    
    /// Number of vertices in this graph.
    /// @return The numbet of vertices.
    int vertices() const;
    
    /// Gets edges data for specific vertice.
    /// @param index The index of vertice.
    /// @return The edges data.
    std::set<size_t> edge(const size_t index) const;
    
    /// Gets vector of indexes of orphaned vertices.
    /// @return The vector of indexes.
    std::vector<size_t> orphans() const;
    
    /// Gets vector of indexes for ordering vertices.
    /// @param rorder The reverse order flag.
    /// @return The vector of indexes.
    template <class U>
    std::vector<size_t> indexes(const bool rorder) const;
    
    /// Gets of map of signature:vertice pairs.
    /// @return The map of signature:vertice pairs.
    template <class U>
    std::map<Signature<U>, T> signatures() const;
    
    /// Gets of global signature.
    /// @return The global signature.
    template <class U>
    Signature<U> global_signature() const;
    
    /// Gets unique factors in this graph.
    /// @return The unique recursion factors.
    std::set<Factor> factors() const;
    
    /// Gets map of factors in graph.
    /// @return The map of factors in graph.
    std::map<Factor, int> map_of_factors() const;
};

template <class T>
Graph<T>::Graph()
    
    : _vertices(std::vector<T>({}))
    
    , _edges(std::vector<std::set<size_t>>({}))
{
    
}

template <class T>
Graph<T>::Graph(const std::vector<T>&                vertices,
                const std::vector<std::set<size_t>>& edges)

    : _vertices(vertices)

    , _edges(edges)
{
    
}

template <class T>
Graph<T>::Graph(const T& vertice)

    : _vertices({vertice,})

    , _edges(std::vector<std::set<size_t>>({{}, }))
{
    
}

template <class T>
const T&
Graph<T>::operator[](const size_t index) const
{
    return _vertices[index]; 
}

template <class T>
bool
Graph<T>::operator==(const Graph<T>& other) const
{
    if (this == &other) return true;

    if (_vertices != other._vertices)
    {
        return false;
    }
    else
    {
        return _edges == other._edges;
    }
}

template <class T>
bool
Graph<T>::operator!=(const Graph<T>& other) const
{
    return !((*this) == other);
}

template <class T>
bool
Graph<T>::operator<(const Graph<T>& other) const
{
    if (_vertices != other._vertices)
    {
        return _vertices < other._vertices;
    }
    else
    {
        return _edges < other._edges;
    }
}

template <class T>
void
Graph<T>::add(const T&     vertice,
              const size_t root)
{
    if (const auto tval = std::find(_vertices.cbegin(), _vertices.cend(), vertice); tval != _vertices.cend())
    {
        _edges[root].insert(static_cast<int>(tval - _vertices.cbegin()));
    }
    else
    {
        _edges[root].insert(static_cast<size_t>(_vertices.size()));
        
        _vertices.push_back(vertice);
        
        _edges.push_back(std::set<size_t>());
    }
}

template <class T>
void
Graph<T>::add(const T& vertice,
              const T& root)
{
    if (const auto tval = std::find(_vertices.cbegin(), _vertices.cend(), root); tval != _vertices.cend())
    {
        const auto idx = tval - _vertices.cbegin();
        
        if (const auto rval = std::find(_vertices.cbegin(), _vertices.cend(), vertice); rval != _vertices.cend())
        {
            _edges[idx].insert(static_cast<int>(rval - _vertices.cbegin()));
        }
        else
        {
            _edges[idx].insert(static_cast<int>(_vertices.size()));
            
            _vertices.push_back(vertice);
            
            _edges.push_back(std::set<size_t>());
        }
    }
}

template <class T>
void
Graph<T>::replace(const T&     vertice,
                  const size_t index)
{
    _vertices[index] = vertice; 
}

template <class T>
void
Graph<T>::merge(const size_t ivertice,
                const size_t jvertice)
{
    // update vertices data
    
    gen::merge<T>(_vertices[ivertice], _vertices[jvertice]);
    
    _vertices.erase(_vertices.begin() + jvertice);
    
    // update edges data
    
    _edges[ivertice].insert(_edges[jvertice].begin(), _edges[jvertice].end());
    
    _edges[ivertice].erase(jvertice);
    
    _edges.erase(_edges.begin() + jvertice);
    
    // update indexing scheme in edges data
    
    const auto nedges = _edges.size();
    
    for (size_t i = 0; i < nedges; i++)
    {
        if (_edges[i].erase(jvertice) > 0)
        {
            if (i < ivertice )
            {
                _edges[i].insert(ivertice);
            }
        }

       for (size_t j = jvertice + 1; j < nedges + 1; j++)
       {
            if (_edges[i].erase(j) > 0)
            {
                _edges[i].insert(j - 1);
            }
       }
    }
}

template <class T>
Graph<T>
Graph<T>::invert() const
{
    std::vector<T> new_vertices(_vertices.crbegin(), _vertices.crend());
    
    std::vector<std::set<size_t>> new_edges(_vertices.size(), std::set<size_t>());
    
    const auto nverts = vertices();
    
    for (int i = 0; i < nverts; i++)
    {
        const auto idx = nverts - i - 1;
        
        for (int j = 0; j < nverts; j++)
        {
            if (const auto tval = _edges[j].find(idx); tval != _edges[j].cend())
            {
                new_edges[i].insert(nverts - j - 1);
            }
        }
    }
    
    return Graph(new_vertices, new_edges);
}

template <class T>
template <class U>
void
Graph<T>::sort(const bool rorder)
{
    // intialize new vertices and edges
    
    std::vector<T> new_vertices(_vertices.size(), T());
    
    std::vector<std::set<size_t>> new_edges(_vertices.size(), std::set<size_t>());
    
    // set up ordering indexes
    
    auto vecids = indexes<U>(rorder);
    
    // reconstruct graph
    
    const auto nverts = vertices();
    
    for (size_t i = 0; i < nverts; i++)
    {
        const auto idx = vecids[i];
        
        new_vertices[idx] = _vertices[i];
        
        for (const auto j : _edges[i])
        {
            new_edges[idx].insert(vecids[j]);
        }
    }
    
    // assign new vertices and edges
    
    _vertices = new_vertices;
    
    _edges = new_edges;
}

template <class T>
void
Graph<T>::reduce()
{
    bool need_merge = true;
    
    while (need_merge)
    {
        need_merge = false;
                
        for (int32_t i = 0; i < vertices(); i++)
        {
            auto pair = std::make_tuple(i, i);
            
            for (int32_t j = i + 1; j < vertices(); j++)
            {
                if (gen::similar(_vertices[i], _vertices[j]))
                {
                    pair = std::make_tuple(i, j);
                    
                    break;
                }
            }
            
            if (const auto [idx, jdx] = pair; idx != jdx)
            {
                merge(idx, jdx);
                
                need_merge = true;
                
                break;
            }
        }
    }
}

template <class T>
int
Graph<T>::vertices() const
{
    return static_cast<int>(_vertices.size());
}

template <class T>
std::set<size_t>
Graph<T>::edge(const size_t index) const
{
    return _edges[index]; 
}

template <class T>
std::vector<size_t>
Graph<T>::orphans() const
{
    std::vector<size_t> indexes;
    
    for (size_t i = 0; i < _vertices.size(); i++)
    {
        if (_edges[i].empty())
        {
            indexes.push_back(i);
        }
    }
    
    return indexes;
}

template <class T>
template <class U>
std::vector<size_t>
Graph<T>::indexes(const bool rorder) const
{
    std::set<U> svalues;
    
    for (const auto& tvert : _vertices)
    {
        if (const auto tval = gen::base<T, U>(tvert))
        {
            svalues.insert(*tval);
        }
    }
    
    std::vector<size_t> vecids;
    
    if (_vertices.size() == svalues.size())
    {
        std::vector<U> vvalues;
        
        if (rorder)
        {
            vvalues = std::vector<U>(svalues.rbegin(), svalues.rend());
        }
        else
        {
            vvalues = std::vector<U>(svalues.begin(), svalues.end());
        }
        
        for (const auto& tvert : _vertices)
        {
            if (const auto tval = gen::base<T, U>(tvert))
            {
                auto idx = std::find(vvalues.begin(), vvalues.end(), *tval);
                
                vecids.push_back(static_cast<size_t>(idx - vvalues.begin()));
            }
        }
    }
    
    return vecids;
}

template <class T>
template <class U>
std::map<Signature<U>, T>
Graph<T>::signatures() const
{
    std::map<Signature<U>, T> vmap;
    
    for (const auto& vert : _vertices)
    {
        const Signature<U> tsign = vert.signature();
        
        if (vmap.find(tsign) == vmap.end())
        {
            auto mvert = T(vert);
            
            mvert.reduce();
            
            vmap[tsign] = mvert;
        }
    }
    
    return vmap; 
}

template <class T>
template <class U>
Signature<U>
Graph<T>::global_signature() const
{
    Signature<U> tsign;
    
    for (const auto& vert : _vertices)
    {
        tsign.merge(vert.signature());
    }
    
    return tsign;
}

template <class T>
std::set<Factor>
Graph<T>::factors() const
{
    std::set<Factor> sfacts;
    
    for (const auto& vert : _vertices)
    {
        if (const auto facts = vert.factors(); !facts.empty())
        {
            sfacts.insert(facts.cbegin(), facts.cend());
        }
    }
    
    return sfacts;
}

template <class T>
std::map<Factor, int>
Graph<T>::map_of_factors() const
{
    std::map<Factor, int> mfacts;
    
    for (const auto& vert : _vertices)
    {
        for (const auto& fact : vert.map_of_factors())
        {
            if (const auto tkval = mfacts.find(fact.first); tkval != mfacts.end())
            {
                mfacts[tkval->first] = tkval->second + fact.second;
            }
            else
            {
                mfacts[fact.first] = fact.second;
            }
        }
    }
    
    return mfacts;
}

template <class T>
using VGraphs = std::vector<Graph<T>>;

template <class T>
using VDynGraphs = std::vector<Graph<T>*>;

#endif /* graph_hpp */
