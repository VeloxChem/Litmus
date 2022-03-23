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

/// Recursion group class.
template <class T>
class Graph
{
    /// Vertices of graph.
    std::vector<T> _vertices;
    
    /// Adjacency list (edges) of graph vertices.
    std::vector<std::set<int>> _edges;
  
public:
    /// Creates an empty graph.
    Graph();
    
    /// Creates a  graph from the vertices and adjacency list.
    /// @param vertices The vector of vertices.
    /// @param edges The two dimensional adjacency list (edges).
    Graph(const std::vector<T>&             vertices,
          const std::vector<std::set<int>>& edges);
    
    /// Creates a  graph from the single vertice.
    /// @param vertice The vertice to create graph.
    Graph(const T& vertice);
    
    /// Retrieves requested vertice from graph.
    /// @param index The index of vertice.
    /// @return The  requested vertice.
    const T& operator[](const int index) const;
    
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
    void add(const T&  vertice,
             const int root);
    
    /// Adds vertice to graph.
    /// @param vertice The vertice to be added.
    /// @param root The root vertice to which added vertice is connected.
    void add(const T& vertice,
             const T& root);
    
    /// Creates inverted graph from this graph.
    /// @return The inverted graph.
    Graph<T> invert() const;
    
    /// Number of vertices in this graph.
    /// @return The numbet of vertices.
    int vertices() const;
};

template <class T>
Graph<T>::Graph()
    
    : _vertices(std::vector<T>({}))
    
    , _edges(std::vector<std::set<int>>({}))
{
    
}

template <class T>
Graph<T>::Graph(const std::vector<T>&             vertices,
                const std::vector<std::set<int>>& edges)

    : _vertices(vertices)

    , _edges(edges)
{
    
}

template <class T>
Graph<T>::Graph(const T& vertice)

    : _vertices({vertice,})

, _edges(std::vector<std::set<int>>({{}, }))
{
    
}

template <class T>
const T&
Graph<T>::operator[](const int index) const
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
Graph<T>::add(const T&  vertice,
              const int root)
{
    if (const auto tval = std::find(_vertices.cbegin(), _vertices.cend(), vertice); tval != _vertices.cend())
    {
        _edges[root].insert(static_cast<int>(tval - _vertices.cbegin()));
    }
    else
    {
        _edges[root].insert(static_cast<int>(_vertices.size()));
        
        _vertices.push_back(vertice);
        
        _edges.push_back(std::set<int>());
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
            
            _edges.push_back(std::set<int>());
        }
    }
}

template <class T>
Graph<T>
Graph<T>::invert() const
{
    std::vector<T> new_vertices(_vertices.crbegin(), _vertices.crend());
    
    std::vector<std::set<int>> new_edges(_vertices.size(), std::set<int>());
    
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
int
Graph<T>::vertices() const
{
    return static_cast<int>(_vertices.size());
}

#endif /* graph_hpp */
