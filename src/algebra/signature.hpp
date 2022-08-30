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

#ifndef signature_hpp
#define signature_hpp

#include <set>
#include <string>
#include <optional>

#include "factor.hpp"

/// Signature class.
template <class T>
class Signature
{
    /// Set of global params.
    std::set<T> _glob_params;
    
    /// Set of output params.
    std::set<T> _out_params;
    
    /// Set of input params.
    std::set<T> _inp_params;
    
    /// Set of input factors.
    std::set<Factor> _factors;
  
public:
    /// Creates an empty signature.
    Signature();
    
    /// Creates a  signature from the set of input/output params and recursion factors.
    /// @param glob_params The set of global params.
    /// @param out_params The set of output params.
    /// @param inp_params The set of input params.
    /// @param factors The set of recursion factors.
    Signature(const std::set<T>&      glob_params,
              const std::set<T>&      out_params,
              const std::set<T>&      inp_params,
              const std::set<Factor>& factors);
    
    /// Compares this signature with other signature.
    /// @param other The other signature to compare.
    /// @return true if signatures are equal, false otherwise.
    bool operator==(const Signature<T>& other) const;
    
    /// Compares this signature with other signature.
    /// @param other The other signature to compare.
    /// @return true if signatures are not equal, false otherwise.
    bool operator!=(const Signature<T>& other) const;
    
    /// Compares this signature with other signature.
    /// @param other The other signature to compare.
    /// @return true if this signature is less than other signature, false otherwise.
    bool operator<(const Signature<T>& other) const;
    
    /// Merges other signature into this signature.
    /// @param other The other signature to merge.
    void merge(const Signature<T>& other);
    
    /// Adds  parameter to this signature.
    /// @param param The parameter to add.
    /// @param destination The destination of added parameter.
    void add(const T&           param,
             const std::string& destination);
    
    /// Adds  factor to this signature.
    /// @param factor The factor to add.
    void add(const Factor& factor);
    
    /// Gets  base integral of uniform signature.
    /// @return The base integral of uniform signature.
    template <class U>
    std::optional<U> base() const;
    
    /// Gets  set of unique recursion expansion integrals used in signature.
    /// @return The set of unique recursion expansion integrals.
    template <class U>
    std::set<U> expansion() const;
    
    /// Gets  set of unique integral components for given integral insignature.
    /// @return The set of unique integral components.
    template <class U>
    std::set<T> expansion_components(const U& integral) const; 
    
    /// Gets number of factors in this signature.
    /// @return The number of factors in  signature.
    size_t nfactors() const;
    
    /// Gets number of  parameters in this signature.
    /// @param destination The destination of counted parameters.
    /// @return The number of parameters in this signature.
    size_t nparams(const std::string& destination) const;
    
    /// Gets set of factors in this signature.
    /// @return The set of factors in  signature.
    std::set<Factor> factors() const;
    
    /// Gets set of factors with the given name in this signature.
    /// @param name The name of factor.
    /// @return The set of factors with specific name in  signature.
    std::set<Factor> factors(const std::string& name) const;
    
    /// Gets set of factor names in this signature.
    /// @return The set of factors in  signature.
    std::set<std::string> factor_names()const;
    
    /// Gets parameters in this signature.
    /// @param destination The destination of parameters.
    /// @return The selected parameters in this signature.
    std::set<T> params(const std::string& destination) const;
    
    /// Gets parameters in this signature.
    /// @param integral The integral to match parameter components. 
    /// @param destination The destination of parameters.
    /// @return The selected parameters in this signature.
    template <class U>
    std::set<T> params(const U&           integral,
                       const std::string& destination) const;
};

template <class T>
Signature<T>::Signature()
    
    : _glob_params(std::set<T>({}))

    , _out_params(std::set<T>({}))

    , _inp_params(std::set<T>({}))

    , _factors(std::set<Factor>({}))
{
    
}

template <class T>
Signature<T>::Signature(const std::set<T>&      glob_params,
                        const std::set<T>&      out_params,
                        const std::set<T>&      inp_params,
                        const std::set<Factor>& factors)
    
    : _glob_params(glob_params)

    , _out_params(out_params)

    , _inp_params(inp_params)

    , _factors(factors)
{
    
}

template <class T>
bool
Signature<T>::operator==(const Signature<T>& other) const
{
    if (this == &other) return true;

    if (_glob_params != other._glob_params)
    {
        return false;
    }
    else if (_out_params != other._out_params)
    {
        return false;
    }
    else if (_inp_params != other._inp_params)
    {
        return false;
    }
    else
    {
        return _factors == other._factors;
    }
}

template <class T>
bool
Signature<T>::operator!=(const Signature<T>& other) const
{
    return !((*this) == other);
}

template <class T>
bool
Signature<T>::operator<(const Signature<T>& other) const
{
    if (_glob_params != other._glob_params)
    {
        return _glob_params < other._glob_params;
    }
    else if (_out_params != other._out_params)
    {
        return _out_params < other._out_params;
    }
    else if (_inp_params != other._inp_params)
    {
        return _inp_params < other._inp_params;
    }
    else
    {
        return _factors < other._factors;
    }
}

template <class T>
void
Signature<T>::merge(const Signature<T>& other)
{
    _out_params.insert(other._out_params.cbegin(), other._out_params.cend());
    
    _inp_params.insert(other._inp_params.cbegin(), other._inp_params.cend());
    
    _factors.insert(other._factors.cbegin(), other._factors.cend());
}

template <class T>
void
Signature<T>::add(const T&           param,
                  const std::string& destination)
{
    if (destination == "glob")
    {
        _glob_params.insert(param);
    }
    
    if (destination == "inp")
    {
        _inp_params.insert(param);
    }
    
    if (destination == "out")
    {
        _out_params.insert(param);
    }
}

template <class T>
void
Signature<T>::add(const Factor& factor)
{
    _factors.insert(factor); 
}

template <class T>
template <class U>
std::optional<U>
Signature<T>::base() const
{
    std::set<U> tints;
    
    for (const auto& tval : _out_params)
    {
        tints.insert(U(tval));
    }
    
    if (tints.size() == 1)
    {
        return *(tints.begin());
    }
    else
    {
        return std::nullopt;
    }
}

template <class T>
template <class U>
std::set<U>
Signature<T>::expansion() const
{
    std::set<U> tints;
    
    for (const auto& tval : _inp_params)
    {
        tints.insert(U(tval));
    }
    
    return tints;
}

template <class T>
template <class U>
std::set<T>
Signature<T>::expansion_components(const U& integral) const
{
    std::set<T> tints;
    
    for (const auto& tval : _inp_params)
    {
        if (integral == U(tval))
        {
            tints.insert(tval);
        }
    }
    
    return tints;
}

template <class T>
size_t
Signature<T>::nfactors() const
{
    return _factors.size();
}

template <class T>
size_t
Signature<T>::nparams(const std::string& destination) const
{
    if (destination == "glob")
    {
        return _glob_params.size();
    }
    
    if (destination == "inp")
    {
        return _inp_params.size();
    }
    
    if (destination == "out")
    {
        return _out_params.size();
    }
    
    return 0;
}

template <class T>
std::set<Factor>
Signature<T>::factors() const
{
    return _factors;
}

template <class T>
std::set<Factor>
Signature<T>::factors(const std::string& name) const
{
    std::set<Factor> facts;
    
    for (const auto& tval : _factors)
    {
        if (tval.name() == name)
        {
            facts.insert(tval); 
        }
    }
    
    return facts; 
}

template <class T>
std::set<std::string>
Signature<T>::factor_names() const
{
    std::set<std::string> labels;
    
    for (const auto fact : _factors)
    {
        labels.insert(fact.name());
    }
    
    return labels;
}

template <class T>
std::set<T>
Signature<T>::params(const std::string& destination) const
{
    if (destination == "glob")
    {
        return _glob_params;
    }
    
    if (destination == "inp")
    {
        return _inp_params;
    }
    
    if (destination == "out")
    {
        return _out_params;
    }
    
    return std::set<T>();
}

template <class T>
template <class U>
std::set<T>
Signature<T>::params(const U&           integral,
                     const std::string& destination) const
{
    std::set<T> tcomps;
    
    for (const auto& tcomp : params(destination))
    {
        if (U(tcomp) == integral)
        {
            tcomps.insert(tcomp);
        }
    }
    
    return tcomps;
}

#endif /* signature_hpp */
