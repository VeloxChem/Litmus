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

#include "factor.hpp"

/// Signature class.
template <class T>
class Signature
{
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
    /// @param out_params The set of output params.
    /// @param inp_params The set of input params.
    /// @param factors The set of recursion factors.
    Signature(const std::set<T>&      out_params,
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
    
    /// Adds  parameter to this signature.
    /// @param param The parameter to add.
    /// @param destination The destination of added parameter.
    void add(const T&           param,
             const std::string& destination);
    
    /// Adds  factor to this signature.
    /// @param factor The factor to add.
    void add(const Factor& factor);
};

template <class T>
Signature<T>::Signature()
    
    : _out_params(std::set<T>({}))

    , _inp_params(std::set<T>({}))

    , _factors(std::set<Factor>({}))
{
    
}

template <class T>
Signature<T>::Signature(const std::set<T>&      out_params,
                        const std::set<T>&      inp_params,
                        const std::set<Factor>& factors)
    
    : _out_params(out_params)

    , _inp_params(inp_params)

    , _factors(factors)
{
    
}

template <class T>
bool
Signature<T>::operator==(const Signature<T>& other) const
{
    if (this == &other) return true;

    if (_out_params != other._out_params)
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
    if (_out_params != other._out_params)
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
Signature<T>::add(const T&           param,
                  const std::string& destination)
{
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

#endif /* signature_hpp */
