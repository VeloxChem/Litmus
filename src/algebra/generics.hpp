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

#ifndef generics_hpp
#define generics_hpp

#include <string>
#include <optional>

#include "signature.hpp"

namespace gen
{
/// Template function for merging two objects into one.
/// @param lhsobj the object to merge into.
/// @param rhsobj the object to be merged.
template <class T>
inline void
merge(      T& lhsobj,
      const T& rhsobj)
{
    lhsobj.merge(rhsobj);
}

/// Template function for merging two objects into one.
/// @param lhsobj the object to merge into.
/// @param rhsobj the object to be merged.
template<>
inline void
merge(      std::string& lhsobj,
      const std::string& rhsobj)
{
    lhsobj.append(rhsobj);
}

/// Template function for checking if two objects are similar.
/// @param lhsobj the first object to compare.
/// @param rhsobj the second object to compare.
/// @return True if two objects are similar, false otherwise.
template <class T>
inline bool
similar(const T& lhsobj,
        const T& rhsobj)
{
    return lhsobj.similar(rhsobj);
}

/// Template function for checking if two objects are similar.
/// @param lhsobj the first object to compare.
/// @param rhsobj the second object to compare.
/// @return True if two objects are similar, false otherwise.
template<>
inline bool
similar(const std::string& lhsobj,
        const std::string& rhsobj)
{
    return lhsobj == rhsobj;
}


/// Template function for getting base value of object.
/// @param obj the object.
/// @return the base value of this object.
template <class T, class U>
inline std::optional<U>
base(const T& obj)
{
    return obj.template base<U>();
}

/// Template function for getting base value of object.
/// @param obj the object.
/// @return the base value of this object.
template <> 
inline std::optional<std::string>
base(const std::string& obj)
{
    if (obj.empty())
    {
        return std::nullopt;
    }
    else
    {
        return obj;
    }
}

}

#endif /* generics_hpp */
