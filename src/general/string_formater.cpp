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

#include "string_formater.hpp"

#include <cctype>

#include <algorithm>
#include <iterator>

namespace fstr {  // fstr namespace

std::string
upcase(const std::string& source)
{
    std::string str;
    
    std::transform(source.cbegin(),
                   source.cend(),
                   std::back_inserter(str),
                   [] (unsigned char c) { return std::toupper(c); });
    
    return str;
}

std::string
lowercase(const std::string& source)
{
    std::string str;
    
    std::transform(source.cbegin(),
                   source.cend(),
                   std::back_inserter(str),
                   [] (unsigned char c) { return std::tolower(c); });
    
    return str;
}

}  // namespace fstr
