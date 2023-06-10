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

#ifndef string_formater_hpp
#define string_formater_hpp

#include <cstdint>

#include <string>
#include <vector>

namespace fstr {  // fstr namespace

/**
 Creates uppercased string from string.
 
 @param source the string.
 @return the uppercased string.
 */
std::string
upcase(const std::string& source);

/**
 Creates lowercased string from string.
 
 @param source the string.
 @return the lowercased string.
 */
std::string
lowercase(const std::string& source);

} // fstr namespace

#endif /* string_formater_hpp */
