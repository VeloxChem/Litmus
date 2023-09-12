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

#ifndef t4c_utils_hpp
#define t4c_utils_hpp

#include <string>

#include "t4c_defs.hpp"

namespace t4c { // t4c namespace

/// Gets standart capitalized label of integral.
/// @param integral The base four center integral.
/// @return The standart capitalized label of integral.
std::string integral_label(const I4CIntegral& integral);

} // t4c namespace

#endif /* t4c_utils_hpp */
