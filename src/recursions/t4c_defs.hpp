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

#ifndef t4c_defs_hpp
#define t4c_defs_hpp

#include <utility>

#include "recursion_group.hpp"
#include "recursion_expansion.hpp"
#include "recursion_term.hpp"
#include "integral_component.hpp"
#include "integral.hpp"
#include "two_center_pair_component.hpp"
#include "two_center_pair.hpp"

using T2CPair = TwoCenterPairComponent;

using T4CIntegral = IntegralComponent<T2CPair, T2CPair>;

using VT4CIntegrals = VIntegralComponents<T2CPair, T2CPair>;

using ST4CIntegrals = SIntegralComponents<T2CPair, T2CPair>;

using R4CTerm = RecursionTerm<T4CIntegral>;

using V4CTerms = VRecursionTerms<T4CIntegral>;

using R4CDist = RecursionExpansion<T4CIntegral>;

using R4Group = RecursionGroup<T4CIntegral>;

using V4Groups = VRecursionGroups<T4CIntegral>;

using I2CPair = TwoCenterPair;

using I4CIntegral = Integral<I2CPair, I2CPair>;

using SI4CIntegrals = SIntegrals<I2CPair, I2CPair>;

using G4Term = std::pair<std::array<int, 4>, I4CIntegral>;

using SG4Terms = std::set<G4Term>;

#endif /* t4c_defs_hpp */
