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

#ifndef t2c_defs_hpp
#define t2c_defs_hpp

#include "container.hpp"
#include "recursion_group.hpp"
#include "recursion_expansion.hpp"
#include "recursion_term.hpp"
#include "integral_component.hpp"
#include "integral.hpp"
#include "one_center.hpp"
#include "one_center_component.hpp"
#include "unique_map.hpp"

using T1CPair = OneCenterComponent;

using T2CIntegral = IntegralComponent<T1CPair, T1CPair>;

using VT2CIntegrals = VIntegralComponents<T1CPair, T1CPair>;

using ST2CIntegrals = SIntegralComponents<T1CPair, T1CPair>;

using R2CTerm = RecursionTerm<T2CIntegral>;

using V2CTerms = VRecursionTerms<T2CIntegral>;

using R2CDist = RecursionExpansion<T2CIntegral>;

using R2Group = RecursionGroup<T2CIntegral>;

using V2Groups = VRecursionGroups<T2CIntegral>;

using R2GroupContainer = Container<R2Group>;

using V2GroupContainers = VContainers<R2Group>; 

using I1CPair = OneCenter;

using I2CIntegral = Integral<I1CPair, I1CPair>;

using VI2CIntegrals = VIntegrals<I1CPair, I1CPair>;

using SI2CIntegrals = SIntegrals<I1CPair, I1CPair>;

using R2CMap = UniqueMap<I2CIntegral, T2CIntegral>; 

#endif /* t2c_defs_hpp */
