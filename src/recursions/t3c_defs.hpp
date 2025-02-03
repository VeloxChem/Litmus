#ifndef t3c_defs_hpp
#define t3c_defs_hpp

#include <utility>

#include "recursion_group.hpp"
#include "recursion_expansion.hpp"
#include "recursion_term.hpp"
#include "integral_component.hpp"
#include "integral.hpp"
#include "one_center.hpp"
#include "one_center_component.hpp"
#include "two_center_pair_component.hpp"
#include "two_center_pair.hpp"

using T1CPair = OneCenterComponent;

using T2CPair = TwoCenterPairComponent;

using T3CIntegral = IntegralComponent<T1CPair, T2CPair>;

using VT3CIntegrals = VIntegralComponents<T1CPair, T2CPair>;

using ST3CIntegrals = SIntegralComponents<T1CPair, T2CPair>;

using R3CTerm = RecursionTerm<T3CIntegral>;

using V3CTerms = VRecursionTerms<T3CIntegral>;

using R3CDist = RecursionExpansion<T3CIntegral>;

using R3Group = RecursionGroup<T3CIntegral>;

using V3Groups = VRecursionGroups<T3CIntegral>;

using I1CPair = OneCenter;

using I2CPair = TwoCenterPair;

using I3CIntegral = Integral<I1CPair, I2CPair>;

using SI3CIntegrals = SIntegrals<I1CPair, I2CPair>;

using G3Term = std::pair<std::array<int, 3>, I3CIntegral>;

using SG3Terms = std::set<G3Term>;

#endif /* t3c_defs_hpp */
