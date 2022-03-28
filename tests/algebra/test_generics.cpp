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

#include "test_generics.hpp"

#include <string>

#include "generics.hpp"
#include "recursion_group.hpp"
#include "recursion_expansion.hpp"
#include "recursion_term.hpp"
#include "integral_component.hpp"
#include "two_center_pair_component.hpp"
#include "integral.hpp"
#include "two_center_pair.hpp"

using T2CPair = TwoCenterPairComponent;

using T4CIntegral = IntegralComponent<T2CPair, T2CPair>;

using R4CTerm = RecursionTerm<T4CIntegral>;

using R4CDist = RecursionExpansion<T4CIntegral>;

using R4Group = RecursionGroup<T4CIntegral>;

using I2CPair = TwoCenterPair;

using I4CIntegral = Integral<I2CPair, I2CPair>;

TEST_F(GenericsTest, MergeForString)
{
    std::string lhs_str("A");
    
    const std::string rhs_str("B");
   
    gen::merge(lhs_str, rhs_str);
    
    EXPECT_EQ(lhs_str, "AB");
}

TEST_F(GenericsTest, MergeForRecursionGroup)
{
    const auto operi = OperatorComponent("1/|r-r'|");
        
    const auto p_x = TensorComponent(1, 0, 0);
        
    const auto p_y = TensorComponent(0, 1, 0);
        
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
        
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
        
    const auto s_0 = TensorComponent(0, 0, 0);
        
    const auto d_xy = TensorComponent(1, 1, 0);
        
    const auto f_yzz = TensorComponent(0, 1, 2);
        
    auto bpair = T2CPair({"GA", "GB"}, {p_x, f_yzz});
        
    auto kpair = T2CPair({"GC", "GD"}, {s_0, d_xy});

    auto t4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
        
    auto r4cint = T4CIntegral(bpair, kpair, operi, 1, {opddc});
        
    const auto pbx = Factor("(P-B)", "pb", p_x);
        
    const auto wpy = Factor("(W-P)", "wp", p_y);
        
    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(3, 7));
        
    const auto r4crta = R4CTerm(r4cint, {{pbx, 1},}, Fraction(1, 3));
        
    const auto r4crtb = R4CTerm(t4cint, {{wpy, 3},}, Fraction(1, 3));
        
    const auto t4cdist = R4CDist(t4crt, {r4crta, r4crtb});
        
    const auto r4cdist = R4CDist(r4crta, {r4crtb});
        
    auto t4group = R4Group({t4cdist});
        
    auto r4group = R4Group({t4cdist});
        
    gen::merge<R4Group>(t4group, r4group);
        
    EXPECT_EQ(t4group, R4Group({t4cdist}));
        
    r4group = R4Group({r4cdist});
        
    gen::merge<R4Group>(t4group, r4group);
        
    EXPECT_EQ(t4group, R4Group({t4cdist, r4cdist}));
        
    gen::merge<R4Group>(t4group, r4group);
        
    EXPECT_EQ(t4group, R4Group({t4cdist, r4cdist}));
}

TEST_F(GenericsTest, SimilarForString)
{
    const std::string lhs_str("A");
    
    const std::string rhs_str("B");
    
    EXPECT_TRUE(gen::similar(lhs_str, lhs_str));
    
    EXPECT_TRUE(gen::similar(rhs_str, rhs_str));
    
    EXPECT_FALSE(gen::similar(lhs_str, rhs_str));
    
    EXPECT_FALSE(gen::similar(rhs_str, lhs_str));
}

TEST_F(GenericsTest, SimilarForRecursionGroup)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    auto bpair = T2CPair({"GA", "GB"}, {p_x, f_yzz});
    
    auto kpair = T2CPair({"GC", "GD"}, {s_0, d_xy});

    auto t4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    const auto pbx = Factor("(P-B)", "pb", p_x);
    
    const auto wpy = Factor("(W-P)", "wp", p_y);
    
    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(3, 7));
    
    const auto r4crta = R4CTerm(t4cint, {{pbx, 1},}, Fraction(1, 3));
    
    const auto r4crtb = R4CTerm(t4cint, {{wpy, 3},}, Fraction(1, 3));
    
    const auto t4cdist = R4CDist(t4crt, {r4crta, r4crtb});
    
    const auto r4cdist = R4CDist(r4crta, {r4crtb});
    
    const auto t4group = R4Group({t4cdist, r4cdist});
    
    EXPECT_TRUE(gen::similar(t4group, R4Group({t4cdist, r4cdist})));
    
    EXPECT_TRUE(gen::similar(t4group, R4Group({t4cdist,})));

    EXPECT_TRUE(gen::similar(t4group, R4Group({r4cdist,})));
}
