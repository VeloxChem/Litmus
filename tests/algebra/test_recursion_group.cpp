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

#include "test_recursion_group.hpp"

#include "recursion_group.hpp"
#include "recursion_expansion.hpp"
#include "recursion_term.hpp"
#include "integral_component.hpp"
#include "two_center_pair_component.hpp"
#include "integral.hpp"
#include "two_center_pair.hpp"
#include "setters.hpp"

using T2CPair = TwoCenterPairComponent;

using T4CIntegral = IntegralComponent<T2CPair, T2CPair>;

using R4CTerm = RecursionTerm<T4CIntegral>;

using R4CDist = RecursionExpansion<T4CIntegral>;

using R4Group = RecursionGroup<T4CIntegral>;

using I2CPair = TwoCenterPair;

using I4CIntegral = Integral<I2CPair, I2CPair>;

TEST_F(RecursionGroupTest, Constructor)
{
    EXPECT_EQ(R4Group(), R4Group(VRecursionExpansions<T4CIntegral>({})));
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    const auto bpair = T2CPair({"GA", "GB"}, {p_x, f_yzz});
    
    const auto kpair = T2CPair({"GC", "GD"}, {s_0, d_xy});
    
    auto t4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    const auto pbx = Factor("(P-B)", "pb", p_x);
    
    const auto wpy = Factor("(W-P)", "wp", p_y);
    
    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    const auto r4crta = R4CTerm(t4cint, {{pbx, 1},}, Fraction(1, 3));
    
    const auto r4crtb = R4CTerm(t4cint, {{wpy, 3},}, Fraction(1, 3));
    
    const auto t4cdist = R4CDist(t4crt, {r4crta, r4crtb});
    
    const auto r4cdist = R4CDist(r4crta, {r4crtb});
    
    EXPECT_EQ(R4Group({t4crt, r4cdist}), R4Group({t4crt, r4cdist}));
}

TEST_F(RecursionGroupTest, OperatorBrackets)
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
    
    EXPECT_EQ(t4group[0], t4cdist);
    
    EXPECT_EQ(t4group[1], r4cdist);
}

TEST_F(RecursionGroupTest, OperatorEqual)
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
    
    EXPECT_TRUE(R4Group({t4cdist, r4cdist}) == R4Group({t4cdist, r4cdist}));
}

TEST_F(RecursionGroupTest, OperatorNotEqual)
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
    
    EXPECT_TRUE(R4Group({t4cdist, r4cdist}) != R4Group({r4cdist, r4cdist}));
    
    EXPECT_TRUE(R4Group({t4cdist, r4cdist}) != R4Group({t4cdist}));
}

TEST_F(RecursionGroupTest, OperatorLess)
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
    
    EXPECT_TRUE(R4Group({t4cdist, r4cdist}) < R4Group({t4cdist, t4cdist}));
    
    EXPECT_FALSE(R4Group({t4cdist, r4cdist}) < R4Group({t4cdist}));
}

TEST_F(RecursionGroupTest, Add)
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
    
    auto t4group = R4Group();
    
    EXPECT_EQ(t4group, R4Group(VRecursionExpansions<T4CIntegral>({})));
 
    t4group.add(r4cdist);
    
    EXPECT_EQ(t4group, R4Group({r4cdist}));
    
    t4group.add(t4cdist);
    
    EXPECT_EQ(t4group, R4Group({r4cdist, t4cdist}));
}

TEST_F(RecursionGroupTest, Expansions)
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
    
    EXPECT_EQ(t4group.expansions(), 2);
}

TEST_F(RecursionGroupTest, SplitTerms)
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
    
    const auto t4group = R4Group({t4cdist, r4cdist});

    auto mterms = t4group.split_terms<I4CIntegral>();
    
    EXPECT_EQ(mterms.size(), 2);
    
    EXPECT_EQ(mterms[0], VRecursionTerms<T4CIntegral>({R4CTerm(r4cint),}));
    
    EXPECT_EQ(mterms[1], VRecursionTerms<T4CIntegral>({R4CTerm(t4cint),}));
}

TEST_F(RecursionGroupTest, Roots)
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
    
    const auto t4group = R4Group({t4cdist, r4cdist});

    auto vterms = t4group.roots();
    
    EXPECT_EQ(vterms.size(), 2);
    
    EXPECT_EQ(vterms[0], R4CTerm(t4cint));
    
    EXPECT_EQ(vterms[1], R4CTerm(r4cint));
}

TEST_F(RecursionGroupTest, Empty)
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
    
    auto t4cdist = R4CDist(t4crt);
    
    auto r4cdist = R4CDist(r4crta, {r4crtb,});
    
    auto t4group = R4Group({t4cdist, r4cdist});
    
    EXPECT_FALSE(t4group.empty());
    
    r4cdist = R4CDist(r4crta);
    
    t4group = R4Group({t4cdist, r4cdist});
    
    EXPECT_TRUE(t4group.empty());
}
