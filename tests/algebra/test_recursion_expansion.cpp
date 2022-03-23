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

#include "test_recursion_expansion.hpp"

#include "recursion_expansion.hpp"
#include "recursion_term.hpp"
#include "integral_component.hpp"
#include "two_center_pair_component.hpp"
#include "setters.hpp"

using T2CPair = TwoCenterPairComponent;

using T4CIntegral = IntegralComponent<T2CPair, T2CPair>;

using R4CTerm = RecursionTerm<T4CIntegral>;

using R4CDist = RecursionExpansion<T4CIntegral>;

using S4CInts = SIntegralComponents<T2CPair, T2CPair>;

TEST_F(RecursionExpansionTest, Constructor)
{
    EXPECT_EQ(R4CDist(), R4CDist(R4CTerm(), {}));
    
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
    
    EXPECT_EQ(R4CDist(t4crt), R4CDist(t4crt, {}));
}

TEST_F(RecursionExpansionTest, OperatorBrackets)
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
    
    EXPECT_EQ(t4cdist[0], r4crta);
    
    EXPECT_EQ(t4cdist[1], r4crtb);
}

TEST_F(RecursionExpansionTest, OperatorEqual)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    const auto bpair = T2CPair({"GA", "GB"}, {p_x, f_yzz});
    
    const auto kpair = T2CPair({"GC", "GD"}, {s_0, d_xy});

    const auto t4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    const auto pbx = Factor("(P-B)", "pb", p_x);
    
    const auto wpy = Factor("(W-P)", "wp", p_y);
    
    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},});
    
    EXPECT_TRUE(R4CDist(t4crt) == R4CDist(t4crt, {}));
    
    const auto r4crta = R4CTerm(t4cint, {{pbx, 1},}, Fraction(1, 3));
    
    const auto r4cint = T4CIntegral(bpair, kpair, opddr, 2, {opddr, opddc});
    
    const auto r4crtb = R4CTerm(t4cint, {{wpy, 3},}, Fraction(1, 3));

    EXPECT_TRUE(R4CDist(t4crt, {r4crta, r4crtb}) == R4CDist(t4crt, {r4crta, r4crtb}));
}

TEST_F(RecursionExpansionTest, OperatorNotEqual)
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
    
    EXPECT_TRUE(t4cdist != R4CDist(t4crt));
    
    EXPECT_TRUE(t4cdist != R4CDist(r4crta, {r4crta, r4crtb}));
    
    EXPECT_TRUE(t4cdist != R4CDist(t4crt, {r4crta}));
}

TEST_F(RecursionExpansionTest, OperatorLess)
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
    
    EXPECT_FALSE(t4cdist < t4cdist);
    
    EXPECT_FALSE(t4cdist < R4CDist(t4crt));
    
    EXPECT_TRUE(t4cdist < R4CDist(t4crt, {r4crta, r4crtb, r4crtb}));
    
    EXPECT_FALSE(t4cdist < R4CDist(t4crt, {r4crta}));
}

TEST_F(RecursionExpansionTest, Add)
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
    
    EXPECT_EQ(t4cdist, R4CDist(t4crt));
    
    t4cdist.add(r4crta);
    
    EXPECT_EQ(t4cdist, R4CDist(t4crt, {r4crta,}));
    
    t4cdist.add(r4crtb);
    
    EXPECT_EQ(t4cdist, R4CDist(t4crt, {r4crta, r4crtb}));
}

TEST_F(RecursionExpansionTest, Root)
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
    
    EXPECT_EQ(t4cdist.root(), t4crt);
}

TEST_F(RecursionExpansionTest, Terms)
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
    
    EXPECT_EQ(t4cdist.terms(), 2);
}

TEST_F(RecursionExpansionTest, UniqueIntegrals)
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
    
    auto r4cint = T4CIntegral(bpair, kpair, opddr, 1, {opddr, opddc});
    
    const auto pbx = Factor("(P-B)", "pb", p_x);
    
    const auto wpy = Factor("(W-P)", "wp", p_y);
    
    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(3, 7));
    
    const auto r4crta = R4CTerm(r4cint, {{pbx, 1},}, Fraction(1, 3));
    
    const auto r4crtb = R4CTerm(t4cint, {{wpy, 3},}, Fraction(1, 3));
    
    auto t4cdist = R4CDist(t4crt, {r4crta, r4crtb});
    
    EXPECT_EQ(t4cdist.unique_integrals(), S4CInts({t4cint, r4cint}));
    
    t4cdist = R4CDist(r4crta, {t4crt, r4crtb});
    
    EXPECT_EQ(t4cdist.unique_integrals(), S4CInts({t4cint}));
}

TEST_F(RecursionExpansionTest, CountNewIntegrals)
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
    
    auto r4cint = T4CIntegral(bpair, kpair, opddr, 1, {opddr, opddc});
    
    const auto pbx = Factor("(P-B)", "pb", p_x);
    
    const auto wpy = Factor("(W-P)", "wp", p_y);
    
    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(3, 7));
    
    const auto r4crta = R4CTerm(r4cint, {{pbx, 1},}, Fraction(1, 3));
    
    const auto r4crtb = R4CTerm(t4cint, {{wpy, 3},}, Fraction(1, 3));
    
    auto t4cdist = R4CDist(t4crt, {r4crta, r4crtb});
    
    EXPECT_EQ(t4cdist.count_new_integrals({}), 2);
    
    EXPECT_EQ(t4cdist.count_new_integrals({t4cint,}), 1);
    
    EXPECT_EQ(t4cdist.count_new_integrals({r4cint,}), 1);
    
    EXPECT_EQ(t4cdist.count_new_integrals({r4cint, t4cint}), 0);
    
    EXPECT_EQ(t4cdist.count_new_integrals({t4cint, r4cint}), 0);
    
    t4cdist = R4CDist(r4crta, {t4crt, r4crtb});
    
    EXPECT_EQ(t4cdist.count_new_integrals({}), 1);
    
    EXPECT_EQ(t4cdist.count_new_integrals({t4cint,}), 0);
    
    EXPECT_EQ(t4cdist.count_new_integrals({r4cint,}), 1);
}

TEST_F(RecursionExpansionTest, Auxilary)
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
    
    EXPECT_TRUE(t4cdist.auxilary(2));
    
    EXPECT_FALSE(t4cdist.auxilary(0));
    
    EXPECT_FALSE(t4cdist.auxilary(1));
    
    EXPECT_FALSE(t4cdist.auxilary(3));
}
