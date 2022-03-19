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

#include "test_recursion_term.hpp"

#include "recursion_term.hpp"
#include "integral_component.hpp"
#include "two_center_pair_component.hpp"
#include "setters.hpp"

using T2CPair = TwoCenterPairComponent;

using T4CIntegral = IntegralComponent<T2CPair, T2CPair>;

using R4CTerm = RecursionTerm<T4CIntegral>;

TEST_F(RecursionTermTest, Constructor)
{
    auto t4cint = T4CIntegral();
    
    const auto fact = Fraction(1);
    
    EXPECT_EQ(R4CTerm(t4cint, {}, fact), R4CTerm());
    
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
    
    t4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    const auto pbx = Factor("(P-B)", "pb", p_x);
    
    const auto wpy = Factor("(W-P)", "wp", p_y);
    
    EXPECT_EQ(R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}),
              R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1)));
}

TEST_F(RecursionTermTest, OperatorBracket)
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

    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt[0], p_x);
    
    EXPECT_EQ(t4crt[1], f_yzz);
    
    EXPECT_EQ(t4crt[2], s_0);
    
    EXPECT_EQ(t4crt[3], d_xy);
}

TEST_F(RecursionTermTest, OperatorEqual)
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
    
    const auto lhsrt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},});
    
    const auto rhsrt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1));
    
    EXPECT_TRUE(lhsrt == rhsrt);
}

TEST_F(RecursionTermTest, OperatorNotEqual)
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
    
    const auto lhsrt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(3, 7));
    
    t4cint = T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc});
    
    auto rhsrt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(3, 7));
    
    EXPECT_TRUE(lhsrt != rhsrt);
    
    t4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    rhsrt = R4CTerm(t4cint, {{pbx, 2}, {wpy, 2},}, Fraction(3, 7));
    
    EXPECT_TRUE(lhsrt != rhsrt);
    
    rhsrt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(3, 5));
    
    EXPECT_TRUE(lhsrt != rhsrt);
}

TEST_F(RecursionTermTest, OperatorLess)
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
    
    const auto lhsrt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(3, 7));
    
    EXPECT_FALSE(lhsrt < lhsrt);
    
    t4cint = T4CIntegral(bpair, kpair, operi, 5, {opddr, opddc});
    
    auto rhsrt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(3, 7));
    
    EXPECT_TRUE(lhsrt < rhsrt);
    
    t4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    rhsrt = R4CTerm(t4cint, {{pbx, 2}, {wpy, 2},}, Fraction(3, 7));
    
    EXPECT_TRUE(lhsrt < rhsrt);
    
    rhsrt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(3, 5));
    
    EXPECT_TRUE(lhsrt < rhsrt);
}

TEST_F(RecursionTermTest, Bra)
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

    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.bra<T2CPair>(), bpair);
}

TEST_F(RecursionTermTest, Ket)
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

    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.ket<T2CPair>(), kpair);
}

TEST_F(RecursionTermTest, Order)
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

    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.order(), 2);
}

TEST_F(RecursionTermTest, Prefixes)
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

    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.prefixes(), VOperatorComponents({opddr, opddc}));
}

TEST_F(RecursionTermTest, Integral)
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

    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.integral(), t4cint);
}

TEST_F(RecursionTermTest, Prefactor)
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

    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.prefactor(), Fraction(1, 3));
}

TEST_F(RecursionTermTest, Factors)
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

    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.factors(), std::set<Factor>({pbx, wpy}));
}

TEST_F(RecursionTermTest, FactorOrder)
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

    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.factor_order(pbx), 1);
    
    EXPECT_EQ(t4crt.factor_order(wpy), 2);
    
    const auto eta = Factor("1/eta", "fxi", TensorComponent(0, 0, 0));
    
    EXPECT_EQ(t4crt.factor_order(eta), 0);
}

TEST_F(RecursionTermTest, Label)
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

    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.label(), "y_x_x_yzz_0_xy");
    
    EXPECT_EQ(t4crt.label(true), "y_x_x_yzz_0_xy_2");
}

TEST_F(RecursionTermTest, Replace)
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

    auto t4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    const auto pbx = Factor("(P-B)", "pb", p_x);
    
    const auto wpy = Factor("(W-P)", "wp", p_y);

    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    t4cint = T4CIntegral(bpair, kpair, opddr, 2, {opddr, opddc});
    
    const auto r4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.replace(opddr), r4crt);
}

TEST_F(RecursionTermTest, Shift)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yz = TensorComponent(0, 1, 1);
    
    const auto d_zz = TensorComponent(0, 0, 2);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    auto bpair = T2CPair({"GA", "GB"}, {p_x, f_yzz});
    
    auto kpair = T2CPair({"GC", "GD"}, {s_0, d_xy});

    const auto t4cint = T4CIntegral(bpair, kpair, operi);
    
    const auto pbx = Factor("(P-B)", "pb", p_x);
    
    const auto wpy = Factor("(W-P)", "wp", p_y);

    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    bpair = T2CPair({"GA", "GB"}, {s_0, f_yzz});
    
    auto r4cint = T4CIntegral(bpair, kpair, operi);
    
    auto r4crt = R4CTerm(r4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.shift('x', -1, 0), r4crt);
    
    bpair = T2CPair({"GA", "GB"}, {p_x, d_yz});
    
    r4cint = T4CIntegral(bpair, kpair, operi);
    
    r4crt = R4CTerm(r4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.shift('z', -1, 1), r4crt);
    
    bpair = T2CPair({"GA", "GB"}, {p_x, d_zz});
    
    r4cint = T4CIntegral(bpair, kpair, operi);
    
    r4crt = R4CTerm(r4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.shift('y', -1, 1), r4crt);
    
    bpair = T2CPair({"GA", "GB"}, {p_x, f_yzz});
    
    kpair = T2CPair({"GC", "GD"}, {s_0, p_x});
    
    r4cint = T4CIntegral(bpair, kpair, operi);
    
    r4crt = R4CTerm(r4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.shift('y', -1, 3), r4crt);
    
    kpair = T2CPair({"GC", "GD"}, {s_0, p_y});
    
    r4cint = T4CIntegral(bpair, kpair, operi);
    
    r4crt = R4CTerm(r4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.shift('x', -1, 3), r4crt);
    
    EXPECT_FALSE(t4crt.shift('x', -2, 0));
    
    EXPECT_FALSE(t4crt.shift('y', -1, 0));
    
    EXPECT_FALSE(t4crt.shift('z', -1, 0));
    
    EXPECT_FALSE(t4crt.shift('x', -1, 1));
    
    EXPECT_FALSE(t4crt.shift('y', -2, 1));
    
    EXPECT_FALSE(t4crt.shift('z', -3, 1));
    
    EXPECT_FALSE(t4crt.shift('x', -1, 2));
    
    EXPECT_FALSE(t4crt.shift('y', -1, 2));
    
    EXPECT_FALSE(t4crt.shift('z', -1, 2));
    
    EXPECT_FALSE(t4crt.shift('x', -2, 3));
    
    EXPECT_FALSE(t4crt.shift('y', -2, 3));
    
    EXPECT_FALSE(t4crt.shift('z', -1, 3));
}

TEST_F(RecursionTermTest, ShiftPrefix)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    const auto opddr0 = OperatorComponent("d/dr", s_0, "bra", 1);
    
    const auto opddc0 = OperatorComponent("d/dC", s_0, "ket", 0);
    
    auto bpair = T2CPair({"GA", "GB"}, {p_x, f_yzz});
    
    auto kpair = T2CPair({"GC", "GD"}, {s_0, d_xy});

    const auto t4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    const auto pbx = Factor("(P-B)", "pb", p_x);
    
    const auto wpy = Factor("(W-P)", "wp", p_y);
    
    const auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    auto r4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr0, opddc});
    
    auto r4crt = R4CTerm(r4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.shift_prefix('y', -1, 0), r4crt);
    
    r4cint = T4CIntegral(bpair, kpair, operi, 2, {opddc});
    
    r4crt = R4CTerm(r4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.shift_prefix('y', -1, 0, true), r4crt);
    
    r4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc0});
    
    r4crt = R4CTerm(r4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.shift_prefix('x', -1, 1), r4crt);
    
    r4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr});
    
    r4crt = R4CTerm(r4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    EXPECT_EQ(t4crt.shift_prefix('x', -1, 1, true), r4crt);
    
    EXPECT_FALSE(t4crt.shift_prefix('x', -1, 0));
    
    EXPECT_FALSE(t4crt.shift_prefix('y', -2, 0));
    
    EXPECT_FALSE(t4crt.shift_prefix('z', -1, 0));
    
    EXPECT_FALSE(t4crt.shift_prefix('x', -1, 0, true));
    
    EXPECT_FALSE(t4crt.shift_prefix('y', -2, 0, true));
    
    EXPECT_FALSE(t4crt.shift_prefix('z', -1, 0, true));
    
    EXPECT_FALSE(t4crt.shift_prefix('x', -2, 1));
    
    EXPECT_FALSE(t4crt.shift_prefix('y', -1, 1));
    
    EXPECT_FALSE(t4crt.shift_prefix('z', -1, 1));
    
    EXPECT_FALSE(t4crt.shift_prefix('x', -2, 1, true));
    
    EXPECT_FALSE(t4crt.shift_prefix('y', -1, 1, true));
    
    EXPECT_FALSE(t4crt.shift_prefix('z', -1, 1, true));
}

TEST_F(RecursionTermTest, Add)
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

    auto t4crt = R4CTerm(t4cint);
    
    auto r4crt = R4CTerm(t4cint, {{pbx, 1}}, Fraction(1, 2));
    
    t4crt.add(pbx, Fraction(1, 2));
    
    r4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 1}}, Fraction(5, 2));
    
    t4crt.add(wpy, Fraction(5));
    
    EXPECT_EQ(t4crt, r4crt);
    
    r4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2}}, Fraction(5, 2));
    
    t4crt.add(wpy);
    
    EXPECT_EQ(t4crt, r4crt);
}

TEST_F(RecursionTermTest, Scale)
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

    auto t4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 3));
    
    const auto r4crt = R4CTerm(t4cint, {{pbx, 1}, {wpy, 2},}, Fraction(1, 2));
    
    t4crt.scale(Fraction(3, 2));
    
    EXPECT_EQ(t4crt, r4crt);
}
