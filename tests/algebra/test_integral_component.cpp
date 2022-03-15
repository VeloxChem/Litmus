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

#include "test_integral_component.hpp"

#include "integral_component.hpp"
#include "two_center_pair_component.hpp"
#include "setters.hpp"

using T2CPair = TwoCenterPairComponent;

using T4CIntegral = IntegralComponent<T2CPair, T2CPair>;

TEST_F(IntegralComponentTest, Constructor)
{
    auto lhsint = T4CIntegral();
    
    auto bpair = T2CPair();
    
    auto kpair = T2CPair();
    
    auto operi = OperatorComponent();
    
    auto rhsint = T4CIntegral(bpair, kpair, operi, 0, {});
    
    EXPECT_EQ(lhsint, rhsint);
    
    operi = OperatorComponent("1/|r-r'|");
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    bpair = T2CPair({"GA", "GB"}, {p_x, f_yzz});
    
    kpair = T2CPair({"GC", "GD"}, {s_0, d_xy});
    
    lhsint = T4CIntegral(bpair, kpair, operi, 0, {});
    
    rhsint = T4CIntegral(bpair, kpair, operi);
    
    EXPECT_EQ(lhsint, rhsint);
    
    lhsint = T4CIntegral(bpair, kpair, operi, 2, {});
    
    rhsint = T4CIntegral(bpair, kpair, operi, 2);
    
    EXPECT_EQ(lhsint, rhsint);
}

TEST_F(IntegralComponentTest, OperatorBracket)
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
    
    EXPECT_EQ(t4cint[0], p_x);
    
    EXPECT_EQ(t4cint[1], f_yzz);
    
    EXPECT_EQ(t4cint[2], s_0);
    
    EXPECT_EQ(t4cint[3], d_xy);
}

TEST_F(IntegralComponentTest, OperatorEqual)
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

    const auto lhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    const auto rhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
        
    EXPECT_TRUE(lhsint == rhsint);
}

TEST_F(IntegralComponentTest, OperatorNotEqual)
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

    const auto lhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    bpair = T2CPair({"GB", "GB"}, {p_x, f_yzz});
    
    auto rhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != rhsint);
    
    bpair = T2CPair({"GA", "GB"}, {p_x, p_x});
    
    rhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != rhsint);
    
    bpair = T2CPair({"GA", "GB"}, {p_x, f_yzz});
    
    kpair = T2CPair({"GC", "LA"}, {s_0, d_xy});
    
    rhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != rhsint);
    
    kpair = T2CPair({"GC", "GD"}, {p_x, d_xy});
    
    rhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != rhsint);
    
    kpair = T2CPair({"GC", "GD"}, {s_0, d_xy});
    
    rhsint = T4CIntegral(bpair, kpair, opddr, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != rhsint);
    
    rhsint = T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != rhsint);
    
    rhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr});
    
    EXPECT_TRUE(lhsint != rhsint);
}

TEST_F(IntegralComponentTest, OperatorLess)
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

    const auto lhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_FALSE(lhsint < lhsint);
    
    bpair = T2CPair({"GB", "GB"}, {p_x, f_yzz});
    
    auto rhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint < rhsint);
    
    bpair = T2CPair({"GA", "GB"}, {p_x, p_x});
    
    rhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint < rhsint);
    
    bpair = T2CPair({"GA", "GB"}, {p_x, f_yzz});
    
    kpair = T2CPair({"GC", "LA"}, {s_0, d_xy});
    
    rhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint < rhsint);
    
    kpair = T2CPair({"GC", "GD"}, {p_x, d_xy});
    
    rhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint < rhsint);
    
    kpair = T2CPair({"GC", "GD"}, {s_0, d_xy});
    
    rhsint = T4CIntegral(bpair, kpair, opddr, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint < rhsint);
    
    rhsint = T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc});
    
    EXPECT_FALSE(lhsint < rhsint);
    
    rhsint = T4CIntegral(bpair, kpair, operi, 2, {opddr});
    
    EXPECT_FALSE(lhsint < rhsint);
}

TEST_F(IntegralComponentTest, Bra)
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
    
    EXPECT_EQ(t4cint.bra(), bpair);
}

TEST_F(IntegralComponentTest, Ket)
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
    
    EXPECT_EQ(t4cint.ket(), kpair);
}

TEST_F(IntegralComponentTest, Integrand)
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
    
    EXPECT_EQ(t4cint.integrand(), operi);
}

TEST_F(IntegralComponentTest, Order)
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
    
    EXPECT_EQ(t4cint.order(), 2);
}

TEST_F(IntegralComponentTest, Prefixes)
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
    
    EXPECT_EQ(t4cint.prefixes(), VOperatorComponents({opddr, opddc}));
}

TEST_F(IntegralComponentTest, Label)
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

    auto t4cint = T4CIntegral(bpair, kpair, operi);
    
    EXPECT_EQ(t4cint.label(), "x_yzz_0_xy");
    
    EXPECT_EQ(t4cint.label(true), "x_yzz_0_xy_0");
    
    t4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_EQ(t4cint.label(), "y_x_x_yzz_0_xy");
    
    EXPECT_EQ(t4cint.label(true), "y_x_x_yzz_0_xy_2");
    
    t4cint = T4CIntegral(bpair, kpair, opddr, 2, {opddr, opddc});
     
    EXPECT_EQ(t4cint.label(), "y_x_y_x_yzz_0_xy");
    
    EXPECT_EQ(t4cint.label(true), "y_x_y_x_yzz_0_xy_2");
}

TEST_F(IntegralComponentTest, Replace)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    auto bpair = T2CPair({"GA", "GB"}, {p_x, f_yzz});
    
    auto kpair = T2CPair({"GC", "GD"}, {s_0, d_xy});

    const auto t4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc});
    
    const auto r4cint = T4CIntegral(bpair, kpair, opddr, 2, {opddr, opddc});
    
    EXPECT_EQ(t4cint.replace(opddr), r4cint);
}

TEST_F(IntegralComponentTest, Shift)
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
    
    bpair = T2CPair({"GA", "GB"}, {s_0, f_yzz});
    
    auto r4cint = T4CIntegral(bpair, kpair, operi);
    
    EXPECT_EQ(t4cint.shift('x', -1, 0), r4cint);
    
    bpair = T2CPair({"GA", "GB"}, {p_x, d_yz});
    
    r4cint = T4CIntegral(bpair, kpair, operi);
    
    EXPECT_EQ(t4cint.shift('z', -1, 1), r4cint);
    
    bpair = T2CPair({"GA", "GB"}, {p_x, d_zz});
    
    r4cint = T4CIntegral(bpair, kpair, operi);
    
    EXPECT_EQ(t4cint.shift('y', -1, 1), r4cint);
    
    bpair = T2CPair({"GA", "GB"}, {p_x, f_yzz});
    
    kpair = T2CPair({"GC", "GD"}, {s_0, p_x});
    
    r4cint = T4CIntegral(bpair, kpair, operi);
    
    EXPECT_EQ(t4cint.shift('y', -1, 3), r4cint);
    
    kpair = T2CPair({"GC", "GD"}, {s_0, p_y});
    
    r4cint = T4CIntegral(bpair, kpair, operi);
    
    EXPECT_EQ(t4cint.shift('x', -1, 3), r4cint);
    
    EXPECT_FALSE(t4cint.shift('x', -2, 0));
    
    EXPECT_FALSE(t4cint.shift('y', -1, 0));
    
    EXPECT_FALSE(t4cint.shift('z', -1, 0));
    
    EXPECT_FALSE(t4cint.shift('x', -1, 1));
    
    EXPECT_FALSE(t4cint.shift('y', -2, 1));
    
    EXPECT_FALSE(t4cint.shift('z', -3, 1));
    
    EXPECT_FALSE(t4cint.shift('x', -1, 2));
    
    EXPECT_FALSE(t4cint.shift('y', -1, 2));
    
    EXPECT_FALSE(t4cint.shift('z', -1, 2));
    
    EXPECT_FALSE(t4cint.shift('x', -2, 3));
    
    EXPECT_FALSE(t4cint.shift('y', -2, 3));
    
    EXPECT_FALSE(t4cint.shift('z', -1, 3));
}

TEST_F(IntegralComponentTest, ShiftPrefix)
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
    
    auto r4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr0, opddc});
    
    EXPECT_EQ(t4cint.shift_prefix('y', -1, 0), r4cint);
    
    r4cint = T4CIntegral(bpair, kpair, operi, 2, {opddc});
    
    EXPECT_EQ(t4cint.shift_prefix('y', -1, 0, true), r4cint);
    
    r4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr, opddc0});
    
    EXPECT_EQ(t4cint.shift_prefix('x', -1, 1), r4cint);
    
    r4cint = T4CIntegral(bpair, kpair, operi, 2, {opddr});
    
    EXPECT_EQ(t4cint.shift_prefix('x', -1, 1, true), r4cint);
    
    EXPECT_FALSE(t4cint.shift_prefix('x', -1, 0));
    
    EXPECT_FALSE(t4cint.shift_prefix('y', -2, 0));
    
    EXPECT_FALSE(t4cint.shift_prefix('z', -1, 0));
    
    EXPECT_FALSE(t4cint.shift_prefix('x', -1, 0, true));
    
    EXPECT_FALSE(t4cint.shift_prefix('y', -2, 0, true));
    
    EXPECT_FALSE(t4cint.shift_prefix('z', -1, 0, true));
    
    EXPECT_FALSE(t4cint.shift_prefix('x', -2, 1));
    
    EXPECT_FALSE(t4cint.shift_prefix('y', -1, 1));
    
    EXPECT_FALSE(t4cint.shift_prefix('z', -1, 1));
    
    EXPECT_FALSE(t4cint.shift_prefix('x', -2, 1, true));
    
    EXPECT_FALSE(t4cint.shift_prefix('y', -1, 1, true));
    
    EXPECT_FALSE(t4cint.shift_prefix('z', -1, 1, true));
}

