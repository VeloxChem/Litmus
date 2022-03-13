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

#include "test_four_center_integral_component.hpp"

#include "four_center_integral_component.hpp"
#include "setters.hpp"

TEST_F(FourCenterIntegralComponentTest, Constructor)
{
    auto lhsint = FourCenterIntegralComponent();
    
    auto bpair = TwoCenterPairComponent();
    
    auto kpair = TwoCenterPairComponent();
    
    auto operi = OperatorComponent();
    
    auto rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 0, {});
    
    EXPECT_EQ(lhsint, rhsint);
    
    operi = OperatorComponent("1/|r-r'|");
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy});
    
    lhsint = FourCenterIntegralComponent(bpair, kpair, operi, 0, {});
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, operi);
    
    EXPECT_EQ(lhsint, rhsint);
    
    lhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {});
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2);
    
    EXPECT_EQ(lhsint, rhsint);
}

TEST_F(FourCenterIntegralComponentTest, OperatorEqual)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    const auto bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    const auto kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy});

    const auto lhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    const auto rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
        
    EXPECT_TRUE(lhsint == rhsint);
}

TEST_F(FourCenterIntegralComponentTest, OperatorNotEqual)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    auto bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    auto kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy});

    const auto lhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    bpair = TwoCenterPairComponent({"GB", "GB"}, {p_x, f_yzz});
    
    auto rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != rhsint);
    
    bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, p_x});
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != rhsint);
    
    bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    kpair = TwoCenterPairComponent({"GC", "LA"}, {s_0, d_xy});
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != rhsint);
    
    kpair = TwoCenterPairComponent({"GC", "GD"}, {p_x, d_xy});
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != rhsint);
    
    kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy});
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, opddr, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != rhsint);
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 1, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != rhsint);
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr});
    
    EXPECT_TRUE(lhsint != rhsint);
}

TEST_F(FourCenterIntegralComponentTest, OperatorLess)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    auto bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    auto kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy});

    const auto lhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_FALSE(lhsint < lhsint);
    
    bpair = TwoCenterPairComponent({"GB", "GB"}, {p_x, f_yzz});
    
    auto rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint < rhsint);
    
    bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, p_x});
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint < rhsint);
    
    bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    kpair = TwoCenterPairComponent({"GC", "LA"}, {s_0, d_xy});
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint < rhsint);
    
    kpair = TwoCenterPairComponent({"GC", "GD"}, {p_x, d_xy});
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint < rhsint);
    
    kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy});
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, opddr, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint < rhsint);
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 1, {opddr, opddc});
    
    EXPECT_FALSE(lhsint < rhsint);
    
    rhsint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr});
    
    EXPECT_FALSE(lhsint < rhsint);
}

TEST_F(FourCenterIntegralComponentTest, BraPair)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    const auto bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    const auto kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy});

    const auto t4cint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_EQ(t4cint.bra_pair(), bpair);
}

TEST_F(FourCenterIntegralComponentTest, KetPair)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    const auto bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    const auto kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy});

    const auto t4cint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_EQ(t4cint.ket_pair(), kpair);
}

TEST_F(FourCenterIntegralComponentTest, Integrand)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    const auto bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    const auto kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy});

    const auto t4cint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_EQ(t4cint.integrand(), operi);
}

TEST_F(FourCenterIntegralComponentTest, Order)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    const auto bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    const auto kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy});

    const auto t4cint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_EQ(t4cint.order(), 2);
}

TEST_F(FourCenterIntegralComponentTest, Prefixes)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    const auto bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    const auto kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy});

    const auto t4cint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_EQ(t4cint.prefixes(), VOperatorComponents({opddr, opddc}));
}

TEST_F(FourCenterIntegralComponentTest, Label)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto opddr = OperatorComponent("d/dr", p_y, "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", p_x, "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    const auto bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    const auto kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy});

    auto t4cint = FourCenterIntegralComponent(bpair, kpair, operi);
    
    EXPECT_EQ(t4cint.label(), "x_yzz_0_xy");
    
    EXPECT_EQ(t4cint.label(true), "x_yzz_0_xy_0");
    
    t4cint = FourCenterIntegralComponent(bpair, kpair, operi, 2, {opddr, opddc});
    
    EXPECT_EQ(t4cint.label(), "y_x_x_yzz_0_xy");
    
    EXPECT_EQ(t4cint.label(true), "y_x_x_yzz_0_xy_2");
    
    t4cint = FourCenterIntegralComponent(bpair, kpair, opddr, 2, {opddr, opddc});
     
    EXPECT_EQ(t4cint.label(), "y_x_y_x_yzz_0_xy");
    
    EXPECT_EQ(t4cint.label(true), "y_x_y_x_yzz_0_xy_2");
}

TEST_F(FourCenterIntegralComponentTest, Shift)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yz = TensorComponent(0, 1, 1);
    
    const auto d_zz = TensorComponent(0, 0, 2);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    auto bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    auto kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy});

    const auto t4cint = FourCenterIntegralComponent(bpair, kpair, operi);
    
    bpair = TwoCenterPairComponent({"GA", "GB"}, {s_0, f_yzz});
    
    auto r4cint = FourCenterIntegralComponent(bpair, kpair, operi);
    
    EXPECT_EQ(t4cint.shift('x', -1, 0), r4cint);
    
    bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz});
    
    r4cint = FourCenterIntegralComponent(bpair, kpair, operi);
    
    EXPECT_EQ(t4cint.shift('z', -1, 1), r4cint);
    
    bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, d_zz});
    
    r4cint = FourCenterIntegralComponent(bpair, kpair, operi);
    
    EXPECT_EQ(t4cint.shift('y', -1, 1), r4cint);
    
    bpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz});
    
    kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, p_x});
    
    r4cint = FourCenterIntegralComponent(bpair, kpair, operi);
    
    EXPECT_EQ(t4cint.shift('y', -1, 3), r4cint);
    
    kpair = TwoCenterPairComponent({"GC", "GD"}, {s_0, p_y});
    
    r4cint = FourCenterIntegralComponent(bpair, kpair, operi);
    
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

