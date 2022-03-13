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
    EXPECT_EQ(FourCenterIntegralComponent(),
              FourCenterIntegralComponent(TwoCenterPairComponent(),
                                          TwoCenterPairComponent(),
                                          OperatorComponent(), 0, {}));
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);

    EXPECT_EQ(FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                          TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                          operi, 0, {}),
              FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                          TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                          operi));
    
    EXPECT_EQ(FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                          TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                          operi, 2, {}),
              FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                          TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                          operi, 2));

}

TEST_F(FourCenterIntegralComponentTest, OperatorEqual)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto opddr = OperatorComponent("d/dr", TensorComponent(0, 1, 0), "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", TensorComponent(1, 0, 0), "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);

    const auto lhsint = FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                    TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                    operi, 2, {opddr, opddc});
        
    EXPECT_TRUE(lhsint == FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                      TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                      operi, 2, {opddr, opddc}));
}

TEST_F(FourCenterIntegralComponentTest, OperatorNotEqual)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto opddr = OperatorComponent("d/dr", TensorComponent(0, 1, 0), "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", TensorComponent(1, 0, 0), "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);

    const auto lhsint = FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                    TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                    operi, 2, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != FourCenterIntegralComponent(TwoCenterPairComponent({"GB", "GB"}, {p_x, f_yzz}),
                                                      TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                      operi, 2, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint != FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, p_x}),
                                                      TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                      operi, 2, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint != FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                      TwoCenterPairComponent({"GC", "LA"}, {s_0, d_xy}),
                                                      operi, 2, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint != FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                      TwoCenterPairComponent({"GC", "GD"}, {p_x, d_xy}),
                                                      operi, 2, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint != FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                      TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                      opddr, 2, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint != FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                      TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                      operi, 1, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint != FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                      TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                      operi, 2, {opddr}));
}

TEST_F(FourCenterIntegralComponentTest, OperatorLess)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto opddr = OperatorComponent("d/dr", TensorComponent(0, 1, 0), "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", TensorComponent(1, 0, 0), "ket", 0);
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);

    const auto lhsint = FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                    TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                    operi, 2, {opddr, opddc});
    
    EXPECT_FALSE(lhsint < lhsint);
    
    EXPECT_TRUE(lhsint < FourCenterIntegralComponent(TwoCenterPairComponent({"GB", "GB"}, {p_x, f_yzz}),
                                                      TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                      operi, 2, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint < FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, p_x}),
                                                      TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                      operi, 2, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint < FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                      TwoCenterPairComponent({"GC", "LA"}, {s_0, d_xy}),
                                                      operi, 2, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint < FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                      TwoCenterPairComponent({"GC", "GD"}, {p_x, d_xy}),
                                                      operi, 2, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint < FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                      TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                      opddr, 2, {opddr, opddc}));
    
    EXPECT_FALSE(lhsint < FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                      TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                      operi, 1, {opddr, opddc}));
    
    EXPECT_FALSE(lhsint < FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                                      TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                                      operi, 2, {opddr}));
}

TEST_F(FourCenterIntegralComponentTest, ToString)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);

    auto t4cint = FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                              TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                              operi);
    
    EXPECT_EQ(t4cint.to_string(),
              "{GA:(1,0,0);GB:(0,1,2)}{1/|r-r'|:(0,0,0)}[none:-1]{GC:(0,0,0);GD:(1,1,0)}^(0)");
    
    const auto opddr = OperatorComponent("d/dr", TensorComponent(0, 1, 0), "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", TensorComponent(1, 0, 0), "ket", 0);
    
    t4cint = FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                         TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                         operi, 2, {opddr, opddc});
    
    EXPECT_EQ(t4cint.to_string(),
              "[{d/dr:(0,1,0)}[bra:1];{d/dC:(1,0,0)}[ket:0];]{GA:(1,0,0);GB:(0,1,2)}{1/|r-r'|:(0,0,0)}[none:-1]{GC:(0,0,0);GD:(1,1,0)}^(2)");
}

TEST_F(FourCenterIntegralComponentTest, Label)
{
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto f_yzz = TensorComponent(0, 1, 2);

    auto t4cint = FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                              TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                              operi);
    
    EXPECT_EQ(t4cint.label(), "x_yzz_0_xy");
    
    EXPECT_EQ(t4cint.label(true), "x_yzz_0_xy_0");
    
    const auto opddr = OperatorComponent("d/dr", TensorComponent(0, 1, 0), "bra", 1);
    
    const auto opddc = OperatorComponent("d/dC", TensorComponent(1, 0, 0), "ket", 0);
    
    t4cint = FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                         TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                         operi, 2, {opddr, opddc});
    
    
    EXPECT_EQ(t4cint.label(), "y_x_x_yzz_0_xy");
    
    EXPECT_EQ(t4cint.label(true), "y_x_x_yzz_0_xy_2");
    
    
    t4cint = FourCenterIntegralComponent(TwoCenterPairComponent({"GA", "GB"}, {p_x, f_yzz}),
                                         TwoCenterPairComponent({"GC", "GD"}, {s_0, d_xy}),
                                         opddr, 2, {opddr, opddc});
    
    
    EXPECT_EQ(t4cint.label(), "y_x_y_x_yzz_0_xy");
    
    EXPECT_EQ(t4cint.label(true), "y_x_y_x_yzz_0_xy_2");
    
     
}
