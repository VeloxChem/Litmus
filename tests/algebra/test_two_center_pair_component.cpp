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

#include "test_two_center_pair_component.hpp"

#include "two_center_pair_component.hpp"

TEST_F(TwoCenterPairComponentTest, Constructor)
{
    const auto scomp = TensorComponent(0, 0, 0);
    
    EXPECT_EQ(TwoCenterPairComponent(), TwoCenterPairComponent({"", ""}, {scomp, scomp}));
}

TEST_F(TwoCenterPairComponentTest, OperatorBrackets)
{
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_yz = TensorComponent(0, 1, 1);
    
    const auto tpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz});
    
    EXPECT_EQ(tpair[0], p_x);
    
    EXPECT_EQ(tpair[1], d_yz);
}

TEST_F(TwoCenterPairComponentTest, OperatorEqual)
{
    const auto scomp = TensorComponent(0, 0, 0);
    
    EXPECT_TRUE(TwoCenterPairComponent() == TwoCenterPairComponent({"", ""}, {scomp, scomp}));
}

TEST_F(TwoCenterPairComponentTest, OperatorNotEqual)
{
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_yz = TensorComponent(0, 1, 1);
    
    EXPECT_TRUE(TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz}) != TwoCenterPairComponent({"LA", "GB"}, {p_x, d_yz}));
    
    EXPECT_TRUE(TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz}) != TwoCenterPairComponent({"GA", "GB"}, {p_x, p_x}));
}

TEST_F(TwoCenterPairComponentTest, OperatorLess)
{
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_yz = TensorComponent(0, 1, 1);
    
    EXPECT_FALSE(TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz}) < TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz}));
    
    EXPECT_TRUE(TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz}) < TwoCenterPairComponent({"LA", "GB"}, {p_x, d_yz}));
    
    EXPECT_TRUE(TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz}) < TwoCenterPairComponent({"GA", "LA"}, {p_x, d_yz}));
    
    EXPECT_FALSE(TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz}) < TwoCenterPairComponent({"GA", "GB"}, {d_yz, p_x}));
    
    EXPECT_TRUE(TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz}) < TwoCenterPairComponent({"GA", "GB"}, {p_x, p_x}));
    
    EXPECT_FALSE(TwoCenterPairComponent({"GA", "GB"}, {p_x, p_x}) < TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz}));
}

TEST_F(TwoCenterPairComponentTest, Names)
{
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_yz = TensorComponent(0, 1, 1);
    
    const auto tpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz});
    
    EXPECT_EQ(tpair.names()[0], "GA");
    
    EXPECT_EQ(tpair.names()[1], "GB");
}

TEST_F(TwoCenterPairComponentTest, Shapes)
{
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_yz = TensorComponent(0, 1, 1);
    
    const auto tpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz});
    
    tpair.shapes();
    
    EXPECT_EQ(tpair.shapes()[0], p_x);
    
    EXPECT_EQ(tpair.shapes()[1], d_yz);
}

TEST_F(TwoCenterPairComponentTest, ToString)
{
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_yz = TensorComponent(0, 1, 1);
    
    const auto tpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz});
    
    EXPECT_EQ(tpair.to_string(), "{GA:(1,0,0);GB:(0,1,1)}");
}

TEST_F(TwoCenterPairComponentTest, Label)
{
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_yz = TensorComponent(0, 1, 1);
    
    const auto tpair = TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz});
    
    EXPECT_EQ(tpair.label(), "x_yz");
}


TEST_F(TwoCenterPairComponentTest, Shift)
{
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto p_z = TensorComponent(0, 0, 1);
    
    const auto d_yz = TensorComponent(0, 1, 1);
    
    const auto t_x_yz = TwoCenterPairComponent({"GA", "GB"}, {p_x, d_yz});
    
    const auto t_0_yz = TwoCenterPairComponent({"GA", "GB"}, {s_0, d_yz});
    
    const auto t_x_y = TwoCenterPairComponent({"GA", "GB"}, {p_x, p_y});
    
    const auto t_x_z = TwoCenterPairComponent({"GA", "GB"}, {p_x, p_z});
    
    EXPECT_EQ(t_x_yz.shift('x', -1, 0), t_0_yz);
    
    EXPECT_EQ(t_x_yz.shift('y', -1, 1), t_x_z);
    
    EXPECT_EQ(t_x_yz.shift('z', -1, 1), t_x_y);
    
    EXPECT_FALSE(t_x_yz.shift('x', -2, 0));
    
    EXPECT_FALSE(t_x_yz.shift('y', -1, 0));
    
    EXPECT_FALSE(t_x_yz.shift('z', -1, 0));
    
    EXPECT_FALSE(t_x_yz.shift('x', -1, 1));
    
    EXPECT_FALSE(t_x_yz.shift('y', -2, 1));
    
    EXPECT_FALSE(t_x_yz.shift('z', -2, 1));
}
