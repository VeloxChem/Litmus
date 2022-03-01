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

#include "test_operator_component.hpp"

#include "operator_component.hpp"

#include "setters.hpp"

TEST_F(OperatorComponentTest, Constructor)
{
    const auto scomp = TensorComponent(0, 0, 0);
    
    EXPECT_EQ(OperatorComponent(), OperatorComponent("", scomp, "none", -1));
    
    EXPECT_EQ(OperatorComponent("1/r"), OperatorComponent("1/r", scomp, "none", -1));
}

TEST_F(OperatorComponentTest, OperatorBrackets)
{
    const auto opcomp = OperatorComponent("X", TensorComponent(2, 1, 5), "ket", 2);
    
    EXPECT_EQ(opcomp['x'], 2);
    
    EXPECT_EQ(opcomp['y'], 1);
    
    EXPECT_EQ(opcomp['z'], 5);
    
    EXPECT_EQ(opcomp['g'], -1);
}

TEST_F(OperatorComponentTest, OperatorEqual)
{
    const auto scomp = TensorComponent(0, 0, 0);
    
    EXPECT_TRUE(OperatorComponent() == OperatorComponent("", scomp, "none", -1));
    
    EXPECT_TRUE(OperatorComponent("1/r") == OperatorComponent("1/r", scomp, "none", -1));
}

TEST_F(OperatorComponentTest, OperatorNotEqual)
{
    const auto scomp = TensorComponent(0, 0, 0);
    
    const auto pcomp = TensorComponent(0, 0, 1);
    
    EXPECT_TRUE(OperatorComponent("1/r") != OperatorComponent("r^2", scomp, "none", -1));
    
    EXPECT_TRUE(OperatorComponent("1/r") != OperatorComponent("1/r", pcomp, "none", -1));
    
    EXPECT_TRUE(OperatorComponent("1/r") != OperatorComponent("1/r", scomp, "bra", -1));
    
    EXPECT_TRUE(OperatorComponent("1/r") != OperatorComponent("1/r", scomp, "none", 2));
}

TEST_F(OperatorComponentTest, OperatorLess)
{
    const auto scomp = TensorComponent(0, 0, 0);
    
    const auto pcomp = TensorComponent(0, 0, 1);
    
    EXPECT_FALSE(OperatorComponent("1/r") < OperatorComponent("1/r", scomp, "none", -1));
    
    EXPECT_TRUE(OperatorComponent("1/r") < OperatorComponent("1/r", scomp, "none", 0));
    
    EXPECT_FALSE(OperatorComponent("1/r") < OperatorComponent("1/r", scomp, "bra", -1));
    
    EXPECT_TRUE(OperatorComponent("1/r") < OperatorComponent("1/r", pcomp, "none", -1));
    
    EXPECT_TRUE(OperatorComponent("1/r") < OperatorComponent("d/dr", scomp, "none", 0));
}

TEST_F(OperatorComponentTest, Name)
{
    const auto pcomp = TensorComponent(0, 0, 1);
    
    const auto opval = OperatorComponent("d/dr", pcomp, "ket", 1);
    
    EXPECT_EQ(opval.name(), "d/dr");
}

TEST_F(OperatorComponentTest, Shape)
{
    const auto pcomp = TensorComponent(0, 0, 1);
    
    const auto opval = OperatorComponent("d/dr", pcomp, "ket", 1);
    
    EXPECT_EQ(opval.shape(), pcomp);
}

TEST_F(OperatorComponentTest, Target)
{
    const auto pcomp = TensorComponent(0, 0, 1);
    
    const auto opval = OperatorComponent("d/dr", pcomp, "ket", 1);
    
    EXPECT_EQ(opval.target(), "ket");
}

TEST_F(OperatorComponentTest, Center)
{
    const auto pcomp = TensorComponent(0, 0, 1);
    
    const auto opval = OperatorComponent("d/dr", pcomp, "ket", 1);
    
    EXPECT_EQ(opval.center(), 1);
}

TEST_F(OperatorComponentTest, ToString)
{
    const auto pcomp = TensorComponent(0, 0, 1);
    
    const auto opval = OperatorComponent("d/dr", pcomp, "ket", 2);
    
    EXPECT_EQ(opval.to_string(), "{d/dr:(0,0,1)}[ket:2]");
}

TEST_F(OperatorComponentTest, Label)
{
    const auto gcomp = TensorComponent(1, 2, 1);
    
    const auto opval = OperatorComponent("d/dr", gcomp, "ket", 2);
    
    EXPECT_EQ(opval.label(), "xyyz");
}

TEST_F(OperatorComponentTest, Shift)
{
    const auto rxxy = OperatorComponent("r^n", TensorComponent(2, 1, 0), "ket", 2);
    
    const auto rxyz = OperatorComponent("r^n", TensorComponent(1, 1, 1), "ket", 2);
    
    const auto rxx = OperatorComponent("r^n", TensorComponent(2, 0, 0), "ket", 2);
    
    const auto rxy = OperatorComponent("r^n", TensorComponent(1, 1, 0), "ket", 2);
    
    const auto rxz = OperatorComponent("r^n", TensorComponent(1, 0, 1), "ket", 2);
    
    const auto ryz = OperatorComponent("r^n", TensorComponent(0, 1, 1), "ket", 2);
    
    const auto r0 = OperatorComponent("r^n", TensorComponent(0, 0, 0), "ket", 2);
    
    // without no scalar
    
    EXPECT_EQ(rxxy.shift('x', -1), rxy);
    
    EXPECT_EQ(rxxy.shift('y', -1), rxx);
    
    EXPECT_EQ(rxyz.shift('x', -1), ryz);
    
    EXPECT_EQ(rxyz.shift('y', -1), rxz);
    
    EXPECT_EQ(rxyz.shift('z', -1), rxy);
    
    EXPECT_EQ(rxx.shift('x', -2), r0);
    
    // with no scalar
    
    EXPECT_EQ(rxxy.shift('x', -1, true), rxy);
    
    EXPECT_EQ(rxxy.shift('y', -1, true), rxx);
    
    EXPECT_EQ(rxyz.shift('x', -1, true), ryz);
    
    EXPECT_EQ(rxyz.shift('y', -1, true), rxz);
    
    EXPECT_EQ(rxyz.shift('z', -1, true), rxy);
    
    EXPECT_FALSE(rxx.shift('x', -2, true));
    
    // non-valid shift
    
    EXPECT_FALSE(rxxy.shift('x', -3));
    
    EXPECT_FALSE(rxxy.shift('y', -2));
    
    EXPECT_FALSE(rxxy.shift('z', -1));
    
    EXPECT_FALSE(rxyz.shift('x', -2));
    
    EXPECT_FALSE(rxyz.shift('y', -2));
    
    EXPECT_FALSE(rxyz.shift('z', -2));
}
