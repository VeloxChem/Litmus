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

#include "test_four_center_integral.hpp"

#include "four_center_integral.hpp"
#include "setters.hpp"

TEST_F(FourCenterIntegralTest, Constructor)
{
    EXPECT_EQ(FourCenterIntegral(),
              FourCenterIntegral(TwoCenterPair(), TwoCenterPair(), Operator(), 0, {}));
    
    const auto operi = Operator("1/|r-r'|");
    
    EXPECT_EQ(FourCenterIntegral(1, 2, 3, 4, operi),
              FourCenterIntegral(TwoCenterPair("GA", 1, "GB", 2),
                                 TwoCenterPair("GC", 3, "GD", 4),
                                 operi, 0, {}));
    
    EXPECT_EQ(FourCenterIntegral(1, 2, 3, 4, operi, 1),
              FourCenterIntegral(TwoCenterPair("GA", 1, "GB", 2),
                                 TwoCenterPair("GC", 3, "GD", 4),
                                 operi, 1, {}));
    
    const auto opddr = Operator("d/dr", Tensor(1), "bra", 1);
    
    const auto opddc = Operator("d/dC", Tensor(1), "ket", 0);
    
    EXPECT_EQ(FourCenterIntegral(1, 2, 3, 4, operi, 1, {opddr, opddc}),
              FourCenterIntegral(TwoCenterPair("GA", 1, "GB", 2),
                                 TwoCenterPair("GC", 3, "GD", 4),
                                 operi, 1, {opddr, opddc}));
}

TEST_F(FourCenterIntegralTest, OperatorEqual)
{
    const auto operi = Operator("1/|r-r'|");
    
    const auto opddr = Operator("d/dr", Tensor(1), "bra", 1);
    
    const auto opddc = Operator("d/dC", Tensor(1), "ket", 0);
    
    const auto lhsint = FourCenterIntegral(1, 2, 3, 4, operi, 1, {opddr, opddc});
    
    const auto rhsint = FourCenterIntegral(TwoCenterPair("GA", 1, "GB", 2),
                                           TwoCenterPair("GC", 3, "GD", 4),
                                           operi, 1, {opddr, opddc});
    
    EXPECT_TRUE(lhsint == rhsint);
}

TEST_F(FourCenterIntegralTest, OperatorNotEqual)
{
    const auto operi = Operator("1/|r-r'|");
    
    const auto opddr = Operator("d/dr", Tensor(1), "bra", 1);
    
    const auto opddc = Operator("d/dC", Tensor(1), "ket", 0);
    
    const auto lhsint = FourCenterIntegral(1, 2, 3, 4, operi, 1, {opddr, opddc});
    
    EXPECT_TRUE(lhsint != FourCenterIntegral(0, 2, 3, 4, operi, 1, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint != FourCenterIntegral(1, 4, 3, 4, operi, 1, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint != FourCenterIntegral(1, 2, 2, 4, operi, 1, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint != FourCenterIntegral(1, 2, 3, 2, operi, 1, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint != FourCenterIntegral(1, 2, 3, 4, opddr, 1, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint != FourCenterIntegral(1, 2, 3, 4, operi, 0, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint != FourCenterIntegral(1, 2, 3, 4, operi, 1, {opddr, opddr}));
}

TEST_F(FourCenterIntegralTest, OperatorLess)
{
    const auto operi = Operator("1/|r-r'|");
    
    const auto opddr = Operator("d/dr", Tensor(1), "bra", 1);
    
    const auto opddc = Operator("d/dC", Tensor(1), "ket", 0);
    
    const auto lhsint = FourCenterIntegral(1, 2, 3, 4, operi, 1, {opddr, opddc});
    
    EXPECT_FALSE(lhsint < lhsint);
    
    EXPECT_TRUE(lhsint < FourCenterIntegral(2, 2, 3, 4, operi, 1, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint < FourCenterIntegral(1, 3, 3, 4, operi, 1, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint < FourCenterIntegral(1, 2, 4, 4, operi, 1, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint < FourCenterIntegral(1, 2, 3, 5, operi, 1, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint < FourCenterIntegral(1, 2, 3, 4, opddr, 1, {opddr, opddc}));
    
    EXPECT_FALSE(lhsint < FourCenterIntegral(1, 2, 3, 4, operi, 0, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint < FourCenterIntegral(1, 2, 3, 4, operi, 1, {opddr, opddr}));
}

TEST_F(FourCenterIntegralTest, ToString)
{
    const auto operi = Operator("1/|r-r'|");
        
    auto t4cint = FourCenterIntegral(1, 2, 3, 4, operi);
    
    EXPECT_EQ(t4cint.to_string(),
              "{GA:(1);GB:(2)}{1/|r-r'|:(0)}[none:-1]{GC:(3);GD:(4)}^(0)");
    
    const auto opddr = Operator("d/dr", Tensor(1), "bra", 1);
    
    const auto opddc = Operator("d/dC", Tensor(1), "ket", 0);
    
    t4cint = FourCenterIntegral(1, 2, 3, 4, operi, 1, {opddr, opddc});
    
    EXPECT_EQ(t4cint.to_string(),
              "[{d/dr:(1)}[bra:1];{d/dC:(1)}[ket:0];]{GA:(1);GB:(2)}{1/|r-r'|:(0)}[none:-1]{GC:(3);GD:(4)}^(1)");
}

TEST_F(FourCenterIntegralTest, Label)
{
    const auto operi = Operator("1/|r-r'|");
        
    auto t4cint = FourCenterIntegral(1, 2, 3, 4, operi);
    
    EXPECT_EQ(t4cint.label(), "PDFG");
    
    EXPECT_EQ(t4cint.label(true), "PDFG_0");
    
    t4cint = FourCenterIntegral(1, 2, 3, 4, operi, 2);
    
    EXPECT_EQ(t4cint.label(), "PDFG");
    
    EXPECT_EQ(t4cint.label(true), "PDFG_2");
    
    const auto opddr = Operator("d/dr", Tensor(1), "bra", 1);
    
    const auto opddc = Operator("d/dC", Tensor(1), "ket", 0);
    
    t4cint = FourCenterIntegral(1, 2, 3, 4, operi, 1, {opddr, opddc});
    
    EXPECT_EQ(t4cint.label(), "PP_PDFG");
    
    EXPECT_EQ(t4cint.label(true), "PP_PDFG_1");
}

