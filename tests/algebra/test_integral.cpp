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

#include "test_integral.hpp"

#include "two_center_pair.hpp"
#include "integral.hpp"
#include "setters.hpp"

using T2CPairComp = TwoCenterPairComponent;

using T4CIntegralComp = IntegralComponent<T2CPairComp, T2CPairComp>;

using T2CPair = TwoCenterPair;

using T4CIntegral = Integral<T2CPair, T2CPair>;

TEST_F(IntegralTest, Constructor)
{
    EXPECT_EQ(T4CIntegral(), T4CIntegral(T2CPair(), T2CPair(), Operator(), 0, {}));
    
    const auto operi = Operator("1/|r-r'|");
    
    const auto bpair = T2CPair("GA", 1, "GB", 2);
    
    const auto kpair = T2CPair("GC", 3, "GD", 4);
                        
    EXPECT_EQ(T4CIntegral(bpair, kpair, operi),
              T4CIntegral(bpair, kpair, operi, 0, {}));
    
    EXPECT_EQ(T4CIntegral(bpair, kpair, operi, 1),
              T4CIntegral(bpair, kpair, operi, 1, {}));
    
    const auto opddr = Operator("d/dr", Tensor(1), "bra", 1);
    
    const auto opddc = Operator("d/dC", Tensor(1), "ket", 0);
    
    const auto lhsint = T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc});
    
    const auto opceri = OperatorComponent("1/|r-r'|");
    
    const auto opcddr = OperatorComponent("d/dr", TensorComponent(1, 0, 0), "bra", 1);
    
    const auto opcddc = OperatorComponent("d/dC", TensorComponent(0, 1, 0), "ket", 0);
    
    const auto bcomp = T2CPairComp({"GA", "GB"}, {TensorComponent(1, 0, 0), TensorComponent(1, 0, 1)});
    
    const auto kcomp = T2CPairComp({"GC", "GD"}, {TensorComponent(1, 2, 0), TensorComponent(1, 2, 1)});
    
    const auto rhsint = T4CIntegral(T4CIntegralComp(bcomp, kcomp, opceri, 1, {opcddr, opcddc}));
    
    EXPECT_EQ(lhsint, rhsint);
}

TEST_F(IntegralTest, OperatorEqual)
{
    const auto operi = Operator("1/|r-r'|");
    
    const auto opddr = Operator("d/dr", Tensor(1), "bra", 1);
    
    const auto opddc = Operator("d/dC", Tensor(1), "ket", 0);
    
    const auto bpair = T2CPair("GA", 1, "GB", 2);
    
    const auto kpair = T2CPair("GC", 3, "GD", 4);
    
    const auto lhsint = T4CIntegral(bpair, kpair, operi);
    
    const auto rhsint = T4CIntegral(bpair, kpair, operi, 0, {});
    
    EXPECT_TRUE(lhsint == rhsint);
}

TEST_F(IntegralTest, OperatorNotEqual)
{
    const auto operi = Operator("1/|r-r'|");
    
    const auto opddr = Operator("d/dr", Tensor(1), "bra", 1);
    
    const auto opddc = Operator("d/dC", Tensor(1), "ket", 0);
    
    auto bpair = T2CPair("GA", 1, "GB", 2);
    
    auto kpair = T2CPair("GC", 3, "GD", 4);
    
    const auto lhsint = T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc});
    
    bpair = T2CPair("GA", 1, "GB", 4);
    
    EXPECT_TRUE(lhsint != T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    bpair = T2CPair("GA", 1, "GB", 3);
    
    EXPECT_TRUE(lhsint != T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    bpair = T2CPair("LA", 1, "GB", 2);
    
    EXPECT_TRUE(lhsint != T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    bpair = T2CPair("GA", 1, "LB", 2);
    
    EXPECT_TRUE(lhsint != T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    bpair = T2CPair("GA", 1, "GB", 2);
    
    kpair = T2CPair("GC", 0, "GD", 4);
    
    EXPECT_TRUE(lhsint != T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    kpair = T2CPair("GC", 3, "GD", 3);
    
    EXPECT_TRUE(lhsint != T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    kpair = T2CPair("LC", 3, "GD", 4);
    
    EXPECT_TRUE(lhsint != T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    kpair = T2CPair("GC", 3, "LD", 4);
    
    EXPECT_TRUE(lhsint != T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    kpair = T2CPair("GC", 3, "GD", 4);
    
    EXPECT_TRUE(lhsint != T4CIntegral(bpair, kpair, operi, 0, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint != T4CIntegral(bpair, kpair, operi, 1, {opddr, opddr}));
}

TEST_F(IntegralTest, OperatorLess)
{
    const auto operi = Operator("1/|r-r'|");
    
    const auto opddr = Operator("d/dr", Tensor(1), "bra", 1);
    
    const auto opddc = Operator("d/dC", Tensor(1), "ket", 0);
    
    auto bpair = T2CPair("GA", 1, "GB", 2);
    
    auto kpair = T2CPair("GC", 3, "GD", 4);
    
    const auto lhsint = T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc});
    
    EXPECT_FALSE(lhsint < lhsint);
    
    bpair = T2CPair("GA", 1, "GB", 4);
    
    EXPECT_TRUE(lhsint < T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    bpair = T2CPair("GA", 1, "GB", 3);
    
    EXPECT_TRUE(lhsint < T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    bpair = T2CPair("LA", 1, "GB", 2);
    
    EXPECT_TRUE(lhsint < T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    bpair = T2CPair("GA", 1, "LB", 2);
    
    EXPECT_TRUE(lhsint < T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    bpair = T2CPair("GA", 1, "GB", 2);
    
    kpair = T2CPair("GC", 0, "GD", 4);
    
    EXPECT_FALSE(lhsint < T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    kpair = T2CPair("GC", 3, "GD", 3);
    
    EXPECT_FALSE(lhsint < T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    kpair = T2CPair("LC", 3, "GD", 4);
    
    EXPECT_TRUE(lhsint < T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    kpair = T2CPair("GC", 3, "LD", 4);
    
    EXPECT_TRUE(lhsint < T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc}));
    
    kpair = T2CPair("GC", 3, "GD", 4);
    
    EXPECT_TRUE(lhsint < T4CIntegral(bpair, kpair, operi, 3, {opddr, opddc}));
    
    EXPECT_TRUE(lhsint < T4CIntegral(bpair, kpair, operi, 1, {opddr, opddr}));
}

TEST_F(IntegralTest, Label)
{
    const auto operi = Operator("1/|r-r'|");
    
    auto bpair = T2CPair("GA", 1, "GB", 2);
    
    auto kpair = T2CPair("GC", 3, "GD", 4);
        
    auto t4cint = T4CIntegral(bpair, kpair, operi);
    
    EXPECT_EQ(t4cint.label(), "PDFG");
    
    EXPECT_EQ(t4cint.label(true), "PDFG_0");
    
    t4cint = T4CIntegral(bpair, kpair, operi, 2);
    
    EXPECT_EQ(t4cint.label(), "PDFG");
    
    EXPECT_EQ(t4cint.label(true), "PDFG_2");
    
    const auto opddr = Operator("d/dr", Tensor(1), "bra", 1);
    
    const auto opddc = Operator("d/dC", Tensor(1), "ket", 0);
    
    t4cint = T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc});
    
    EXPECT_EQ(t4cint.label(), "PDFG");
    
    EXPECT_EQ(t4cint.label(true), "PDFG_1");
}

TEST_F(IntegralTest, Components)
{
    const auto operi = Operator("1/|r-r'|");
    
    const auto opddr = Operator("d/dr", Tensor(1), "bra", 1);
    
    const auto opddc = Operator("d/dC", Tensor(1), "ket", 0);
    
    const auto bpair = T2CPair("GA", 1, "GB", 2);
    
    const auto kpair = T2CPair("GC", 0, "GD", 3);
    
    const auto t4cint = T4CIntegral(bpair, kpair, operi, 1, {opddr, opddc});
    
    const auto vt4comps = t4cint.components<T2CPairComp, T2CPairComp>();
    
    EXPECT_EQ(vt4comps.size(), 1620);
    
    int idx = 0;
    
    for (const auto& drcomp : opddr.components())
    {
        for (const auto& dccomp : opddc.components())
        {
            for (const auto& opcomp : operi.components())
            {
                for (const auto& bcomp : bpair.components())
                {
                    for (const auto& kcomp : kpair.components())
                    {
                        const auto rhsint = T4CIntegralComp(bcomp, kcomp, opcomp, 1, {drcomp, dccomp});
                        
                        EXPECT_EQ(vt4comps[idx], rhsint);
                        
                        idx++;
                    }
                }
            }
        }
    }
}
    
TEST_F(IntegralTest, DiagComponents)
{
    const auto operi = Operator("1/|r-r'|");
    
    const auto opddr = Operator("d/dr", Tensor(1), "bra", 1);
    
    const auto opddc = Operator("d/dC", Tensor(1), "ket", 0);
    
    const auto bpair = T2CPair("GA", 1, "GB", 2);
    
    const auto t4cint = T4CIntegral(bpair, bpair, operi, 1, {opddr, opddc});
    
    const auto vt4comps = t4cint.diag_components<T2CPairComp, T2CPairComp>();
    
    EXPECT_EQ(vt4comps.size(), 162);
    
    int idx = 0;
    
    for (const auto& drcomp : opddr.components())
    {
        for (const auto& dccomp : opddc.components())
        {
            for (const auto& opcomp : operi.components())
            {
                for (const auto& bcomp : bpair.components())
                {
                    const auto rhsint = T4CIntegralComp(bcomp, bcomp, opcomp, 1, {drcomp, dccomp});
                        
                    EXPECT_EQ(vt4comps[idx], rhsint);
                        
                    idx++;
                }
            }
        }
    }
}
