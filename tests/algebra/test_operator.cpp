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

#include "test_operator.hpp"

#include "operator.hpp"
#include "setters.hpp"

TEST_F(OperatorTest, Constructor)
{
    EXPECT_EQ(Operator(), Operator("", Tensor(0), "none", -1));
    
    EXPECT_EQ(Operator("1/r"), Operator("1/r", Tensor(0), "none", -1));
    
    const auto opval = Operator("rxr", Tensor(2), "bra", 2);
    
    for (const auto& tcomp : gset::tensor_components(2))
    {
        EXPECT_EQ(opval, Operator(OperatorComponent("rxr", tcomp, "bra", 2)));
    }
}

TEST_F(OperatorTest, OperatorEqual)
{
    EXPECT_TRUE(Operator("1/r") == Operator("1/r", Tensor(0), "none", -1));
}

TEST_F(OperatorTest, OperatorNotEqual)
{
    EXPECT_TRUE(Operator("1/r") != Operator("r^2", Tensor(0), "none", -1));
    
    EXPECT_TRUE(Operator("1/r") != Operator("1/r", Tensor(1), "none", -1));
    
    EXPECT_TRUE(Operator("1/r") != Operator("1/r", Tensor(0), "bra", -1));
    
    EXPECT_TRUE(Operator("1/r") != Operator("1/r", Tensor(0), "none", 2));
}

TEST_F(OperatorTest, OperatorLess)
{
    EXPECT_FALSE(Operator("1/r") < Operator("1/r", Tensor(0), "none", -1));
    
    EXPECT_TRUE(Operator("1/r") < Operator("1/r", Tensor(0), "none", 0));
    
    EXPECT_FALSE(Operator("1/r") < Operator("1/r", Tensor(0), "bra", -1));
    
    EXPECT_TRUE(Operator("1/r") < Operator("1/r", Tensor(1), "none", -1));
    
    EXPECT_TRUE(Operator("1/r") < Operator("d/dr", Tensor(0), "none", 0));
}

TEST_F(OperatorTest, ToString)
{
    const auto opval = Operator("X", Tensor(3), "bra", 2);
    
    EXPECT_EQ(opval.to_string(), "{X:(3)}[bra:2]");
}

TEST_F(OperatorTest, Label)
{
    const auto opval = Operator("X", Tensor(3), "bra", 2);
    
    EXPECT_EQ(opval.label(), "F");
}

TEST_F(OperatorTest, Components)
{
    const auto opval = Operator("rxr", Tensor(2), "bra", 2);
    
    const auto opcomps = opval.components();
    
    const auto dcomps = gset::tensor_components(2);
    
    EXPECT_EQ(opcomps.size(), 6);
    
    for (size_t i = 0; i < 6; i++)
    {
        EXPECT_EQ(opcomps[i], OperatorComponent("rxr", dcomps[i], "bra", 2));
    }
}

