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

#include "test_factor.hpp"

#include "factor.hpp"

TEST_F(FactorTest, Constructor)
{
    const auto scomp = TensorComponent(0, 0, 0);
    
    EXPECT_EQ(Factor(), Factor("", "", scomp));
    
    EXPECT_EQ(Factor("(P-B)", "pb"), Factor("(P-B)", "pb", scomp));
}

TEST_F(FactorTest, OperatorEqual)
{
    const auto scomp = TensorComponent(0, 0, 0);
    
    EXPECT_TRUE(Factor("(P-B)", "pb") == Factor("(P-B)", "pb", scomp));
}

TEST_F(FactorTest, OperatorNotEqual)
{
    const auto scomp = TensorComponent(0, 0, 0);
    
    const auto pcomp = TensorComponent(0, 0, 1);
    
    EXPECT_TRUE(Factor("(P-B)", "pb") != Factor("(P+B)", "pb", scomp));
    
    EXPECT_TRUE(Factor("(P-B)", "pb") != Factor("(P-B)", "rpb", scomp));
    
    EXPECT_TRUE(Factor("(P-B)", "pb") != Factor("(P-B)", "pb", pcomp));
}

TEST_F(FactorTest, OperatorLess)
{
    const auto scomp = TensorComponent(0, 0, 0);
    
    const auto pcomp = TensorComponent(0, 0, 1);
    
    EXPECT_FALSE(Factor("(P-B)", "pb") < Factor("(P+B)", "pb", scomp));
    
    EXPECT_FALSE(Factor("(P-B)", "pb") < Factor("(P-B)", "rpb", scomp));
    
    EXPECT_TRUE(Factor("(P-B)", "pb") < Factor("(P-B)", "pb", pcomp));
}

TEST_F(FactorTest, ToString)
{
    const auto pcomp = TensorComponent(0, 0, 1);
    
    const auto fact = Factor("(P-B)", "pb", pcomp);
    
    EXPECT_EQ(fact.to_string(), "{(P-B)(pb):(0,0,1)}");
}

TEST_F(FactorTest, Label)
{
    const auto scomp = TensorComponent(0, 0, 0);
    
    const auto pcomp = TensorComponent(0, 1, 0);
    
    auto fact = Factor("1/eta", "fz", scomp);
    
    EXPECT_EQ(fact.label(), "fz");
    
    fact = Factor("(P-B)", "pb", pcomp);
    
    EXPECT_EQ(fact.label(), "pb_y");
}
