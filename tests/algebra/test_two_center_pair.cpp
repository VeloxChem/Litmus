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

#include "test_two_center_pair.hpp"

#include "two_center_pair.hpp"
#include "setters.hpp"

TEST_F(TwoCenterPairTest, Constructor)
{
    auto lhspair = TwoCenterPair();
    
    auto rhspair = TwoCenterPair({"", ""}, {Tensor(0), Tensor(0)});
    
    EXPECT_EQ(lhspair, rhspair);
    
    lhspair = TwoCenterPair("GA", 2, "GB", 1);
    
    rhspair = TwoCenterPair({"GA", "GB"}, {Tensor(2), Tensor(1)});
    
    EXPECT_EQ(lhspair, rhspair);
    
    for (const auto& t2pcomp : gset::two_center_pair_components("GA", 2, "GB", 1))
    {
        EXPECT_EQ(lhspair, TwoCenterPair(t2pcomp));
    }
}

TEST_F(TwoCenterPairTest, OperatorEqual)
{
    const auto lhspair = TwoCenterPair("GA", 2, "GB", 1);
    
    const auto rhspair = TwoCenterPair({"GA", "GB"}, {Tensor(2), Tensor(1)});
    
    EXPECT_TRUE( rhspair == lhspair);
}

TEST_F(TwoCenterPairTest, OperatorNotEqual)
{
    EXPECT_TRUE(TwoCenterPair("GA", 2, "GB", 1) != TwoCenterPair("LA", 2, "GB", 1));
    
    EXPECT_TRUE(TwoCenterPair("GA", 2, "GB", 1) != TwoCenterPair("GA", 0, "GB", 1));
   
    EXPECT_TRUE(TwoCenterPair("GA", 2, "GB", 1) != TwoCenterPair("GA", 2, "GA", 1));
    
    EXPECT_TRUE(TwoCenterPair("GA", 2, "GB", 1) != TwoCenterPair("GA", 2, "GB", 0));
}

TEST_F(TwoCenterPairTest, OperatorLess)
{
    EXPECT_FALSE(TwoCenterPair("GA", 2, "GB", 1) < TwoCenterPair("GA", 2, "GB", 1));
    
    EXPECT_TRUE(TwoCenterPair("GA", 2, "GB", 1) < TwoCenterPair("LA", 2, "GB", 1));
    
    EXPECT_TRUE(TwoCenterPair("GA", 2, "GB", 1) < TwoCenterPair("GA", 3, "GB", 1));
    
    EXPECT_TRUE(TwoCenterPair("GA", 2, "GA", 1) < TwoCenterPair("GA", 2, "GB", 1));
   
    EXPECT_FALSE(TwoCenterPair("GA", 2, "GB", 1) < TwoCenterPair("GA", 2, "GA", 1));
    
    EXPECT_FALSE(TwoCenterPair("GA", 2, "GB", 1) < TwoCenterPair("GA", 2, "GB", 0));
}

TEST_F(TwoCenterPairTest, ToString)
{
    const auto tpair = TwoCenterPair("GA", 2, "GB", 1);
    
    EXPECT_EQ(tpair.to_string(), "{GA:(2);GB:(1)}");
}

TEST_F(TwoCenterPairTest, Label)
{
    const auto tpair = TwoCenterPair("GA", 2, "GB", 1);
    
    EXPECT_EQ(tpair.label(), "DP");
}

TEST_F(TwoCenterPairTest, Components)
{
    const auto tpair = TwoCenterPair("GA", 2, "GB", 1);

    const auto t2pcomps = tpair.components();
    
    const auto dpcomps = gset::two_center_pair_components("GA", 2, "GB", 1);
    
    EXPECT_EQ(t2pcomps.size(), 18);
    
    for (size_t i = 0; i < 18; i++)
    {
        EXPECT_EQ(dpcomps[i], t2pcomps[i]);
    }
}
