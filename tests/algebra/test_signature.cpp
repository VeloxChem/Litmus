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

#include "test_signature.hpp"

#include "signature.hpp"

using IntSign = Signature<int>;

TEST_F(SignatureTest, Constructor)
{
    const auto tsign = IntSign();
    
    const auto rsign = IntSign({}, {}, {});
    
    EXPECT_EQ(tsign, rsign);
}

TEST_F(SignatureTest, OperatorEqual)
{
    const auto pbx = Factor("(P-B)", "pb");
    
    const auto wpy = Factor("(W-P)", "wp");
    
    const auto lhsrt = IntSign({1, 3}, {2, 4}, {pbx, wpy});
    
    const auto rhsrt = IntSign({3, 1}, {4, 2}, {wpy, pbx});
    
    EXPECT_TRUE(lhsrt == rhsrt);
}

TEST_F(SignatureTest, OperatorNotEqual)
{
    const auto pbx = Factor("(P-B)", "pb");
    
    const auto wpy = Factor("(W-P)", "wp");
    
    const auto lhsrt = IntSign({1, 3}, {2, 4}, {pbx, wpy});
    
    auto rhsrt = IntSign({3, 1}, {4, 2}, {wpy, pbx});
    
    EXPECT_FALSE(lhsrt != rhsrt);
    
    rhsrt = IntSign({3, 2}, {4, 2}, {wpy, pbx});
    
    EXPECT_TRUE(lhsrt != rhsrt);
    
    rhsrt = IntSign({3, 1}, {4, 2, 5}, {wpy, pbx});
    
    EXPECT_TRUE(lhsrt != rhsrt);
    
    rhsrt = IntSign({3, 1}, {4, 2}, {wpy,});
    
    EXPECT_TRUE(lhsrt != rhsrt);
}

TEST_F(SignatureTest, OperatorLess)
{
    const auto pbx = Factor("(P-B)", "pb");
    
    const auto wpy = Factor("(W-P)", "wp");
    
    const auto lhsrt = IntSign({1, 3}, {2, 4}, {pbx, wpy});
    
    auto rhsrt = IntSign({3, 1}, {4, 2}, {wpy, pbx});
    
    EXPECT_FALSE(lhsrt < rhsrt);
    
    rhsrt = IntSign({3, 2}, {4, 2}, {wpy, pbx});
    
    EXPECT_TRUE(lhsrt < rhsrt);
    
    rhsrt = IntSign({3, 1}, {4, 2, 5}, {wpy, pbx});
    
    EXPECT_TRUE(lhsrt < rhsrt);
    
    rhsrt = IntSign({3, 1}, {4, 2}, {wpy,});
    
    EXPECT_TRUE(lhsrt < rhsrt);
}

TEST_F(SignatureTest, Merge)
{
    const auto pbx = Factor("(P-B)", "pb");
    
    const auto wpy = Factor("(W-P)", "wp");
    
    auto lhsrt = IntSign({1, 3}, {2, 4}, {pbx,});
    
    auto rhsrt = IntSign({3, 2}, {1, 8, 2}, {pbx, wpy,});
    
    lhsrt.merge(rhsrt);
    
    rhsrt = IntSign({1, 2, 3}, {1, 2, 4, 8}, {pbx, wpy,});
    
    EXPECT_EQ(lhsrt, rhsrt);
}

TEST_F(SignatureTest, AddParam)
{
    const auto pbx = Factor("(P-B)", "pb");
    
    const auto wpy = Factor("(W-P)", "wp");
    
    auto lhsrt = IntSign({}, {}, {pbx, wpy});
    
    auto rhsrt = IntSign({}, {}, {wpy, pbx});
    
    EXPECT_EQ(lhsrt, rhsrt);
    
    lhsrt.add(1, "inp");
    
    rhsrt = IntSign({}, {1,}, {wpy, pbx});
    
    EXPECT_EQ(lhsrt, rhsrt);
    
    lhsrt.add(4, "inp");
    
    rhsrt = IntSign({}, {1, 4}, {wpy, pbx});
    
    EXPECT_EQ(lhsrt, rhsrt);
    
    lhsrt.add(3, "out");
    
    rhsrt = IntSign({3,}, {1, 4}, {wpy, pbx});
    
    EXPECT_EQ(lhsrt, rhsrt);
    
    lhsrt.add(5, "out");
    
    rhsrt = IntSign({3, 5}, {1, 4}, {wpy, pbx});
    
    EXPECT_EQ(lhsrt, rhsrt);
}

TEST_F(SignatureTest, AddFactor)
{
    const auto pbx = Factor("(P-B)", "pb");
    
    const auto wpy = Factor("(W-P)", "wp");
    
    auto lhsrt = IntSign({1, 3}, {2, 4}, {});
    
    auto rhsrt = IntSign({3, 1}, {4, 2}, {});
    
    EXPECT_EQ(lhsrt, rhsrt);
    
    lhsrt.add(pbx);
    
    rhsrt = IntSign({3, 1}, {4, 2}, {pbx});
    
    EXPECT_EQ(lhsrt, rhsrt);
    
    lhsrt.add(pbx);
    
    rhsrt = IntSign({3, 1}, {4, 2}, {pbx});
    
    EXPECT_EQ(lhsrt, rhsrt);
    
    lhsrt.add(wpy);
    
    rhsrt = IntSign({3, 1}, {4, 2}, {pbx, wpy});
    
    EXPECT_EQ(lhsrt, rhsrt);
}

TEST_F(SignatureTest, NParams)
{
    const auto pbx = Factor("(P-B)", "pb");
    
    const auto wpy = Factor("(W-P)", "wp");
    
    auto rhsrt = IntSign({1, 3, 5}, {2, 4, 7, 8}, {pbx, wpy});
    
    EXPECT_EQ(rhsrt.nparams("out"), 3);
    
    EXPECT_EQ(rhsrt.nparams("inp"), 4);
}

TEST_F(SignatureTest, NFactors)
{
    const auto pbx = Factor("(P-B)", "pb");
    
    const auto wpy = Factor("(W-P)", "wp");
    
    auto rhsrt = IntSign({1, 3, 5}, {2, 4, 7, 8}, {pbx, wpy});
    
    EXPECT_EQ(rhsrt.nfactors(), 2);
}

TEST_F(SignatureTest, Params)
{
    const auto pbx = Factor("(P-B)", "pb");
    
    const auto wpy = Factor("(W-P)", "wp");
    
    auto rhsrt = IntSign({1, 3, 5}, {2, 4, 7, 8}, {pbx, wpy});
    
    EXPECT_EQ(rhsrt.params("out"), std::set<int>({1, 3, 5}));
    
    EXPECT_EQ(rhsrt.params("inp"), std::set<int>({2, 4, 7, 8}));
}

TEST_F(SignatureTest, Factors)
{
    const auto pbx = Factor("(P-B)", "pb");
    
    const auto wpy = Factor("(W-P)", "wp");
    
    auto rhsrt = IntSign({1, 3, 5}, {2, 4, 7, 8}, {pbx, wpy});
    
    EXPECT_EQ(rhsrt.factors(), std::set<Factor>({pbx, wpy}));
}

