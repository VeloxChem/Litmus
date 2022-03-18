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

#include "test_graph.hpp"

#include "graph.hpp"

TEST_F(GraphTest, Constructor)
{
    EXPECT_EQ(Graph<int>(), Graph<int>({}, {}));
    
    EXPECT_EQ(Graph<int>(1), Graph<int>({1,}, {{}, }));
    
    const auto lhsg = Graph<int>({1, }, {{2, 1}, });
    
    const auto rhsg = Graph<int>({1, }, {{1, 2}, });
    
    EXPECT_EQ(lhsg, rhsg);
}

TEST_F(GraphTest, OperatorEqual)
{
    const auto lhsg = Graph<int>({1, 7, }, {{2, 1}, {7, 3}, });
 
    const auto rhsg = Graph<int>({1, 7, }, {{1, 2}, {3, 7}, });
    
    EXPECT_TRUE(lhsg == rhsg);
}

TEST_F(GraphTest, OperatorNotEqual)
{
    const auto lhsg = Graph<int>({1, 7, }, {{2, 1}, {7, 3}, });
 
    auto rhsg = Graph<int>({2, 7, }, {{1, 2}, {3, 7}, });
    
    EXPECT_TRUE(lhsg != rhsg);
    
    rhsg = Graph<int>({2, }, {{1, 2}});
    
    EXPECT_TRUE(lhsg != rhsg);
    
    rhsg = Graph<int>({2, 7, }, {{1, 2}, {1, 7}, });
    
    EXPECT_TRUE(lhsg != rhsg);
}

TEST_F(GraphTest, OperatorLess)
{
    const auto lhsg = Graph<int>({1, 7, }, {{2, 1}, {7, 3}, });
    
    EXPECT_FALSE(lhsg < lhsg);
 
    auto rhsg = Graph<int>({2, 7, }, {{1, 2}, {3, 7}, });
    
    EXPECT_TRUE(lhsg < rhsg);
    
    rhsg = Graph<int>({1, 7, }, {{1, 2}, {4, 7}, });
    
    EXPECT_TRUE(lhsg < rhsg);
}

TEST_F(GraphTest, Add)
{
    auto tg = Graph<std::string>("A");
    
    auto rg = Graph<std::string>({"A", }, {{}, });
    
    EXPECT_EQ(tg, rg);
    
    tg.add("B", 0);
    
    rg = Graph<std::string>({"A", "B", },
                            {{1}, {},});
    
    EXPECT_EQ(tg, rg);
    
    tg.add("C", 0);
    
    rg = Graph<std::string>({"A",     "B", "C",},
                            {{1, 2,}, {},  {},});
    
    EXPECT_EQ(tg, rg);
    
    tg.add("D", 1);
    
    rg = Graph<std::string>({"A",     "B", "C", "D",},
                            {{1, 2,}, {3}, {},  {}, });
    
    EXPECT_EQ(tg, rg);
    
    tg.add("E", 1);
    
    rg = Graph<std::string>({"A",     "B",    "C", "D", "E",},
                            {{1, 2,}, {3, 4},  {},  {},  {},});
    
    EXPECT_EQ(tg, rg);
    
    tg.add("E", 2);
    
    rg = Graph<std::string>({"A",     "B",    "C", "D", "E",},
                            {{1, 2,}, {3, 4}, {4},  {},  {},});
    
    EXPECT_EQ(tg, rg);
}

TEST_F(GraphTest, AddWithoutIndex)
{
    auto tg = Graph<std::string>("A");
    
    auto rg = Graph<std::string>({"A", }, {{}, });
    
    EXPECT_EQ(tg, rg);
    
    tg.add("B", "A");
    
    rg = Graph<std::string>({"A", "B", },
                            {{1}, {},});
    
    EXPECT_EQ(tg, rg);
    
    tg.add("C", "A");
    
    rg = Graph<std::string>({"A",   "B", "C",},
                            {{1, 2}, {},  {},});
    
    EXPECT_EQ(tg, rg);
    
    tg.add("D", "B");
    
    rg = Graph<std::string>({"A",    "B", "C", "D",},
                            {{1, 2}, {3}, {},  {}, });
    
    EXPECT_EQ(tg, rg);
    
    tg.add("E", "B");
    
    rg = Graph<std::string>({"A",    "B",    "C", "D", "E",},
                            {{1, 2}, {3, 4},  {},  {},  {},});
    
    EXPECT_EQ(tg, rg);
    
    tg.add("E", "C");
    
    rg = Graph<std::string>({"A",    "B",    "C", "D", "E",},
                            {{1, 2}, {3, 4}, {4},  {},  {},});
    
    EXPECT_EQ(tg, rg);
}

TEST_F(GraphTest, Invert)
{
    const auto tg = Graph<std::string>({"A", "B", "C", "D", "E",},
                                       {{1, 2}, {3, 4}, {4}, {}, {},});
    
    const auto rg = Graph<std::string>({"E", "D", "C", "B", "A",},
                                       {{2, 3}, {3}, {4}, {4}, {},});

    EXPECT_EQ(tg.invert(), rg);
    
    EXPECT_EQ(rg.invert(), tg);
}

