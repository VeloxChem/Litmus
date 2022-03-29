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

TEST_F(GraphTest, OperatorBrackets)
{
    const auto tg = Graph<std::string>({"A",    "B",    "C", "D", "E",},
                                       {{1, 2}, {3, 4}, {4},  {},  {},});
    
    EXPECT_EQ(tg[0], "A");
    
    EXPECT_EQ(tg[1], "B");
    
    EXPECT_EQ(tg[2], "C");
    
    EXPECT_EQ(tg[3], "D");
    
    EXPECT_EQ(tg[4], "E");
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

TEST_F(GraphTest, Replace)
{
    auto tg = Graph<std::string>({"A", "B", "C", "D", "E",},
                                 {{1, 2}, {3, 4}, {4}, {}, {},});

    tg.replace("F", 0);
    
    tg.replace("X", 2);
    
    const auto rg = Graph<std::string>({"F", "B", "X", "D", "E",},
                                       {{1, 2}, {3, 4}, {4}, {}, {},});
    
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

TEST_F(GraphTest, Vertices)
{
    const auto tg = Graph<std::string>({"A", "B", "C", "D", "E",},
                                       {{1, 2}, {3, 4}, {4}, {}, {},});

    EXPECT_EQ(tg.vertices(), 5);
}

TEST_F(GraphTest, Edge)
{
    const auto tg = Graph<std::string>({"A", "B", "C", "D", "E",},
                                       {{1, 2}, {3, 4}, {4}, {}, {},});

    EXPECT_EQ(tg.edge(0), std::set<int>({1, 2}));
    
    EXPECT_EQ(tg.edge(1), std::set<int>({3, 4}));
    
    EXPECT_EQ(tg.edge(2), std::set<int>({4}));
    
    EXPECT_EQ(tg.edge(3), std::set<int>({}));
    
    EXPECT_EQ(tg.edge(4), std::set<int>({}));
    
}

TEST_F(GraphTest, Orpahns)
{
    const auto tg = Graph<std::string>({"A", "B", "C", "D", "E",},
                                       {{1, 2}, {3, 4}, {4}, {}, {},});
    
    const auto rg = Graph<std::string>({"E", "D", "C", "B", "A",},
                                       {{2, 3}, {3}, {4}, {4}, {},});

    EXPECT_EQ(tg.orphans(), std::vector<int>({3, 4,}));
    
    EXPECT_EQ(rg.orphans(), std::vector<int>({4,}));
}

TEST_F(GraphTest, Merge)
{
    auto tg = Graph<std::string>({"A", "B", "C", "D", "E",},
                                 {{1, 2, 3}, {2, 3, 4}, {3, 4}, {4}, {},});

    tg.merge(1, 3);
    
    EXPECT_EQ(tg.vertices(), 4);
    
    auto rg = Graph<std::string>({"A", "BD", "C", "E",},
                                 {{1, 2}, {2, 3}, {3}, {},});
    
    EXPECT_EQ(tg, rg);
}

TEST_F(GraphTest, Similar)
{
    auto tg = Graph<std::string>({"A", "B", "C", "B", "BB",},
                                 {{1, 2, 3}, {2, 3, 4}, {3, 4}, {4}, {},});

    tg.reduce();
    
    const auto rg = Graph<std::string>({"A", "BBBB", "C",},
                                       {{1, 2}, {2}, {},});
    
    EXPECT_EQ(tg, rg);
}

TEST_F(GraphTest, Indexes)
{
    const auto tg = Graph<std::string>({"A", "C", "B", "E", "D",},
                                       {{1, 2, 3}, {2, 3, 4}, {3, 4}, {4}, {},});

    EXPECT_EQ(tg.indexes<std::string>(), std::vector<int>({0, 2, 1, 4, 3}));
}

TEST_F(GraphTest, Sort)
{
    auto tg = Graph<std::string>({"A", "C", "B", "E", "D",},
                                 {{1, 2, 3}, {2, 3, 4}, {3, 4}, {4}, {},});
    
    tg.sort<std::string>(true);
    
    
    const auto nverts = tg.vertices();
    
    for (int i = 0; i < nverts; i++)
    {
        std::cout << "Vertice: " << i << std::endl;
    
        std::cout << "@value -> " << tg[i] << std::endl;
   
        std::cout << "@edges -> ";
        
        for (const auto& tval : tg.edge(i))
        {
            std::cout << tval << " ";
        }
    
        std::cout << std::endl;
    }
}

