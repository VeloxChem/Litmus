// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "two_center_pair.hpp"

TEST(TwoCenterPairTest, DefaultConstruction)
{
    const TwoCenterPair pair;

    EXPECT_EQ(pair.centers(), 2);
    EXPECT_EQ(pair[0], 0);
    EXPECT_EQ(pair[1], 0);
}

TEST(TwoCenterPairTest, AngularMomentumConstructionAndIndexing)
{
    const TwoCenterPair pair("a", 1, "b", 2);

    EXPECT_EQ(pair[0], 1);
    EXPECT_EQ(pair[1], 2);
}

TEST(TwoCenterPairTest, ConstructFromComponentUsesOrders)
{
    const TwoCenterPairComponent comp({"a", "b"},
                                      {TensorComponent(1, 1, 0), TensorComponent(0, 0, 1)});
    const TwoCenterPair pair(comp);

    EXPECT_EQ(pair[0], 2);
    EXPECT_EQ(pair[1], 1);
}

TEST(TwoCenterPairTest, Equality)
{
    const TwoCenterPair pair("a", 1, "b", 2);

    EXPECT_EQ(pair, TwoCenterPair("a", 1, "b", 2));
    EXPECT_NE(pair, TwoCenterPair("a", 1, "b", 3));
    EXPECT_NE(pair, TwoCenterPair("a", 1, "c", 2));
}

TEST(TwoCenterPairTest, OrderingByNamesThenShapes)
{
    EXPECT_TRUE(TwoCenterPair("a", 5, "b", 0) < TwoCenterPair("a", 0, "c", 0));  // names first
    EXPECT_TRUE(TwoCenterPair("a", 1, "b", 0) < TwoCenterPair("a", 1, "b", 1));  // then shapes
}

TEST(TwoCenterPairTest, Label)
{
    EXPECT_EQ(TwoCenterPair("a", 1, "b", 0).label(), "PS");
    EXPECT_EQ(TwoCenterPair("a", 2, "b", 1).label(), "DP");
}

TEST(TwoCenterPairTest, ToString)
{
    EXPECT_EQ(TwoCenterPair("a", 1, "b", 2).to_string(), "{a:(1);b:(2)}");
}

TEST(TwoCenterPairTest, ShiftTargetsSelectedCenter)
{
    const TwoCenterPair pair("a", 1, "b", 2);

    const auto first = pair.shift(1, 0);
    ASSERT_TRUE(first.has_value());
    EXPECT_EQ((*first)[0], 2);
    EXPECT_EQ((*first)[1], 2);

    const auto second = pair.shift(-1, 1);
    ASSERT_TRUE(second.has_value());
    EXPECT_EQ((*second)[1], 1);

    EXPECT_FALSE(pair.shift(-2, 0).has_value());  // below zero
    EXPECT_FALSE(pair.shift(1, 2).has_value());   // no third center
}

TEST(TwoCenterPairTest, ComponentCountIsCartesianProduct)
{
    // P (3) x P (3) = 9 paired components.
    EXPECT_EQ(TwoCenterPair("a", 1, "b", 1).components().size(), 9u);
    // S (1) x D (6) = 6.
    EXPECT_EQ(TwoCenterPair("a", 0, "b", 2).components().size(), 6u);
}
