// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "one_center.hpp"

TEST(OneCenterTest, DefaultConstruction)
{
    const OneCenter center;

    EXPECT_EQ(center.shape(), 0);
    EXPECT_EQ(center.centers(), 1);
}

TEST(OneCenterTest, ConstructFromAngularMomentum)
{
    const OneCenter center("a", 2);

    EXPECT_EQ(center.shape(), 2);
    EXPECT_EQ(center[0], 2);  // operator[] returns the tensor order
}

TEST(OneCenterTest, ConstructFromComponentUsesItsOrder)
{
    const OneCenter center(OneCenterComponent("a", TensorComponent(1, 1, 0)));

    EXPECT_EQ(center.shape(), 2);
    EXPECT_EQ(center.label(), "D");
}

TEST(OneCenterTest, Equality)
{
    EXPECT_EQ(OneCenter("a", 1), OneCenter("a", Tensor(1)));
    EXPECT_NE(OneCenter("a", 1), OneCenter("a", 2));
    EXPECT_NE(OneCenter("a", 1), OneCenter("b", 1));
}

TEST(OneCenterTest, OrderingByNameThenShape)
{
    EXPECT_TRUE(OneCenter("a", 5) < OneCenter("b", 0));
    EXPECT_TRUE(OneCenter("a", 1) < OneCenter("a", 2));
    EXPECT_FALSE(OneCenter("a", 2) < OneCenter("a", 1));
}

TEST(OneCenterTest, Label)
{
    EXPECT_EQ(OneCenter("a", 0).label(), "S");
    EXPECT_EQ(OneCenter("a", 1).label(), "P");
}

TEST(OneCenterTest, ToString)
{
    EXPECT_EQ(OneCenter("a", 1).to_string(), "{a:(1)}");
}

TEST(OneCenterTest, ShiftOnlyAffectsCenterZero)
{
    const OneCenter center("a", 1);

    const auto raised = center.shift(1, 0);
    ASSERT_TRUE(raised.has_value());
    EXPECT_EQ(raised->shape(), 2);

    EXPECT_FALSE(center.shift(-2, 0).has_value());  // below zero
    EXPECT_FALSE(center.shift(1, 1).has_value());   // no second center
}

TEST(OneCenterTest, ComponentCountMatchesShell)
{
    EXPECT_EQ(OneCenter("a", 0).components().size(), 1u);
    EXPECT_EQ(OneCenter("a", 1).components().size(), 3u);
    EXPECT_EQ(OneCenter("a", 2).components().size(), 6u);
}
