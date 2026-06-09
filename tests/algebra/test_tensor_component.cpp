// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "tensor_component.hpp"

TEST(TensorComponentTest, OrderIsSumOfAxialValues)
{
    EXPECT_EQ(TensorComponent(2, 1, 0).order(), 3);
    EXPECT_EQ(TensorComponent(0, 0, 0).order(), 0);
}

TEST(TensorComponentTest, AxisAccess)
{
    const TensorComponent tcomp(2, 1, 0);

    EXPECT_EQ(tcomp['x'], 2);
    EXPECT_EQ(tcomp['y'], 1);
    EXPECT_EQ(tcomp['z'], 0);
    EXPECT_EQ(tcomp['w'], -1);  // invalid axis
}

TEST(TensorComponentTest, Equality)
{
    EXPECT_EQ(TensorComponent(1, 1, 0), TensorComponent(1, 1, 0));
    EXPECT_NE(TensorComponent(1, 1, 0), TensorComponent(2, 0, 0));
}

TEST(TensorComponentTest, LexicographicOrdering)
{
    EXPECT_TRUE(TensorComponent(0, 0, 0) < TensorComponent(1, 0, 0));
    EXPECT_TRUE(TensorComponent(1, 0, 0) < TensorComponent(1, 1, 0));
    EXPECT_TRUE(TensorComponent(1, 1, 0) < TensorComponent(1, 1, 1));
    EXPECT_FALSE(TensorComponent(1, 0, 0) < TensorComponent(1, 0, 0));
}

TEST(TensorComponentTest, SimilarMeansSameOrder)
{
    EXPECT_TRUE(TensorComponent(1, 1, 0).similar(TensorComponent(2, 0, 0)));
    EXPECT_FALSE(TensorComponent(1, 0, 0).similar(TensorComponent(2, 0, 0)));
}

TEST(TensorComponentTest, MaximumAndPrimary)
{
    EXPECT_EQ(TensorComponent(2, 1, 0).maximum(), 2);
    EXPECT_EQ(TensorComponent(0, 3, 1).maximum(), 3);

    EXPECT_EQ(TensorComponent(0, 2, 0).primary(), 'y');
    EXPECT_EQ(TensorComponent(0, 0, 0).primary(), 'x');  // default
}

TEST(TensorComponentTest, Label)
{
    EXPECT_EQ(TensorComponent(0, 0, 0).label(), "0");
    EXPECT_EQ(TensorComponent(1, 0, 0).label(), "x");
    EXPECT_EQ(TensorComponent(2, 1, 0).label(), "xxy");
    EXPECT_EQ(TensorComponent(1, 1, 1).label(), "xyz");
}

TEST(TensorComponentTest, ToString)
{
    EXPECT_EQ(TensorComponent(2, 1, 0).to_string(), "(2,1,0)");
}

TEST(TensorComponentTest, ShiftRaisesAndLowers)
{
    const auto raised = TensorComponent(1, 0, 0).shift('x', 1);

    ASSERT_TRUE(raised.has_value());
    EXPECT_EQ(*raised, TensorComponent(2, 0, 0));

    const auto lowered = TensorComponent(1, 0, 0).shift('x', -1);

    ASSERT_TRUE(lowered.has_value());
    EXPECT_EQ(*lowered, TensorComponent(0, 0, 0));
}

TEST(TensorComponentTest, ShiftBelowZeroReturnsNullopt)
{
    EXPECT_FALSE(TensorComponent(0, 0, 0).shift('x', -1).has_value());
    EXPECT_FALSE(TensorComponent(0, 1, 0).shift('z', -1).has_value());
}
