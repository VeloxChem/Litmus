// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "operator.hpp"

TEST(OperatorTest, DefaultConstruction)
{
    const Operator op;

    EXPECT_EQ(op.name(), "");
    EXPECT_EQ(op.shape().order(), 0);
}

TEST(OperatorTest, NameAndShape)
{
    const Operator op("Overlap", Tensor(1));

    EXPECT_EQ(op.name(), "Overlap");
    EXPECT_EQ(op.shape().order(), 1);
    EXPECT_EQ(op.label(), "P");  // label delegates to the tensor shell label
}

TEST(OperatorTest, Equality)
{
    EXPECT_EQ(Operator("A", Tensor(0)), Operator("A", Tensor(0)));
    EXPECT_NE(Operator("A", Tensor(0)), Operator("A", Tensor(1)));
    EXPECT_NE(Operator("A", Tensor(0)), Operator("B", Tensor(0)));
}

TEST(OperatorTest, OrderingByNameThenShape)
{
    EXPECT_TRUE(Operator("A", Tensor(0)) < Operator("B", Tensor(0)));
    EXPECT_TRUE(Operator("A", Tensor(0)) < Operator("A", Tensor(1)));
    EXPECT_FALSE(Operator("A", Tensor(1)) < Operator("A", Tensor(0)));
}

TEST(OperatorTest, ShiftRaisesAndLowersShape)
{
    const auto raised = Operator("X", Tensor(1)).shift(1);

    ASSERT_TRUE(raised.has_value());
    EXPECT_EQ(raised->shape().order(), 2);

    EXPECT_FALSE(Operator("X", Tensor(1)).shift(-2).has_value());
}

TEST(OperatorTest, ComponentCountMatchesShape)
{
    EXPECT_EQ(Operator("X", Tensor(0)).components().size(), 1u);
    EXPECT_EQ(Operator("X", Tensor(1)).components().size(), 3u);
}
