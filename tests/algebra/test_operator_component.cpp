// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "operator_component.hpp"

TEST(OperatorComponentTest, DefaultConstruction)
{
    const OperatorComponent op;

    EXPECT_EQ(op.name(), "");
    EXPECT_EQ(op.shape(), TensorComponent(0, 0, 0));
    EXPECT_EQ(op.target(), "none");
    EXPECT_EQ(op.center(), -1);
}

TEST(OperatorComponentTest, FullConstructionAndAccessors)
{
    const OperatorComponent op("d/dA", TensorComponent(1, 0, 0), "bra", 0);

    EXPECT_EQ(op.name(), "d/dA");
    EXPECT_EQ(op.shape(), TensorComponent(1, 0, 0));
    EXPECT_EQ(op.target(), "bra");
    EXPECT_EQ(op.center(), 0);
    EXPECT_EQ(op['x'], 1);
    EXPECT_EQ(op['y'], 0);
}

TEST(OperatorComponentTest, Equality)
{
    const OperatorComponent op("A", TensorComponent(1, 0, 0), "bra", 0);

    EXPECT_EQ(op, OperatorComponent("A", TensorComponent(1, 0, 0), "bra", 0));
    EXPECT_NE(op, OperatorComponent("B", TensorComponent(1, 0, 0), "bra", 0));
    EXPECT_NE(op, OperatorComponent("A", TensorComponent(0, 1, 0), "bra", 0));
    EXPECT_NE(op, OperatorComponent("A", TensorComponent(1, 0, 0), "ket", 0));
    EXPECT_NE(op, OperatorComponent("A", TensorComponent(1, 0, 0), "bra", 1));
}

TEST(OperatorComponentTest, OrderingPrecedence)
{
    // name dominates, then shape, then target, then center.
    EXPECT_TRUE(OperatorComponent("A") < OperatorComponent("B"));
    EXPECT_TRUE(OperatorComponent("A", TensorComponent(0, 0, 0)) <
                OperatorComponent("A", TensorComponent(1, 0, 0)));
    EXPECT_TRUE(OperatorComponent("A", TensorComponent(1, 0, 0), "bra") <
                OperatorComponent("A", TensorComponent(1, 0, 0), "ket"));
    EXPECT_TRUE(OperatorComponent("A", TensorComponent(1, 0, 0), "bra", 0) <
                OperatorComponent("A", TensorComponent(1, 0, 0), "bra", 1));
}

TEST(OperatorComponentTest, SimilarityIgnoresAxialSplit)
{
    const OperatorComponent op("A", TensorComponent(1, 0, 0), "bra", 0);

    EXPECT_TRUE(op.similar(OperatorComponent("A", TensorComponent(0, 1, 0), "bra", 0)));
    EXPECT_FALSE(op.similar(OperatorComponent("A", TensorComponent(2, 0, 0), "bra", 0)));
    EXPECT_FALSE(op.similar(OperatorComponent("A", TensorComponent(1, 0, 0), "ket", 0)));
    EXPECT_FALSE(op.similar(OperatorComponent("A", TensorComponent(1, 0, 0), "bra", 1)));
}

TEST(OperatorComponentTest, ToString)
{
    const OperatorComponent op("d/dA", TensorComponent(1, 0, 0), "bra", 0);

    EXPECT_EQ(op.to_string(), "{d/dA:(1,0,0)}[bra:0]");
}

TEST(OperatorComponentTest, ShiftRaisesShape)
{
    const OperatorComponent op("A", TensorComponent(1, 0, 0));

    const auto raised = op.shift('x', 1);
    ASSERT_TRUE(raised.has_value());
    EXPECT_EQ(raised->shape(), TensorComponent(2, 0, 0));
}

TEST(OperatorComponentTest, ShiftToScalarDroppedWhenNoscalar)
{
    const OperatorComponent op("A", TensorComponent(1, 0, 0));

    EXPECT_FALSE(op.shift('x', -1, /*noscalar=*/true).has_value());

    const auto kept = op.shift('x', -1, /*noscalar=*/false);
    ASSERT_TRUE(kept.has_value());
    EXPECT_EQ(kept->shape(), TensorComponent(0, 0, 0));
}

TEST(OperatorComponentTest, ShiftBelowZeroIsNullopt)
{
    EXPECT_FALSE(OperatorComponent("A", TensorComponent(0, 0, 0)).shift('x', -1).has_value());
}
