// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "tensor.hpp"

TEST(TensorTest, Order)
{
    EXPECT_EQ(Tensor(0).order(), 0);
    EXPECT_EQ(Tensor(3).order(), 3);
}

TEST(TensorTest, ConstructFromComponentUsesItsOrder)
{
    const Tensor tensor(TensorComponent(1, 1, 0));

    EXPECT_EQ(tensor.order(), 2);
}

TEST(TensorTest, ShellLabels)
{
    EXPECT_EQ(Tensor(0).label(), "S");
    EXPECT_EQ(Tensor(1).label(), "P");
    EXPECT_EQ(Tensor(2).label(), "D");
    EXPECT_EQ(Tensor(3).label(), "F");
}

TEST(TensorTest, LabelBoundaryAtSixteen)
{
    // names string "SPDFGHIKLMNOQRTUV" has 17 entries (indices 0..16).
    EXPECT_EQ(Tensor(16).label(), "V");
    EXPECT_EQ(Tensor(17).label(), "l17");
}

TEST(TensorTest, ComponentCounts)
{
    // A shell of order n has (n + 1)(n + 2) / 2 Cartesian components.
    EXPECT_EQ(Tensor(0).components().size(), 1u);
    EXPECT_EQ(Tensor(1).components().size(), 3u);
    EXPECT_EQ(Tensor(2).components().size(), 6u);
    EXPECT_EQ(Tensor(3).components().size(), 10u);
}

TEST(TensorTest, Ordering)
{
    EXPECT_TRUE(Tensor(1) < Tensor(2));
    EXPECT_FALSE(Tensor(2) < Tensor(1));
}

TEST(TensorTest, Equality)
{
    EXPECT_EQ(Tensor(2), Tensor(TensorComponent(1, 1, 0)));
    EXPECT_NE(Tensor(1), Tensor(2));
}
