// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "factor.hpp"

TEST(FactorTest, Equality)
{
    const Factor a("PA", "PA", TensorComponent(1, 0, 0));
    const Factor b("PA", "PA", TensorComponent(1, 0, 0));
    const Factor c("PB", "PB", TensorComponent(1, 0, 0));

    EXPECT_EQ(a, b);
    EXPECT_NE(a, c);
}

TEST(FactorTest, OrderByName)
{
    const Factor a("a", "lbl", TensorComponent(0, 0, 0));
    const Factor b("b", "lbl", TensorComponent(0, 0, 0));

    EXPECT_TRUE(a < b);
    EXPECT_FALSE(b < a);
}

// Regression test: when two factors share a name but differ in label,
// operator< must order them by label. A prior bug compared shapes in this
// branch, so two factors differing only in label compared as equivalent
// (neither < the other), violating strict weak ordering in std::set/std::map.
TEST(FactorTest, OrderByLabelWhenNameMatches)
{
    const Factor alpha("w", "alpha", TensorComponent(0, 0, 0));
    const Factor beta("w", "beta", TensorComponent(0, 0, 0));

    EXPECT_NE(alpha, beta);

    // Exactly one ordering holds (strict weak ordering), following label order.
    EXPECT_TRUE(alpha < beta);
    EXPECT_FALSE(beta < alpha);
}

TEST(FactorTest, OrderByShapeWhenNameAndLabelMatch)
{
    const Factor lo("x", "x", TensorComponent(1, 0, 0));
    const Factor hi("x", "x", TensorComponent(2, 0, 0));

    EXPECT_TRUE(lo < hi);
    EXPECT_FALSE(hi < lo);
}

TEST(FactorTest, Order)
{
    EXPECT_EQ(Factor("p", "p", TensorComponent(2, 1, 0)).order(), 3);
    EXPECT_EQ(Factor("s", "s", TensorComponent(0, 0, 0)).order(), 0);
}
