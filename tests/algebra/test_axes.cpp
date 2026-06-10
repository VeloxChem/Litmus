// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "axes.hpp"

TEST(AxesTest, ToIndex)
{
    EXPECT_EQ(axes::to_index('x'), 0);
    EXPECT_EQ(axes::to_index('y'), 1);
    EXPECT_EQ(axes::to_index('z'), 2);
}

TEST(AxesTest, ToIndexUnknownAxisIsMinusOne)
{
    EXPECT_EQ(axes::to_index('w'), -1);
    EXPECT_EQ(axes::to_index('0'), -1);
}

TEST(AxesTest, ToAxis)
{
    EXPECT_EQ(axes::to_axis(0), 'x');
    EXPECT_EQ(axes::to_axis(1), 'y');
    EXPECT_EQ(axes::to_axis(2), 'z');
}

TEST(AxesTest, ToAxisUnknownIndexIsZeroChar)
{
    EXPECT_EQ(axes::to_axis(3), '0');
    EXPECT_EQ(axes::to_axis(-1), '0');
}

TEST(AxesTest, RoundTrip)
{
    for (const char axis : std::string("xyz"))
    {
        EXPECT_EQ(axes::to_axis(axes::to_index(axis)), axis);
    }
}
