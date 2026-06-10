// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "string_formater.hpp"

TEST(StringFormaterTest, LowercaseConvertsUppercase)
{
    EXPECT_EQ(fstr::lowercase("HELLO"), "hello");
    EXPECT_EQ(fstr::lowercase("MixedCase"), "mixedcase");
}

TEST(StringFormaterTest, LowercaseLeavesNonLettersUntouched)
{
    EXPECT_EQ(fstr::lowercase("Overlap_2c"), "overlap_2c");
    EXPECT_EQ(fstr::lowercase("X12+Y"), "x12+y");
}

TEST(StringFormaterTest, LowercaseIsIdempotentOnLowercase)
{
    EXPECT_EQ(fstr::lowercase("already"), "already");
}

TEST(StringFormaterTest, LowercaseEmptyString)
{
    EXPECT_EQ(fstr::lowercase(""), "");
}
