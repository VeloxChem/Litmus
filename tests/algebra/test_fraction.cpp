// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "fraction.hpp"

TEST(FractionTest, DefaultConstructorIsNan)
{
    const Fraction frac;

    EXPECT_TRUE(frac.isNan());
    EXPECT_EQ(frac.numerator(), 0);
    EXPECT_EQ(frac.denominator(), 0);
}

TEST(FractionTest, IntegerConstructorHasUnitDenominator)
{
    const Fraction frac(5);

    EXPECT_FALSE(frac.isNan());
    EXPECT_EQ(frac.numerator(), 5);
    EXPECT_EQ(frac.denominator(), 1);
}

TEST(FractionTest, ConstructorReducesToLowestTerms)
{
    const Fraction frac(2, 4);

    EXPECT_EQ(frac.numerator(), 1);
    EXPECT_EQ(frac.denominator(), 2);
}

TEST(FractionTest, NegativeDenominatorIsNormalizedToNumerator)
{
    const Fraction frac(3, -6);

    EXPECT_EQ(frac.numerator(), -1);
    EXPECT_EQ(frac.denominator(), 2);
}

TEST(FractionTest, DoubleNegativeReducesToPositive)
{
    const Fraction frac(-2, -4);

    EXPECT_EQ(frac.numerator(), 1);
    EXPECT_EQ(frac.denominator(), 2);
}

TEST(FractionTest, Addition)
{
    const auto frac = Fraction(1, 2) + Fraction(1, 3);

    EXPECT_EQ(frac.numerator(), 5);
    EXPECT_EQ(frac.denominator(), 6);
}

TEST(FractionTest, SubtractionToZeroHasUnitDenominator)
{
    const auto frac = Fraction(1, 2) - Fraction(1, 2);

    EXPECT_EQ(frac.numerator(), 0);
    EXPECT_EQ(frac.denominator(), 1);
    EXPECT_FALSE(frac.isNan());
}

TEST(FractionTest, MultiplicationReduces)
{
    const auto frac = Fraction(2, 3) * Fraction(3, 4);

    EXPECT_EQ(frac.numerator(), 1);
    EXPECT_EQ(frac.denominator(), 2);
}

TEST(FractionTest, Division)
{
    const auto frac = Fraction(1, 2) / Fraction(3, 4);

    EXPECT_EQ(frac.numerator(), 2);
    EXPECT_EQ(frac.denominator(), 3);
}

TEST(FractionTest, AbsoluteValue)
{
    EXPECT_EQ(Fraction(-3, 4).abs(), Fraction(3, 4));
    EXPECT_EQ(Fraction(3, -4).abs(), Fraction(3, 4));
}

TEST(FractionTest, IsNegative)
{
    EXPECT_TRUE(Fraction(-1, 2).is_negative());
    EXPECT_FALSE(Fraction(1, 2).is_negative());
    EXPECT_TRUE(Fraction(1, -2).is_negative());
}

TEST(FractionTest, EqualityAfterReduction)
{
    EXPECT_EQ(Fraction(1, 2), Fraction(2, 4));
    EXPECT_NE(Fraction(1, 2), Fraction(1, 3));
}

TEST(FractionTest, LessThanOrdering)
{
    EXPECT_TRUE(Fraction(1, 3) < Fraction(1, 2));
    EXPECT_FALSE(Fraction(1, 2) < Fraction(1, 3));
    EXPECT_TRUE(Fraction(-1, 2) < Fraction(1, 2));
}

TEST(FractionTest, Label)
{
    EXPECT_EQ(Fraction(3).label(), "3.0");
    EXPECT_EQ(Fraction(1, 2).label(), "1.0 / 2.0");
}
