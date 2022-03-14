// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.
// E-mail: rinkevic@kth.se
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "test_fraction.hpp"

#include "fraction.hpp"

TEST_F(FractionTest, Constructor)
{
    EXPECT_EQ(Fraction(), Fraction(0, 0));
    
    EXPECT_EQ(Fraction(2), Fraction(2, 1));
    
    EXPECT_EQ(Fraction(-2), Fraction(-2, 1));
    
    EXPECT_EQ(Fraction(2, -1), Fraction(-2, 1));
    
    EXPECT_EQ(Fraction(1, 2), Fraction(3, 6));
}

TEST_F(FractionTest, OperatorEqual)
{
    EXPECT_TRUE(Fraction(3) == Fraction(6, 2));
}

TEST_F(FractionTest, OperatorNotEqual)
{
    EXPECT_TRUE(Fraction(3) != Fraction(5, 2));
}

TEST_F(FractionTest, OperatorLess)
{
    EXPECT_FALSE(Fraction(1, 3) < Fraction(1, 4));
    
    EXPECT_TRUE(Fraction(1, 5) < Fraction(1, 4));
}

TEST_F(FractionTest, OperatorAdd)
{
    const auto lhsfrac = Fraction(2, 3) + Fraction(1, 5);
    
    const auto rhsfrac = Fraction(13, 15);
    
    EXPECT_EQ(lhsfrac, rhsfrac);
}

TEST_F(FractionTest, OperatorSubstract)
{
    const auto lhsfrac = Fraction(2, 3) - Fraction(1, 5);
    
    const auto rhsfrac = Fraction(7, 15);
    
    EXPECT_EQ(lhsfrac, rhsfrac);
}

TEST_F(FractionTest, OperatorMultiply)
{
    const auto lhsfrac = Fraction(2, 3) * Fraction(3, 5);
    
    const auto rhsfrac = Fraction(2, 5);
    
    EXPECT_EQ(lhsfrac, rhsfrac);
}

TEST_F(FractionTest, OperatorDivide)
{
    const auto lhsfrac = Fraction(2, 3) / Fraction(3, 5);
    
    const auto rhsfrac = Fraction(10, 9);
    
    EXPECT_EQ(lhsfrac, rhsfrac);
}

TEST_F(FractionTest, Numerator)
{
    auto frac = Fraction(-2, 3);
    
    EXPECT_EQ(frac.numerator(), -2);
    
    frac = Fraction(2, -3);
    
    EXPECT_EQ(frac.numerator(), -2);
    
    frac = Fraction(2, 3);
    
    EXPECT_EQ(frac.numerator(), 2);
    
    frac = Fraction(-2, -3);
    
    EXPECT_EQ(frac.numerator(), 2);
}

TEST_F(FractionTest, Denominator)
{
    auto frac = Fraction(-2, 3);
    
    EXPECT_EQ(frac.denominator(), 3);
    
    frac = Fraction(2, -3);
    
    EXPECT_EQ(frac.denominator(), 3);
    
    frac = Fraction(2, 3);
    
    EXPECT_EQ(frac.denominator(), 3);
    
    frac = Fraction(-2, -3);
    
    EXPECT_EQ(frac.denominator(), 3);
}

TEST_F(FractionTest, IsNegative)
{
    auto frac = Fraction(-2, 3);
    
    EXPECT_TRUE(frac.is_negative());
    
    frac = Fraction(2, -3);
    
    EXPECT_TRUE(frac.is_negative());
    
    frac = Fraction(2, 3);
    
    EXPECT_FALSE(frac.is_negative());
    
    frac = Fraction(-2, -3);
    
    EXPECT_FALSE(frac.is_negative());
}

TEST_F(FractionTest, ToString)
{
    auto frac = Fraction(-2, 3);
    
    EXPECT_EQ(frac.to_string(), "-2/3");
    
    frac = Fraction(2, -3);
    
    EXPECT_EQ(frac.to_string(), "-2/3");
    
    frac = Fraction(2, 3);
    
    EXPECT_EQ(frac.to_string(), "2/3");
    
    frac = Fraction(-2, -3);
    
    EXPECT_EQ(frac.to_string(), "2/3");
}

TEST_F(FractionTest, Label)
{
    auto frac = Fraction(-2);
    
    EXPECT_EQ(frac.label(), "-2.0");
    
    frac = Fraction(2);
    
    EXPECT_EQ(frac.label(), "2.0");
    
    frac = Fraction(2, 3);
    
    EXPECT_EQ(frac.label(), "2.0 / 3.0");
    
    frac = Fraction(2, -3);
    
    EXPECT_EQ(frac.label(), "-2.0 / 3.0");
}
