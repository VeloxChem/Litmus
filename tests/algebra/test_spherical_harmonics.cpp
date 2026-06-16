// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <utility>

#include "spherical_harmonics.hpp"
#include "tensor.hpp"
#include "tensor_component.hpp"

using namespace sphar;

namespace {

/// Looks up the coefficient on a given Cartesian monomial in a harmonic, as a
/// numerical value (zero when the monomial is absent).
double
coeff_of(const VSphericalTerms& terms, const TensorComponent& monomial)
{
    for (const auto& [component, factor] : terms)
    {
        if (component == monomial) return factor.value();
    }

    return 0.0;
}

/// Maps a Cartesian component to its canonical index for order l (the VeloxChem
/// transformation-factor index).
int
cartesian_index(const int l, const TensorComponent& monomial)
{
    int index = 0;

    for (const auto& component : Tensor(l).components())
    {
        if (component == monomial) return index;

        index++;
    }

    return -1;
}

}  // namespace

TEST(SphericalHarmonicsTest, SFactorIsUnity)
{
    const auto terms = spherical_factors(0, 0);

    ASSERT_EQ(terms.size(), 1u);
    EXPECT_EQ(terms[0].first, TensorComponent(0, 0, 0));
    EXPECT_EQ(terms[0].second, SphericalFactor(Fraction(1), 1));
}

TEST(SphericalHarmonicsTest, POrderMatchesVeloxChem)
{
    // component = m + l: 0 -> y, 1 -> z, 2 -> x.
    EXPECT_EQ(spherical_component_factors(1, 0)[0].first, TensorComponent(0, 1, 0));
    EXPECT_EQ(spherical_component_factors(1, 1)[0].first, TensorComponent(0, 0, 1));
    EXPECT_EQ(spherical_component_factors(1, 2)[0].first, TensorComponent(1, 0, 0));

    for (int c = 0; c < 3; c++)
    {
        const auto terms = spherical_component_factors(1, c);

        ASSERT_EQ(terms.size(), 1u);
        EXPECT_EQ(terms[0].second, SphericalFactor(Fraction(1), 1));
    }
}

TEST(SphericalHarmonicsTest, DCoefficientsAreExact)
{
    // S(2,-2) = sqrt(3) xy
    EXPECT_EQ(spherical_factors(2, -2), VSphericalTerms({{TensorComponent(1, 1, 0), {Fraction(1), 3}}}));

    // S(2,-1) = sqrt(3) yz
    EXPECT_EQ(spherical_factors(2, -1), VSphericalTerms({{TensorComponent(0, 1, 1), {Fraction(1), 3}}}));

    // S(2,0) = -1/2 x^2 - 1/2 y^2 + z^2
    EXPECT_EQ(spherical_factors(2, 0),
              VSphericalTerms({{TensorComponent(2, 0, 0), {Fraction(-1, 2), 1}},
                               {TensorComponent(0, 2, 0), {Fraction(-1, 2), 1}},
                               {TensorComponent(0, 0, 2), {Fraction(1), 1}}}));

    // S(2,1) = sqrt(3) xz
    EXPECT_EQ(spherical_factors(2, 1), VSphericalTerms({{TensorComponent(1, 0, 1), {Fraction(1), 3}}}));

    // S(2,2) = sqrt(3)/2 (x^2 - y^2)
    EXPECT_EQ(spherical_factors(2, 2),
              VSphericalTerms({{TensorComponent(2, 0, 0), {Fraction(1, 2), 3}},
                               {TensorComponent(0, 2, 0), {Fraction(-1, 2), 3}}}));
}

TEST(SphericalHarmonicsTest, FMatchesVeloxChemValues)
{
    // VeloxChem F block, with f10 = 0.25 sqrt(10), f15 = sqrt(15), f6 = sqrt(6).
    const double f10 = 0.25 * std::sqrt(10.0);
    const double f15 = std::sqrt(15.0);
    const double f6  = std::sqrt(6.0);

    // component 0 (m=-3): {1, 3 f10}, {6, -f10}
    {
        const auto t = spherical_component_factors(3, 0);
        EXPECT_NEAR(coeff_of(t, TensorComponent(2, 1, 0)), 3.0 * f10, 1.0e-13);  // index 1: x^2 y
        EXPECT_NEAR(coeff_of(t, TensorComponent(0, 3, 0)), -f10, 1.0e-13);       // index 6: y^3
    }
    // component 2 (m=-1): {8, f6}, {1, -0.25 f6}, {6, -0.25 f6}
    {
        const auto t = spherical_component_factors(3, 2);
        EXPECT_NEAR(coeff_of(t, TensorComponent(0, 1, 2)), f6, 1.0e-13);          // index 8: y z^2
        EXPECT_NEAR(coeff_of(t, TensorComponent(2, 1, 0)), -0.25 * f6, 1.0e-13);  // index 1: x^2 y
        EXPECT_NEAR(coeff_of(t, TensorComponent(0, 3, 0)), -0.25 * f6, 1.0e-13);  // index 6: y^3
    }
    // component 3 (m=0): {9, 1.0}, {2, -1.5}, {7, -1.5}
    {
        const auto t = spherical_component_factors(3, 3);
        EXPECT_NEAR(coeff_of(t, TensorComponent(0, 0, 3)), 1.0, 1.0e-13);   // index 9: z^3
        EXPECT_NEAR(coeff_of(t, TensorComponent(2, 0, 1)), -1.5, 1.0e-13);  // index 2: x^2 z
        EXPECT_NEAR(coeff_of(t, TensorComponent(0, 2, 1)), -1.5, 1.0e-13);  // index 7: y^2 z
    }
    // component 5 (m=2): {2, 0.5 f15}, {7, -0.5 f15}
    {
        const auto t = spherical_component_factors(3, 5);
        EXPECT_NEAR(coeff_of(t, TensorComponent(2, 0, 1)), 0.5 * f15, 1.0e-13);   // index 2: x^2 z
        EXPECT_NEAR(coeff_of(t, TensorComponent(0, 2, 1)), -0.5 * f15, 1.0e-13);  // index 7: y^2 z
    }
}

TEST(SphericalHarmonicsTest, GMatchesVeloxChemValues)
{
    const double f35 = 0.5 * std::sqrt(35.0);
    const double f5  = 0.5 * std::sqrt(5.0);

    // component 4 (m=0): pure rationals.
    {
        const auto t = spherical_component_factors(4, 4);
        EXPECT_NEAR(coeff_of(t, TensorComponent(0, 0, 4)), 1.0, 1.0e-13);    // index 14: z^4
        EXPECT_NEAR(coeff_of(t, TensorComponent(4, 0, 0)), 0.375, 1.0e-13);  // index 0: x^4
        EXPECT_NEAR(coeff_of(t, TensorComponent(0, 4, 0)), 0.375, 1.0e-13);  // index 10: y^4
        EXPECT_NEAR(coeff_of(t, TensorComponent(2, 2, 0)), 0.75, 1.0e-13);   // index 3: x^2 y^2
        EXPECT_NEAR(coeff_of(t, TensorComponent(2, 0, 2)), -3.0, 1.0e-13);   // index 5: x^2 z^2
        EXPECT_NEAR(coeff_of(t, TensorComponent(0, 2, 2)), -3.0, 1.0e-13);   // index 12: y^2 z^2
    }
    // component 0 (m=-4): {1, f35}, {6, -f35}
    {
        const auto t = spherical_component_factors(4, 0);
        EXPECT_NEAR(coeff_of(t, TensorComponent(3, 1, 0)), f35, 1.0e-13);   // index 1: x^3 y
        EXPECT_NEAR(coeff_of(t, TensorComponent(1, 3, 0)), -f35, 1.0e-13);  // index 6: x y^3
    }
    // component 8 (m=4): {0, 0.25 f35}, {10, 0.25 f35}, {3, -1.5 f35}
    {
        const auto t = spherical_component_factors(4, 8);
        EXPECT_NEAR(coeff_of(t, TensorComponent(4, 0, 0)), 0.25 * f35, 1.0e-13);  // index 0: x^4
        EXPECT_NEAR(coeff_of(t, TensorComponent(0, 4, 0)), 0.25 * f35, 1.0e-13);  // index 10: y^4
        EXPECT_NEAR(coeff_of(t, TensorComponent(2, 2, 0)), -1.5 * f35, 1.0e-13);  // index 3: x^2 y^2
    }
    // component 6 (m=2): {5, 3 f5}, {12, -3 f5}, {0, -0.5 f5}, {10, 0.5 f5}
    {
        const auto t = spherical_component_factors(4, 6);
        EXPECT_NEAR(coeff_of(t, TensorComponent(2, 0, 2)), 3.0 * f5, 1.0e-13);   // index 5: x^2 z^2
        EXPECT_NEAR(coeff_of(t, TensorComponent(0, 2, 2)), -3.0 * f5, 1.0e-13);  // index 12: y^2 z^2
        EXPECT_NEAR(coeff_of(t, TensorComponent(4, 0, 0)), -0.5 * f5, 1.0e-13);  // index 0: x^4
        EXPECT_NEAR(coeff_of(t, TensorComponent(0, 4, 0)), 0.5 * f5, 1.0e-13);   // index 10: y^4
    }
}

TEST(SphericalHarmonicsTest, ComponentCountAndCanonicalIndices)
{
    for (int l = 0; l <= 5; l++)
    {
        // every spherical component yields a non-empty expansion ...
        for (int c = 0; c <= 2 * l; c++)
        {
            const auto terms = spherical_component_factors(l, c);

            ASSERT_FALSE(terms.empty()) << "l=" << l << " c=" << c;

            // ... whose monomials all have order l and a valid canonical index.
            for (const auto& [component, factor] : terms)
            {
                EXPECT_EQ(component.order(), l);
                EXPECT_GE(cartesian_index(l, component), 0);
                EXPECT_FALSE(factor.is_zero());
            }
        }

        // out-of-range indices give nothing.
        EXPECT_TRUE(spherical_component_factors(l, 2 * l + 1).empty());
    }
}

TEST(SphericalHarmonicsTest, OutOfRangeOrdersAreEmpty)
{
    EXPECT_TRUE(spherical_factors(-1, 0).empty());
    EXPECT_TRUE(spherical_factors(2, 3).empty());
    EXPECT_TRUE(spherical_factors(2, -3).empty());
}
