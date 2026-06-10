// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v2i_ovl_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Two-center overlap integral with bra/ket angular momenta a and b.
I2CIntegral ovl(int a, int b)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("1"), 0, {});
}

bool contains(const SI2CIntegrals& set, const I2CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V2IOverlapDriverTest, IsOverlap)
{
    const V2IOverlapDriver drv;

    EXPECT_TRUE(drv.is_overlap(ovl(1, 0)));
    // Non-overlap operator is rejected.
    EXPECT_FALSE(drv.is_overlap(
        I2CIntegral(OneCenter("a", 1), OneCenter("b", 0), Operator("p", Tensor(1)), 0, {})));
}

TEST(V2IOverlapDriverTest, BraVrrLowersBra)
{
    const V2IOverlapDriver drv;

    // (P|S) reduces to a single (S|S) base integral.
    const auto tints = drv.bra_vrr(ovl(1, 0));

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints, ovl(0, 0)));
}

TEST(V2IOverlapDriverTest, BraVrrOnScalarIsEmpty)
{
    const V2IOverlapDriver drv;

    // Nothing to lower on an (S|S) integral.
    EXPECT_TRUE(drv.bra_vrr(ovl(0, 0)).empty());
}

TEST(V2IOverlapDriverTest, KetVrrLowersKet)
{
    const V2IOverlapDriver drv;

    const auto tints = drv.ket_vrr(ovl(0, 1));

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints, ovl(0, 0)));
}

TEST(V2IOverlapDriverTest, CreateRecursionReducesToBases)
{
    const V2IOverlapDriver drv;

    // (P|S) expands to itself plus the (S|S) base integral.
    const auto tints = drv.create_recursion(SI2CIntegrals({ovl(1, 0)}));

    EXPECT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, ovl(1, 0)));
    EXPECT_TRUE(contains(tints, ovl(0, 0)));

    // Every produced integral is an overlap integral.
    for (const auto& tint : tints)
    {
        EXPECT_EQ(tint.integrand(), Operator("1"));
    }
}

TEST(V2IOverlapDriverTest, CreateRecursionOnBaseIsItself)
{
    const V2IOverlapDriver drv;

    const auto tints = drv.create_recursion(SI2CIntegrals({ovl(0, 0)}));

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints, ovl(0, 0)));
}
