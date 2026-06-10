// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v3i_ovl_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Three-center overlap integral (operator "G(r)") over the two-center types.
I2CIntegral ovl(int a, int b)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("G(r)"), 0, {});
}

bool contains(const SI2CIntegrals& set, const I2CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V3IOverlapDriverTest, IsOverlap)
{
    const V3IOverlapDriver drv;

    EXPECT_TRUE(drv.is_overlap(ovl(1, 0)));
    // Plain two-center overlap operator "1" is not a three-center G(r) integral.
    EXPECT_FALSE(drv.is_overlap(
        I2CIntegral(OneCenter("a", 1), OneCenter("b", 0), Operator("1"), 0, {})));
}

TEST(V3IOverlapDriverTest, BraVrrLowersBra)
{
    const V3IOverlapDriver drv;

    const auto tints = drv.bra_vrr(ovl(1, 0));

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints, ovl(0, 0)));
}

TEST(V3IOverlapDriverTest, KetVrrLowersKet)
{
    const V3IOverlapDriver drv;

    const auto tints = drv.ket_vrr(ovl(0, 1));

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints, ovl(0, 0)));
}

TEST(V3IOverlapDriverTest, CreateRecursionReducesToBases)
{
    const V3IOverlapDriver drv;

    const auto tints = drv.create_recursion(SI2CIntegrals({ovl(1, 0)}));

    // The (P|S) G(r) integral, its (S|S) G(r) base, and the (S|S) plain-overlap
    // auxiliary produced by aux_vrr.
    EXPECT_EQ(tints.size(), 3u);
    EXPECT_TRUE(contains(tints, ovl(1, 0)));
    EXPECT_TRUE(contains(tints, ovl(0, 0)));

    const auto plain = I2CIntegral(OneCenter("a", 0), OneCenter("b", 0), Operator("1"), 0, {});
    EXPECT_TRUE(contains(tints, plain));
}
