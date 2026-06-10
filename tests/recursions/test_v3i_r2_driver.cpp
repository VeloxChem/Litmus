// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v3i_r2_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Three-center r^2 integral (operator "GR2(r)") over the two-center types.
I2CIntegral r2(int a, int b)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("GR2(r)"), 0, {});
}

// Plain three-center overlap integral (operator "G(r)").
I2CIntegral ovl(int a, int b)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("G(r)"), 0, {});
}

bool contains(const SI2CIntegrals& set, const I2CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V3IR2DriverTest, IsR2)
{
    const V3IR2Driver drv;

    EXPECT_TRUE(drv.is_r2(r2(1, 0)));

    // Plain overlap operator "G(r)" is not the r^2 operator.
    EXPECT_FALSE(drv.is_r2(ovl(1, 0)));
    EXPECT_FALSE(drv.is_r2(
        I2CIntegral(OneCenter("a", 1), OneCenter("b", 0), Operator("1"), 0, {})));
}

TEST(V3IR2DriverTest, AuxVrrOnPShell)
{
    const V3IR2Driver drv;

    // (P|P): G(r), lowered-bra (S|P), lowered-ket (P|S), and bra+ket lowered (S|S).
    const auto tints = drv.aux_vrr(r2(1, 1));

    ASSERT_EQ(tints.size(), 4u);
    EXPECT_TRUE(contains(tints, ovl(1, 1)));
    EXPECT_TRUE(contains(tints, ovl(0, 1)));
    EXPECT_TRUE(contains(tints, ovl(1, 0)));
    EXPECT_TRUE(contains(tints, ovl(0, 0)));
}

TEST(V3IR2DriverTest, AuxVrrOnDShellAddsTwiceLoweredTerms)
{
    const V3IR2Driver drv;

    // (D|D) additionally activates the shift(-2, .) terms.
    const auto tints = drv.aux_vrr(r2(2, 2));

    ASSERT_EQ(tints.size(), 6u);
    EXPECT_TRUE(contains(tints, ovl(2, 2)));
    EXPECT_TRUE(contains(tints, ovl(0, 2)));  // shift(-2, 0)
    EXPECT_TRUE(contains(tints, ovl(2, 0)));  // shift(-2, 1)
}

TEST(V3IR2DriverTest, AuxVrrRejectsNonR2)
{
    const V3IR2Driver drv;

    EXPECT_TRUE(drv.aux_vrr(ovl(1, 0)).empty());
}
