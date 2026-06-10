// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v3i_rr2_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Three-center r*r^2 integral (operator "GR.R2(r)", rank-1 tensor) over the
// two-center types.
I2CIntegral rr2(int a, int b)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("GR.R2(r)", Tensor(1)), 0, {});
}

// r^2 integral (operator "GR2(r)").
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

TEST(V3IRR2DriverTest, IsRR2)
{
    const V3IRR2Driver drv;

    EXPECT_TRUE(drv.is_rr2(rr2(1, 0)));

    // r^2 operator and plain overlap operator are not the r*r^2 operator.
    EXPECT_FALSE(drv.is_rr2(r2(1, 0)));
    EXPECT_FALSE(drv.is_rr2(ovl(1, 0)));
}

TEST(V3IRR2DriverTest, AuxVrrOnPShell)
{
    const V3IRR2Driver drv;

    // (P|P): three GR2(r) terms (self, lowered-bra, lowered-ket) plus the
    // three matching G(r) terms.
    const auto tints = drv.aux_vrr(rr2(1, 1));

    ASSERT_EQ(tints.size(), 6u);

    EXPECT_TRUE(contains(tints, r2(1, 1)));
    EXPECT_TRUE(contains(tints, r2(0, 1)));
    EXPECT_TRUE(contains(tints, r2(1, 0)));

    EXPECT_TRUE(contains(tints, ovl(1, 1)));
    EXPECT_TRUE(contains(tints, ovl(0, 1)));
    EXPECT_TRUE(contains(tints, ovl(1, 0)));
}

TEST(V3IRR2DriverTest, AuxVrrOnSShell)
{
    const V3IRR2Driver drv;

    // (S|S): only the two self terms survive (no lowering possible).
    const auto tints = drv.aux_vrr(rr2(0, 0));

    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, r2(0, 0)));
    EXPECT_TRUE(contains(tints, ovl(0, 0)));
}

TEST(V3IRR2DriverTest, AuxVrrRejectsNonRR2)
{
    const V3IRR2Driver drv;

    EXPECT_TRUE(drv.aux_vrr(r2(1, 0)).empty());
}
