// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v2i_linmom_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Two-center linear-momentum integral (operator p) with momenta a and b.
I2CIntegral linmom(int a, int b)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("p", Tensor(1)), 0, {});
}

I2CIntegral ovl(int a, int b)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("1"), 0, {});
}

}  // namespace

TEST(V2ILinearMomentumDriverTest, IsLinmom)
{
    const V2ILinearMomentumDriver drv;

    EXPECT_TRUE(drv.is_linmom(linmom(0, 0)));
    EXPECT_FALSE(drv.is_linmom(ovl(0, 0)));
}

// Regression: op_vrr must lower the p operator to the overlap operator "1" on
// the ket-shifted terms. A discarded const replace() previously left the
// produced integrals carrying the linear-momentum operator p.
TEST(V2ILinearMomentumDriverTest, OpVrrProducesOverlapOperator)
{
    const V2ILinearMomentumDriver drv;

    const auto tints = drv.op_vrr(linmom(0, 1));

    ASSERT_FALSE(tints.empty());
    for (const auto& tint : tints)
    {
        EXPECT_EQ(tint.integrand(), Operator("1"));
        EXPECT_NE(tint.integrand(), Operator("p", Tensor(1)));
    }
}

TEST(V2ILinearMomentumDriverTest, OpVrrTermCounts)
{
    const V2ILinearMomentumDriver drv;

    // ket S: only the ket+1 term exists (ket-1 would be negative).
    EXPECT_EQ(drv.op_vrr(linmom(0, 0)).size(), 1u);

    // ket P: both ket+1 and ket-1 terms exist.
    EXPECT_EQ(drv.op_vrr(linmom(0, 1)).size(), 2u);
}

TEST(V2ILinearMomentumDriverTest, OpVrrOnNonLinmomIsEmpty)
{
    const V2ILinearMomentumDriver drv;

    EXPECT_TRUE(drv.op_vrr(ovl(0, 1)).empty());
}

// Smoke test for the apply/create chain (exercises apply_op_vrr, whose loop is
// now gated on is_linmom): it must terminate and yield the overlap bases.
TEST(V2ILinearMomentumDriverTest, CreateRecursionTerminatesAndYieldsOverlap)
{
    const V2ILinearMomentumDriver drv;

    const auto tints = drv.create_recursion(SI2CIntegrals({linmom(0, 1)}));

    ASSERT_FALSE(tints.empty());
    EXPECT_NE(tints.find(ovl(0, 0)), tints.cend());  // reduced overlap base present
}
