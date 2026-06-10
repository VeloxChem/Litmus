// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v3i_ovl_grad_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Three-center overlap-gradient integral (operator "GX(r)", rank-1 tensor)
// over the two-center types.
I2CIntegral grad(int a, int b)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("GX(r)", Tensor(1)), 0, {});
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

TEST(V3IOverlapGradientDriverTest, IsOverlapGradient)
{
    const V3IOverlapGradientDriver drv;

    EXPECT_TRUE(drv.is_overlap_gradient(grad(1, 0)));

    // Plain overlap operator "G(r)" is not the gradient operator.
    EXPECT_FALSE(drv.is_overlap_gradient(ovl(1, 0)));
    // Scalar operator "1" is not the gradient operator either.
    EXPECT_FALSE(drv.is_overlap_gradient(
        I2CIntegral(OneCenter("a", 1), OneCenter("b", 0), Operator("1"), 0, {})));
}

TEST(V3IOverlapGradientDriverTest, AuxVrrOnPShellBra)
{
    const V3IOverlapGradientDriver drv;

    // (P|S): the G(r) integral and its lowered-bra (S|S) partner.
    const auto tints = drv.aux_vrr(grad(1, 0));

    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, ovl(1, 0)));
    EXPECT_TRUE(contains(tints, ovl(0, 0)));
}

TEST(V3IOverlapGradientDriverTest, AuxVrrOnPShellBraAndKet)
{
    const V3IOverlapGradientDriver drv;

    // (P|P): G(r) plus lowered-bra (S|P) and lowered-ket (P|S).
    const auto tints = drv.aux_vrr(grad(1, 1));

    ASSERT_EQ(tints.size(), 3u);
    EXPECT_TRUE(contains(tints, ovl(1, 1)));
    EXPECT_TRUE(contains(tints, ovl(0, 1)));
    EXPECT_TRUE(contains(tints, ovl(1, 0)));
}

TEST(V3IOverlapGradientDriverTest, AuxVrrRejectsNonGradient)
{
    const V3IOverlapGradientDriver drv;

    EXPECT_TRUE(drv.aux_vrr(ovl(1, 0)).empty());
}
