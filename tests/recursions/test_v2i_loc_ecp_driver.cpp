// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v2i_loc_ecp_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Two-center local-ECP integral (operator "U_L") with momenta a and b.
I2CIntegral ecp(int a, int b, int order = 0)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("U_L"), order, {});
}

I2CIntegral ovl(int a, int b)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("1"), 0, {});
}

bool contains(const SI2CIntegrals& set, const I2CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V2ILocalECPDriverTest, IsLocalEcp)
{
    const V2ILocalECPDriver drv;

    EXPECT_TRUE(drv.is_local_ecp(ecp(1, 0)));
    // Wrong operator (overlap) is rejected.
    EXPECT_FALSE(drv.is_local_ecp(ovl(1, 0)));
}

// full_bra_vrr on (P|S): tval = (S|S); the two lower terms vanish (negative).
TEST(V2ILocalECPDriverTest, FullBraVrrOnPShell)
{
    const V2ILocalECPDriver drv;

    const auto tints = drv.full_bra_vrr(ecp(1, 0));

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints, ecp(0, 0)));
}

// full_bra_vrr on (D|S): tval = (P|S); shift(-1,0) adds (S|S); the ket-lower term
// vanishes (negative ket). Both produced integrals keep the U_L operator.
TEST(V2ILocalECPDriverTest, FullBraVrrOnDShell)
{
    const V2ILocalECPDriver drv;

    const auto tints = drv.full_bra_vrr(ecp(2, 0));

    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, ecp(1, 0)));
    EXPECT_TRUE(contains(tints, ecp(0, 0)));
    for (const auto& tint : tints)
    {
        EXPECT_EQ(tint.integrand(), Operator("U_L"));
    }
}

// plain bra_vrr keeps only the first two recursion terms of full_bra_vrr.
TEST(V2ILocalECPDriverTest, BraVrrOnDShell)
{
    const V2ILocalECPDriver drv;

    const auto tints = drv.bra_vrr(ecp(2, 0));

    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, ecp(1, 0)));
    EXPECT_TRUE(contains(tints, ecp(0, 0)));
}

// ket_vrr on (S|D): tval = (S|P); shift(-1,1) adds (S|S).
TEST(V2ILocalECPDriverTest, KetVrrOnDShell)
{
    const V2ILocalECPDriver drv;

    const auto tints = drv.ket_vrr(ecp(0, 2));

    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, ecp(0, 1)));
    EXPECT_TRUE(contains(tints, ecp(0, 0)));
}

// Guard: scalar (S|S) shell yields nothing to lower.
TEST(V2ILocalECPDriverTest, VrrOnScalarIsEmpty)
{
    const V2ILocalECPDriver drv;

    EXPECT_TRUE(drv.full_bra_vrr(ecp(0, 0)).empty());
    EXPECT_TRUE(drv.bra_vrr(ecp(0, 0)).empty());
    EXPECT_TRUE(drv.ket_vrr(ecp(0, 0)).empty());
}

// Guard: vrr on the wrong operator yields nothing.
TEST(V2ILocalECPDriverTest, VrrOnNonEcpIsEmpty)
{
    const V2ILocalECPDriver drv;

    EXPECT_TRUE(drv.full_bra_vrr(ovl(1, 0)).empty());
    EXPECT_TRUE(drv.bra_vrr(ovl(1, 0)).empty());
    EXPECT_TRUE(drv.ket_vrr(ovl(0, 1)).empty());
}

// create_full_recursion smoke test: non-empty, retains original, reduces to base.
TEST(V2ILocalECPDriverTest, CreateFullRecursionReducesToBase)
{
    const V2ILocalECPDriver drv;

    const auto tints = drv.create_full_recursion(SI2CIntegrals({ecp(1, 1)}));

    ASSERT_FALSE(tints.empty());
    EXPECT_TRUE(contains(tints, ecp(1, 1)));   // original retained
    EXPECT_TRUE(contains(tints, ecp(0, 0)));   // reduced (S|S) base
    for (const auto& tint : tints)
    {
        EXPECT_EQ(tint.integrand(), Operator("U_L"));
    }
}

// create_recursion smoke test (plain variant).
TEST(V2ILocalECPDriverTest, CreateRecursionReducesToBase)
{
    const V2ILocalECPDriver drv;

    const auto tints = drv.create_recursion(SI2CIntegrals({ecp(1, 1)}));

    ASSERT_FALSE(tints.empty());
    EXPECT_TRUE(contains(tints, ecp(1, 1)));
    EXPECT_TRUE(contains(tints, ecp(0, 0)));
}

// A non-ECP integral passes through create_recursion unchanged.
TEST(V2ILocalECPDriverTest, CreateRecursionPassesThroughNonEcp)
{
    const V2ILocalECPDriver drv;

    const auto tints = drv.create_recursion(SI2CIntegrals({ovl(1, 0)}));

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints, ovl(1, 0)));
}
