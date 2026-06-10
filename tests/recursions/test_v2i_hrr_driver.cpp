// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v2i_hrr_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Two-center overlap integral with bra/ket angular momenta a and b. The HRR
// transfer driver is operator-agnostic, so any integrand works here.
I2CIntegral ovl(int a, int b)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("1"), 0, {});
}

bool contains(const SI2CIntegrals& set, const I2CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V2IHRRDriverTest, BraHrrTransfersFromBraToKet)
{
    const V2IHRRDriver drv;

    // (P|S): the bra HRR lowers the bra (S|S) and raises the ket (S|P).
    const auto tints = drv.bra_hrr(ovl(1, 0));

    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, ovl(0, 0)));
    EXPECT_TRUE(contains(tints, ovl(0, 1)));
}

TEST(V2IHRRDriverTest, KetHrrTransfersFromKetToBra)
{
    const V2IHRRDriver drv;

    // (S|P): the ket HRR lowers the ket (S|S) and raises the bra (P|S).
    const auto tints = drv.ket_hrr(ovl(0, 1));

    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, ovl(0, 0)));
    EXPECT_TRUE(contains(tints, ovl(1, 0)));
}

TEST(V2IHRRDriverTest, BraHrrOnScalarBraIsEmpty)
{
    const V2IHRRDriver drv;

    // Nothing to lower on the bra of an (S|S) integral.
    EXPECT_TRUE(drv.bra_hrr(ovl(0, 0)).empty());
}

TEST(V2IHRRDriverTest, KetHrrOnScalarKetIsEmpty)
{
    const V2IHRRDriver drv;

    // Nothing to lower on the ket of an (S|S) integral.
    EXPECT_TRUE(drv.ket_hrr(ovl(0, 0)).empty());
}

TEST(V2IHRRDriverTest, BraHrrFromHigherBra)
{
    const V2IHRRDriver drv;

    // (D|P): lowers bra to (P|P) and raises ket to (P|D).
    const auto tints = drv.bra_hrr(ovl(2, 1));

    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, ovl(1, 1)));
    EXPECT_TRUE(contains(tints, ovl(1, 2)));
}
