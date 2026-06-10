// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v2i_dip_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Two-center dipole integral (operator "r", rank-1) with bra/ket momenta a, b.
I2CIntegral dip(int a, int b)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("r", Tensor(1)), 0, {});
}

// Two-center overlap integral (the base of the dipole recursion).
I2CIntegral ovl(int a, int b)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("1"), 0, {});
}

bool contains(const SI2CIntegrals& set, const I2CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V2IDipoleDriverTest, IsDipole)
{
    const V2IDipoleDriver drv;

    EXPECT_TRUE(drv.is_dipole(dip(1, 0)));
    // Wrong operator (overlap) is rejected.
    EXPECT_FALSE(drv.is_dipole(ovl(1, 0)));
    // Prefixed integral is rejected.
    EXPECT_FALSE(drv.is_dipole(I2CIntegral(OneCenter("a", 1), OneCenter("b", 0),
                                           Operator("r", Tensor(1)), 0,
                                           {Operator("d/dA", Tensor(1))})));
}

TEST(V2IDipoleDriverTest, BraVrrLowersBraAndReducesToOverlap)
{
    const V2IDipoleDriver drv;

    // (P|r|S): lowers bra to (S|r|S) plus the overlap (S|1|S) base term.
    const auto tints = drv.bra_vrr(dip(1, 0));

    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, dip(0, 0)));
    EXPECT_TRUE(contains(tints, ovl(0, 0)));
}

TEST(V2IDipoleDriverTest, KetVrrLowersKetAndReducesToOverlap)
{
    const V2IDipoleDriver drv;

    // (S|r|P): lowers ket to (S|r|S) plus the overlap (S|1|S) base term.
    const auto tints = drv.ket_vrr(dip(0, 1));

    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, dip(0, 0)));
    EXPECT_TRUE(contains(tints, ovl(0, 0)));
}

TEST(V2IDipoleDriverTest, BraVrrOnWrongOperatorIsEmpty)
{
    const V2IDipoleDriver drv;

    // Overlap operator is not a dipole; nothing is produced.
    EXPECT_TRUE(drv.bra_vrr(ovl(1, 0)).empty());
}

TEST(V2IDipoleDriverTest, CreateRecursionContainsDipoleAndOverlapBases)
{
    const V2IDipoleDriver drv;

    const auto tints = drv.create_recursion(SI2CIntegrals({dip(1, 0)}));

    EXPECT_FALSE(tints.empty());
    // Original integral is preserved.
    EXPECT_TRUE(contains(tints, dip(1, 0)));
    // Fully reduced dipole base and overlap base are present.
    EXPECT_TRUE(contains(tints, dip(0, 0)));
    EXPECT_TRUE(contains(tints, ovl(0, 0)));
}

TEST(V2IDipoleDriverTest, CreateRecursionLeavesNonDipoleUntouched)
{
    const V2IDipoleDriver drv;

    // A non-dipole integral passes through unchanged.
    const auto tints = drv.create_recursion(SI2CIntegrals({ovl(1, 0)}));

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints, ovl(1, 0)));
}
