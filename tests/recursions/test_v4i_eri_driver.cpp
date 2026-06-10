// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v4i_eri_driver.hpp"
#include "t4c_defs.hpp"

namespace {

// Four-center ERI integral over bra pair (A, B) and ket pair (C, D).
I4CIntegral eri(int a, int b, int c, int d)
{
    return I4CIntegral(TwoCenterPair("a", a, "b", b),
                       TwoCenterPair("c", c, "d", d),
                       Operator("1/|r-r'|"), 0, {});
}

}  // namespace

TEST(V4IElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const V4IElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.is_electron_repulsion(eri(1, 0, 0, 0)));
    EXPECT_FALSE(drv.is_electron_repulsion(
        I4CIntegral(TwoCenterPair("a", 1, "b", 0), TwoCenterPair("c", 0, "d", 0), Operator("1"), 0, {})));
}

TEST(V4IElectronRepulsionDriverTest, BraHrrMovesAtoB)
{
    const V4IElectronRepulsionDriver drv;

    // (P,S|S,S): HRR lowers A and raises B -> two integrals.
    EXPECT_EQ(drv.bra_hrr(eri(1, 0, 0, 0)).size(), 2u);
}

TEST(V4IElectronRepulsionDriverTest, BraHrrOnScalarAIsEmpty)
{
    const V4IElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.bra_hrr(eri(0, 0, 0, 0)).empty());
}

TEST(V4IElectronRepulsionDriverTest, KetHrrMovesCtoD)
{
    const V4IElectronRepulsionDriver drv;

    // (S,S|P,S): HRR lowers C and raises D -> two integrals.
    EXPECT_EQ(drv.ket_hrr(eri(0, 0, 1, 0)).size(), 2u);
}

TEST(V4IElectronRepulsionDriverTest, KetHrrOnScalarCIsEmpty)
{
    const V4IElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.ket_hrr(eri(0, 0, 0, 0)).empty());
}

TEST(V4IElectronRepulsionDriverTest, HrrRejectsNonEriIntegral)
{
    const V4IElectronRepulsionDriver drv;

    const auto overlap =
        I4CIntegral(TwoCenterPair("a", 1, "b", 0), TwoCenterPair("c", 0, "d", 0), Operator("1"), 0, {});

    EXPECT_TRUE(drv.bra_hrr(overlap).empty());
}
