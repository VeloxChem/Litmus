// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v4i_center_driver.hpp"
#include "t4c_defs.hpp"

namespace {

// Four-center ERI integral with geometrical-derivative prefixes of the given
// orders on centers A, B, C, D.
I4CIntegral geom(int a, int b, int c, int d, int pa, int pb, int pc, int pd)
{
    const VOperators prefixes({Operator("d/dA", Tensor(pa)),
                               Operator("d/dB", Tensor(pb)),
                               Operator("d/dC", Tensor(pc)),
                               Operator("d/dD", Tensor(pd))});

    return I4CIntegral(TwoCenterPair("a", a, "b", b),
                       TwoCenterPair("c", c, "d", d),
                       Operator("1/|r-r'|"), 0, prefixes);
}

}  // namespace

TEST(V4ICenterDriverTest, IsAuxiliary)
{
    const V4ICenterDriver drv;

    // Scalar (order-0) prefix on A is auxiliary; order-1 prefix is not.
    EXPECT_TRUE(drv.is_auxiliary(geom(1, 0, 0, 0, 0, 0, 0, 0), 0));
    EXPECT_FALSE(drv.is_auxiliary(geom(1, 0, 0, 0, 1, 0, 0, 0), 0));

    // Order-1 prefix on B is not auxiliary on index 1, but A (index 0) is.
    EXPECT_TRUE(drv.is_auxiliary(geom(0, 1, 0, 0, 0, 1, 0, 0), 0));
    EXPECT_FALSE(drv.is_auxiliary(geom(0, 1, 0, 0, 0, 1, 0, 0), 1));
}

TEST(V4ICenterDriverTest, BraKetVrrOnAuxiliaryIsEmpty)
{
    const V4ICenterDriver drv;

    // Auxiliary (scalar) prefix on A -> no recursion terms.
    EXPECT_TRUE(drv.bra_ket_vrr(geom(1, 0, 0, 0, 0, 0, 0, 0), 0).empty());
}

TEST(V4ICenterDriverTest, BraKetVrrLowersPrefixAndShiftsCenter)
{
    const V4ICenterDriver drv;

    // Order-1 prefix on A, A momentum 1: prefix reduces to scalar (cleared),
    // then center A is shifted up (-> momentum 2) and down (-> momentum 0).
    const auto tints = drv.bra_ket_vrr(geom(1, 0, 0, 0, 1, 0, 0, 0), 0);

    EXPECT_EQ(tints.size(), 2u);

    const auto up   = I4CIntegral(TwoCenterPair("a", 2, "b", 0),
                                  TwoCenterPair("c", 0, "d", 0),
                                  Operator("1/|r-r'|"), 0, {});
    const auto down = I4CIntegral(TwoCenterPair("a", 0, "b", 0),
                                  TwoCenterPair("c", 0, "d", 0),
                                  Operator("1/|r-r'|"), 0, {});

    EXPECT_TRUE(tints.find(up) != tints.cend());
    EXPECT_TRUE(tints.find(down) != tints.cend());
}

TEST(V4ICenterDriverTest, ApplyBraKetVrrIncludesSeedAndExpansion)
{
    const V4ICenterDriver drv;

    const auto seed  = geom(1, 0, 0, 0, 1, 0, 0, 0);
    const auto tints = drv.apply_bra_ket_vrr(seed);

    // Seed integral plus its two A-center expansion terms.
    EXPECT_EQ(tints.size(), 3u);
    EXPECT_TRUE(tints.find(seed) != tints.cend());

    const auto down = I4CIntegral(TwoCenterPair("a", 0, "b", 0),
                                  TwoCenterPair("c", 0, "d", 0),
                                  Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(tints.find(down) != tints.cend());
}

TEST(V4ICenterDriverTest, ApplyBraKetVrrWithoutPrefixesReturnsSelf)
{
    const V4ICenterDriver drv;

    const auto plain = I4CIntegral(TwoCenterPair("a", 1, "b", 0),
                                   TwoCenterPair("c", 0, "d", 0),
                                   Operator("1/|r-r'|"), 0, {});

    const auto tints = drv.apply_bra_ket_vrr(plain);

    EXPECT_EQ(tints.size(), 1u);
    EXPECT_TRUE(tints.find(plain) != tints.cend());
}
