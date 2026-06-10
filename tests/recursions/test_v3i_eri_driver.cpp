// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v3i_eri_driver.hpp"
#include "t3c_defs.hpp"

namespace {

// Three-center electron-repulsion integral over the three-center types.
// Centers index as 0 = A (bra), 1 = C, 2 = D (ket pair).
I3CIntegral eri(int a, int c, int d)
{
    return I3CIntegral(OneCenter("a", a), TwoCenterPair("c", c, "d", d),
                       Operator("1/|r-r'|"), 0, {});
}

bool contains(const SI3CIntegrals& set, const I3CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V3IElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const V3IElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.is_electron_repulsion(eri(1, 0, 0)));

    // Wrong operator -> not an electron-repulsion integral.
    EXPECT_FALSE(drv.is_electron_repulsion(
        I3CIntegral(OneCenter("a", 1), TwoCenterPair("c", 0, "d", 0), Operator("1"), 0, {})));
}

TEST(V3IElectronRepulsionDriverTest, KetHrrLowersCRaisesD)
{
    const V3IElectronRepulsionDriver drv;

    const auto tints = drv.ket_hrr(eri(1, 1, 0));

    // (A1,C0,D0) plus the raised (A1,C0,D1) term.
    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, eri(1, 0, 0)));
    EXPECT_TRUE(contains(tints, eri(1, 0, 1)));
}

TEST(V3IElectronRepulsionDriverTest, BraVrrOnPShell)
{
    const V3IElectronRepulsionDriver drv;

    // A=1: only the order-raised (A0) term survives (no further lowering).
    const auto tints = drv.bra_vrr(eri(1, 0, 0));

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints,
        I3CIntegral(OneCenter("a", 0), TwoCenterPair("c", 0, "d", 0), Operator("1/|r-r'|"), 1, {})));
}

TEST(V3IElectronRepulsionDriverTest, BraVrrOnDShell)
{
    const V3IElectronRepulsionDriver drv;

    // A=2: order-raised (A1), lowered (A0), and order-raised (A0).
    const auto tints = drv.bra_vrr(eri(2, 0, 0));

    EXPECT_EQ(tints.size(), 3u);
}

TEST(V3IElectronRepulsionDriverTest, KetVrrLowersD)
{
    const V3IElectronRepulsionDriver drv;

    const auto tints = drv.ket_vrr(eri(0, 0, 1));

    // (D0) plus its order-raised partner.
    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, eri(0, 0, 0)));
}

TEST(V3IElectronRepulsionDriverTest, KetHrrRejectsNonEri)
{
    const V3IElectronRepulsionDriver drv;

    const auto nonEri =
        I3CIntegral(OneCenter("a", 1), TwoCenterPair("c", 1, "d", 0), Operator("1"), 0, {});

    EXPECT_TRUE(drv.ket_hrr(nonEri).empty());
    EXPECT_TRUE(drv.bra_vrr(nonEri).empty());
    EXPECT_TRUE(drv.ket_vrr(nonEri).empty());
}

TEST(V3IElectronRepulsionDriverTest, CreateVrrRecursionContainsSeed)
{
    const V3IElectronRepulsionDriver drv;

    const auto tints = drv.create_vrr_recursion(SI3CIntegrals({eri(1, 0, 0)}));

    EXPECT_TRUE(contains(tints, eri(1, 0, 0)));
    EXPECT_GE(tints.size(), 2u);
}
