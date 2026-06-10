// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <vector>

#include "v3i_geom100_eri_driver.hpp"
#include "t3c_defs.hpp"

namespace {

// Geometric-derivative (1,0,0) three-center ERI integral: a first-order
// derivative prefix on center A plus scalar prefixes on C and D.
// Centers index as 0 = A (bra), 1 = C, 2 = D (ket pair).
I3CIntegral geom(int a, int c, int d)
{
    const VOperators prefixes({Operator("d/dR", Tensor(1)),
                               Operator("d/dR", Tensor(0)),
                               Operator("d/dR", Tensor(0))});

    return I3CIntegral(OneCenter("a", a), TwoCenterPair("c", c, "d", d),
                       Operator("1/|r-r'|"), 0, prefixes);
}

bool contains(const SI3CIntegrals& set, const I3CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V3IGeom100ElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const V3IGeom100ElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.is_electron_repulsion(geom(0, 0, 0)));
    EXPECT_EQ(geom(0, 0, 0).prefixes_order(), std::vector<int>({1, 0, 0}));

    // No prefixes -> wrong derivative order.
    EXPECT_FALSE(drv.is_electron_repulsion(
        I3CIntegral(OneCenter("a", 0), TwoCenterPair("c", 0, "d", 0), Operator("1/|r-r'|"), 0, {})));

    // Right derivative order but wrong operator.
    const VOperators prefixes({Operator("d/dR", Tensor(1)),
                               Operator("d/dR", Tensor(0)),
                               Operator("d/dR", Tensor(0))});
    EXPECT_FALSE(drv.is_electron_repulsion(
        I3CIntegral(OneCenter("a", 0), TwoCenterPair("c", 0, "d", 0), Operator("1"), 0, prefixes)));
}

TEST(V3IGeom100ElectronRepulsionDriverTest, BraHrrOnScalarBra)
{
    const V3IGeom100ElectronRepulsionDriver drv;

    // A==0: only the raised, prefix-free term survives.
    const auto tints = drv.bra_hrr(geom(0, 0, 0));

    ASSERT_EQ(tints.size(), 1u);

    const auto raised =
        I3CIntegral(OneCenter("a", 1), TwoCenterPair("c", 0, "d", 0), Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(contains(tints, raised));
    // Resulting terms are prefix-free base integrals.
    EXPECT_TRUE(tints.cbegin()->prefixes().empty());
}

TEST(V3IGeom100ElectronRepulsionDriverTest, BraHrrOnPShellGivesRaiseAndLower)
{
    const V3IGeom100ElectronRepulsionDriver drv;

    // A==1: a raised (A2) and a lowered (A0) prefix-free term.
    const auto tints = drv.bra_hrr(geom(1, 0, 0));

    ASSERT_EQ(tints.size(), 2u);

    const auto raised =
        I3CIntegral(OneCenter("a", 2), TwoCenterPair("c", 0, "d", 0), Operator("1/|r-r'|"), 0, {});
    const auto lowered =
        I3CIntegral(OneCenter("a", 0), TwoCenterPair("c", 0, "d", 0), Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(contains(tints, raised));
    EXPECT_TRUE(contains(tints, lowered));
}

TEST(V3IGeom100ElectronRepulsionDriverTest, BraHrrRejectsNonEri)
{
    const V3IGeom100ElectronRepulsionDriver drv;

    // No prefixes -> wrong derivative order -> empty.
    const auto plain =
        I3CIntegral(OneCenter("a", 1), TwoCenterPair("c", 0, "d", 0), Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(drv.bra_hrr(plain).empty());
}

TEST(V3IGeom100ElectronRepulsionDriverTest, ApplyBraHrrRecursionExpandsPShell)
{
    const V3IGeom100ElectronRepulsionDriver drv;

    // The P-shell derivative expands into prefix-free raise/lower base terms.
    const auto tints = drv.apply_bra_hrr_recursion(geom(1, 0, 0));

    EXPECT_GE(tints.size(), 2u);

    const auto raised =
        I3CIntegral(OneCenter("a", 2), TwoCenterPair("c", 0, "d", 0), Operator("1/|r-r'|"), 0, {});
    const auto lowered =
        I3CIntegral(OneCenter("a", 0), TwoCenterPair("c", 0, "d", 0), Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(contains(tints, raised));
    EXPECT_TRUE(contains(tints, lowered));
}

TEST(V3IGeom100ElectronRepulsionDriverTest, ApplyBraHrrRecursionGuardsOnScalarBra)
{
    const V3IGeom100ElectronRepulsionDriver drv;

    // A==0 short-circuits the driver loop to the empty set.
    EXPECT_TRUE(drv.apply_bra_hrr_recursion(geom(0, 0, 0)).empty());
}
