// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <vector>

#include "v3i_geom010_eri_driver.hpp"
#include "t3c_defs.hpp"

namespace {

// Geometric-derivative (0,1,0) three-center ERI integral: a first-order
// derivative prefix on center C plus scalar prefixes on A and D.
// Centers index as 0 = A (bra), 1 = C, 2 = D (ket pair).
I3CIntegral geom(int a, int c, int d)
{
    const VOperators prefixes({Operator("d/dR", Tensor(0)),
                               Operator("d/dR", Tensor(1)),
                               Operator("d/dR", Tensor(0))});

    return I3CIntegral(OneCenter("a", a), TwoCenterPair("c", c, "d", d),
                       Operator("1/|r-r'|"), 0, prefixes);
}

bool contains(const SI3CIntegrals& set, const I3CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V3IGeom010ElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const V3IGeom010ElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.is_electron_repulsion(geom(0, 1, 0)));
    EXPECT_EQ(geom(0, 1, 0).prefixes_order(), std::vector<int>({0, 1, 0}));

    // No prefixes -> wrong derivative order.
    EXPECT_FALSE(drv.is_electron_repulsion(
        I3CIntegral(OneCenter("a", 0), TwoCenterPair("c", 1, "d", 0), Operator("1/|r-r'|"), 0, {})));

    // Right derivative order but wrong operator.
    const VOperators prefixes({Operator("d/dR", Tensor(0)),
                               Operator("d/dR", Tensor(1)),
                               Operator("d/dR", Tensor(0))});
    EXPECT_FALSE(drv.is_electron_repulsion(
        I3CIntegral(OneCenter("a", 0), TwoCenterPair("c", 1, "d", 0), Operator("1"), 0, prefixes)));
}

TEST(V3IGeom010ElectronRepulsionDriverTest, KetHrrThreeTerms)
{
    const V3IGeom010ElectronRepulsionDriver drv;

    const auto tints = drv.ket_hrr(geom(0, 1, 0));

    // (1) lowered-C prefixed term, (2) prefix-reduced base term, (3) raised-D term.
    ASSERT_EQ(tints.size(), 3u);

    // The prefix-reduced second term is the plain ERI base with no prefixes.
    const auto base =
        I3CIntegral(OneCenter("a", 0), TwoCenterPair("c", 0, "d", 0), Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(contains(tints, base));

    // The lowered-C term keeps the derivative prefix.
    EXPECT_TRUE(contains(tints, geom(0, 0, 0)));
    // The raised-D term keeps the derivative prefix.
    EXPECT_TRUE(contains(tints, geom(0, 0, 1)));
}

TEST(V3IGeom010ElectronRepulsionDriverTest, KetAuxHrrOnScalarC)
{
    const V3IGeom010ElectronRepulsionDriver drv;

    // C==0 with derivative order (0,1,0): base term plus raised-D base term.
    const auto tints = drv.ket_aux_hrr(geom(0, 0, 0));

    ASSERT_EQ(tints.size(), 2u);

    const auto base =
        I3CIntegral(OneCenter("a", 0), TwoCenterPair("c", 0, "d", 0), Operator("1/|r-r'|"), 0, {});
    const auto raised =
        I3CIntegral(OneCenter("a", 0), TwoCenterPair("c", 0, "d", 1), Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(contains(tints, base));
    EXPECT_TRUE(contains(tints, raised));
}

TEST(V3IGeom010ElectronRepulsionDriverTest, KetAuxHrrGuardsOnPositiveC)
{
    const V3IGeom010ElectronRepulsionDriver drv;

    // C > 0 short-circuits to the empty set.
    EXPECT_TRUE(drv.ket_aux_hrr(geom(0, 1, 0)).empty());
}

TEST(V3IGeom010ElectronRepulsionDriverTest, KetHrrRejectsNonEri)
{
    const V3IGeom010ElectronRepulsionDriver drv;

    // No prefixes -> wrong derivative order -> empty.
    const auto plain =
        I3CIntegral(OneCenter("a", 0), TwoCenterPair("c", 1, "d", 0), Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(drv.ket_hrr(plain).empty());
    EXPECT_TRUE(drv.ket_aux_hrr(plain).empty());
}
