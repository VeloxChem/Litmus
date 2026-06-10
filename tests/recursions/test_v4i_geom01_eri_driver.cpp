// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v4i_geom01_eri_driver.hpp"
#include "t4c_defs.hpp"

namespace {

// Four-center ERI with geometrical-derivative prefixes of orders {pa,pb,pc,pd}.
I4CIntegral geom(int a, int b, int c, int d, int pa, int pb, int pc, int pd,
                 const std::string& integrand = "1/|r-r'|")
{
    const VOperators prefixes({Operator("d/dA", Tensor(pa)),
                               Operator("d/dB", Tensor(pb)),
                               Operator("d/dC", Tensor(pc)),
                               Operator("d/dD", Tensor(pd))});

    return I4CIntegral(TwoCenterPair("a", a, "b", b),
                       TwoCenterPair("c", c, "d", d),
                       Operator(integrand), 0, prefixes);
}

}  // namespace

TEST(V4IGeom01ElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const V4IGeom01ElectronRepulsionDriver drv;

    // Required prefix order is {0, 1, 0, 0}.
    EXPECT_TRUE(drv.is_electron_repulsion(geom(0, 1, 0, 0, 0, 1, 0, 0)));

    // Wrong prefix order.
    EXPECT_FALSE(drv.is_electron_repulsion(geom(1, 0, 0, 0, 1, 0, 0, 0)));

    // Wrong operator.
    EXPECT_FALSE(drv.is_electron_repulsion(geom(0, 1, 0, 0, 0, 1, 0, 0, "1")));

    // No prefixes at all.
    EXPECT_FALSE(drv.is_electron_repulsion(
        I4CIntegral(TwoCenterPair("a", 0, "b", 1),
                    TwoCenterPair("c", 0, "d", 0),
                    Operator("1/|r-r'|"), 0, {})));
}

TEST(V4IGeom01ElectronRepulsionDriverTest, BraHrrProducesThreeTerms)
{
    const V4IGeom01ElectronRepulsionDriver drv;

    // (P,S|S,S) with prefix {0,1,0,0}: lower A, base term, and raise B.
    const auto tints = drv.bra_hrr(geom(1, 0, 0, 0, 0, 1, 0, 0));

    EXPECT_EQ(tints.size(), 3u);

    // Base (prefix-less) recursion term.
    const auto base = I4CIntegral(TwoCenterPair("a", 0, "b", 0),
                                  TwoCenterPair("c", 0, "d", 0),
                                  Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(tints.find(base) != tints.cend());

    // Raised-B term keeps the prefix.
    const auto raised = geom(0, 1, 0, 0, 0, 1, 0, 0);
    EXPECT_TRUE(tints.find(raised) != tints.cend());
}

TEST(V4IGeom01ElectronRepulsionDriverTest, BraHrrRejectsNonEri)
{
    const V4IGeom01ElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.bra_hrr(geom(1, 0, 0, 0, 1, 0, 0, 0)).empty());
}

TEST(V4IGeom01ElectronRepulsionDriverTest, ApplyBraHrrRecursionExpands)
{
    const V4IGeom01ElectronRepulsionDriver drv;

    const auto tints = drv.apply_bra_hrr_recursion(geom(1, 0, 0, 0, 0, 1, 0, 0));

    EXPECT_FALSE(tints.empty());

    // The fully reduced base term is part of the expansion.
    const auto base = I4CIntegral(TwoCenterPair("a", 0, "b", 0),
                                  TwoCenterPair("c", 0, "d", 0),
                                  Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(tints.find(base) != tints.cend());
}

TEST(V4IGeom01ElectronRepulsionDriverTest, ApplyBraHrrRecursionOnScalarAIsEmpty)
{
    const V4IGeom01ElectronRepulsionDriver drv;

    // A momentum is zero -> no bra recursion.
    EXPECT_TRUE(drv.apply_bra_hrr_recursion(geom(0, 1, 0, 0, 0, 1, 0, 0)).empty());
}
