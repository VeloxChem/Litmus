// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v4i_geom11_eri_driver.hpp"
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

TEST(V4IGeom11ElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const V4IGeom11ElectronRepulsionDriver drv;

    // Required prefix order is {1, 1, 0, 0}.
    EXPECT_TRUE(drv.is_electron_repulsion(geom(1, 0, 0, 0, 1, 1, 0, 0)));

    // Wrong prefix order.
    EXPECT_FALSE(drv.is_electron_repulsion(geom(1, 0, 0, 0, 1, 0, 0, 0)));

    // Wrong operator.
    EXPECT_FALSE(drv.is_electron_repulsion(geom(1, 0, 0, 0, 1, 1, 0, 0, "1")));
}

TEST(V4IGeom11ElectronRepulsionDriverTest, BraHrrProducesFourTerms)
{
    const V4IGeom11ElectronRepulsionDriver drv;

    // (P,S|S,S) with prefix {1,1,0,0}: lower A, lower prefix A, lower prefix B,
    // and raise B -> four distinct terms.
    const auto tints = drv.bra_hrr(geom(1, 0, 0, 0, 1, 1, 0, 0));

    EXPECT_EQ(tints.size(), 4u);

    // Lower-prefix-on-A term: prefix {0,1,0,0}.
    EXPECT_TRUE(tints.find(geom(0, 0, 0, 0, 0, 1, 0, 0)) != tints.cend());

    // Lower-prefix-on-B term: prefix {1,0,0,0}.
    EXPECT_TRUE(tints.find(geom(0, 0, 0, 0, 1, 0, 0, 0)) != tints.cend());

    // Raised-B term keeps the full prefix.
    EXPECT_TRUE(tints.find(geom(0, 1, 0, 0, 1, 1, 0, 0)) != tints.cend());
}

TEST(V4IGeom11ElectronRepulsionDriverTest, BraHrrRejectsNonEri)
{
    const V4IGeom11ElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.bra_hrr(geom(1, 0, 0, 0, 1, 0, 0, 0)).empty());
}

TEST(V4IGeom11ElectronRepulsionDriverTest, ApplyBraHrrRecursionExpands)
{
    const V4IGeom11ElectronRepulsionDriver drv;

    const auto tints = drv.apply_bra_hrr_recursion(geom(1, 0, 0, 0, 1, 1, 0, 0));

    EXPECT_FALSE(tints.empty());
    EXPECT_GE(tints.size(), 3u);
}

TEST(V4IGeom11ElectronRepulsionDriverTest, ApplyBraHrrRecursionOnScalarAIsEmpty)
{
    const V4IGeom11ElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.apply_bra_hrr_recursion(geom(0, 0, 0, 0, 1, 1, 0, 0)).empty());
}
