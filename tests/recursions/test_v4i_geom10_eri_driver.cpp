// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v4i_geom10_eri_driver.hpp"
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

TEST(V4IGeom10ElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const V4IGeom10ElectronRepulsionDriver drv;

    // Accepted prefix orders: {1,0,0,0}, {0,0,1,0}, {1,0,1,0}.
    EXPECT_TRUE(drv.is_electron_repulsion(geom(1, 0, 0, 0, 1, 0, 0, 0)));
    EXPECT_TRUE(drv.is_electron_repulsion(geom(0, 0, 1, 0, 0, 0, 1, 0)));
    EXPECT_TRUE(drv.is_electron_repulsion(geom(1, 0, 1, 0, 1, 0, 1, 0)));

    // Rejected prefix order.
    EXPECT_FALSE(drv.is_electron_repulsion(geom(0, 1, 0, 0, 0, 1, 0, 0)));

    // Wrong operator.
    EXPECT_FALSE(drv.is_electron_repulsion(geom(1, 0, 0, 0, 1, 0, 0, 0, "1")));
}

TEST(V4IGeom10ElectronRepulsionDriverTest, BraHrrProducesThreeTerms)
{
    const V4IGeom10ElectronRepulsionDriver drv;

    // (P,S|S,S) with prefix {1,0,0,0}: lower A, lower-prefix (-> base), raise B.
    const auto tints = drv.bra_hrr(geom(1, 0, 0, 0, 1, 0, 0, 0));

    EXPECT_EQ(tints.size(), 3u);

    // The lowered-prefix term collapses to scalar {0,0,0,0} -> base integral.
    const auto base = I4CIntegral(TwoCenterPair("a", 0, "b", 0),
                                  TwoCenterPair("c", 0, "d", 0),
                                  Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(tints.find(base) != tints.cend());

    // Raised-B term keeps the prefix.
    EXPECT_TRUE(tints.find(geom(0, 1, 0, 0, 1, 0, 0, 0)) != tints.cend());
}

TEST(V4IGeom10ElectronRepulsionDriverTest, KetHrrProducesThreeTerms)
{
    const V4IGeom10ElectronRepulsionDriver drv;

    // (S,S|P,S) with prefix {0,0,1,0}: lower C, lower-prefix (-> base), raise D.
    const auto tints = drv.ket_hrr(geom(0, 0, 1, 0, 0, 0, 1, 0));

    EXPECT_EQ(tints.size(), 3u);

    const auto base = I4CIntegral(TwoCenterPair("a", 0, "b", 0),
                                  TwoCenterPair("c", 0, "d", 0),
                                  Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(tints.find(base) != tints.cend());

    EXPECT_TRUE(tints.find(geom(0, 0, 0, 1, 0, 0, 1, 0)) != tints.cend());
}

TEST(V4IGeom10ElectronRepulsionDriverTest, BraAuxHrrOnScalarA)
{
    const V4IGeom10ElectronRepulsionDriver drv;

    // A momentum zero, prefix {1,0,0,0}: prefix lowered and B raised.
    const auto tints = drv.bra_aux_hrr(geom(0, 0, 0, 0, 1, 0, 0, 0));

    EXPECT_EQ(tints.size(), 2u);
}

TEST(V4IGeom10ElectronRepulsionDriverTest, BraAuxHrrSkipsNonZeroA)
{
    const V4IGeom10ElectronRepulsionDriver drv;

    // A momentum is non-zero -> bra_aux_hrr returns empty.
    EXPECT_TRUE(drv.bra_aux_hrr(geom(1, 0, 0, 0, 1, 0, 0, 0)).empty());
}

TEST(V4IGeom10ElectronRepulsionDriverTest, KetAuxHrrOnScalarC)
{
    const V4IGeom10ElectronRepulsionDriver drv;

    // C momentum zero, prefix {0,0,1,0}: two base terms.
    const auto tints = drv.ket_aux_hrr(geom(0, 0, 0, 0, 0, 0, 1, 0));

    EXPECT_EQ(tints.size(), 2u);
}

TEST(V4IGeom10ElectronRepulsionDriverTest, BraHrrRejectsNonEri)
{
    const V4IGeom10ElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.bra_hrr(geom(1, 1, 0, 0, 1, 1, 0, 0)).empty());
}

TEST(V4IGeom10ElectronRepulsionDriverTest, ApplyBraHrrRecursionExpands)
{
    const V4IGeom10ElectronRepulsionDriver drv;

    const auto tints = drv.apply_bra_hrr_recursion(geom(1, 0, 0, 0, 1, 0, 0, 0));

    EXPECT_FALSE(tints.empty());

    const auto base = I4CIntegral(TwoCenterPair("a", 0, "b", 0),
                                  TwoCenterPair("c", 0, "d", 0),
                                  Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(tints.find(base) != tints.cend());
}

TEST(V4IGeom10ElectronRepulsionDriverTest, ApplyKetHrrRecursionExpands)
{
    const V4IGeom10ElectronRepulsionDriver drv;

    const auto tints = drv.apply_ket_hrr_recursion(geom(0, 0, 1, 0, 0, 0, 1, 0));

    EXPECT_FALSE(tints.empty());

    const auto base = I4CIntegral(TwoCenterPair("a", 0, "b", 0),
                                  TwoCenterPair("c", 0, "d", 0),
                                  Operator("1/|r-r'|"), 0, {});
    EXPECT_TRUE(tints.find(base) != tints.cend());
}
