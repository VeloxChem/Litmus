// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <vector>

#include "t3c_geom_100_eri_driver.hpp"
#include "t3c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

// Geometric-derivative (1,0,0) three-center ERI term: a first-order derivative
// prefix on center A plus scalar prefixes on C and D.
R3CTerm geom_term(const TensorComponent& a, const TensorComponent& pA)
{
    const VOperatorComponents prefixes({OperatorComponent("d/dA", pA, "bra", 0),
                                        OperatorComponent("d/dC", S, "ket", 1),
                                        OperatorComponent("d/dD", S, "ket", 2)});

    return R3CTerm(T3CIntegral(OneCenterComponent("a", a),
                               TwoCenterPairComponent({"c", "d"}, {S, S}),
                               OperatorComponent("1/|r-r'|"), 0, prefixes));
}

}  // namespace

TEST(T3CGeom100ElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const T3CGeom100ElectronRepulsionDriver drv;

    // Derivative order (1,0,0) over the ERI operator.
    EXPECT_TRUE(drv.is_electron_repulsion(geom_term(S, Px)));
    EXPECT_EQ(geom_term(S, Px).prefixes_order(), std::vector<int>({1, 0, 0}));

    // Scalar prefixes -> wrong derivative order.
    EXPECT_FALSE(drv.is_electron_repulsion(geom_term(S, S)));
}

TEST(T3CGeom100ElectronRepulsionDriverTest, BraHrrConvertsPrefixToMomentum)
{
    const T3CGeom100ElectronRepulsionDriver drv;

    // Scalar bra: the derivative becomes a single raised, prefix-free integral.
    const auto rec = drv.bra_hrr(geom_term(S, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
    EXPECT_TRUE((*rec)[0].integral().prefixes().empty());  // prefixes cleared
}

TEST(T3CGeom100ElectronRepulsionDriverTest, BraHrrOnPShellGivesRaiseAndLower)
{
    const T3CGeom100ElectronRepulsionDriver drv;

    // P bra: the derivative identity gives both a raise and a lower term.
    const auto rec = drv.bra_hrr(geom_term(Px, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
}

TEST(T3CGeom100ElectronRepulsionDriverTest, BraHrrWrongAxisIsNullopt)
{
    const T3CGeom100ElectronRepulsionDriver drv;

    // The derivative prefix has no y component to lower.
    EXPECT_FALSE(drv.bra_hrr(geom_term(S, Px), 'y').has_value());
}

TEST(T3CGeom100ElectronRepulsionDriverTest, BraHrrRejectsNonEriTerm)
{
    const T3CGeom100ElectronRepulsionDriver drv;

    EXPECT_FALSE(drv.bra_hrr(geom_term(S, S), 'x').has_value());  // wrong derivative order
}
