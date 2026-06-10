// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>
#include <vector>

#include "t3c_geom_010_eri_driver.hpp"
#include "t3c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

// Geometric-derivative (0,1,0) three-center ERI term: a first-order derivative
// prefix on ket center C, scalar prefixes on A and D.
R3CTerm geom_term(const TensorComponent& pC)
{
    const VOperatorComponents prefixes({OperatorComponent("d/dA", S, "bra", 0),
                                        OperatorComponent("d/dC", pC, "ket", 1),
                                        OperatorComponent("d/dD", S, "ket", 2)});

    return R3CTerm(T3CIntegral(OneCenterComponent("a", S),
                               TwoCenterPairComponent({"c", "d"}, {S, S}),
                               OperatorComponent("1/|r-r'|"), 0, prefixes));
}

std::set<std::string> factor_names(const R3CDist& dist)
{
    std::set<std::string> names;
    for (const auto& factor : dist.unique_factors())
    {
        names.insert(factor.name());
    }
    return names;
}

}  // namespace

TEST(T3CGeom010ElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const T3CGeom010ElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.is_electron_repulsion(geom_term(Px)));
    EXPECT_EQ(geom_term(Px).prefixes_order(), std::vector<int>({0, 1, 0}));

    EXPECT_FALSE(drv.is_electron_repulsion(geom_term(S)));  // scalar prefix -> wrong order
}

TEST(T3CGeom010ElectronRepulsionDriverTest, KetAuxHrrConvertsPrefixToMomentum)
{
    const T3CGeom010ElectronRepulsionDriver drv;

    // The C derivative becomes the DC transfer term plus the D raise term.
    const auto rec = drv.ket_aux_hrr(geom_term(Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"DC"}));
    EXPECT_TRUE((*rec)[0].integral().prefixes().empty());  // prefixes cleared
}

TEST(T3CGeom010ElectronRepulsionDriverTest, KetAuxHrrWrongAxisIsNullopt)
{
    const T3CGeom010ElectronRepulsionDriver drv;

    EXPECT_FALSE(drv.ket_aux_hrr(geom_term(Px), 'y').has_value());
}

TEST(T3CGeom010ElectronRepulsionDriverTest, KetAuxHrrRejectsNonEriTerm)
{
    const T3CGeom010ElectronRepulsionDriver drv;

    EXPECT_FALSE(drv.ket_aux_hrr(geom_term(S), 'x').has_value());  // wrong derivative order
}
