// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>
#include <vector>

#include "t4c_geom_11_hrr_eri_driver.hpp"
#include "t4c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

// Geometric-derivative (1,1,0,0) four-center ERI term: first-order derivative
// prefixes on centers A and B, scalar prefixes on C and D.
R4CTerm geom_term(const TensorComponent& pA, const TensorComponent& pB)
{
    const VOperatorComponents prefixes({OperatorComponent("d/dA", pA, "bra", 0),
                                        OperatorComponent("d/dB", pB, "bra", 1),
                                        OperatorComponent("d/dC", S, "ket", 2),
                                        OperatorComponent("d/dD", S, "ket", 3)});

    return R4CTerm(T4CIntegral(TwoCenterPairComponent({"a", "b"}, {S, S}),
                               TwoCenterPairComponent({"c", "d"}, {S, S}),
                               OperatorComponent("1/|r-r'|"), 0, prefixes));
}

std::set<std::string> factor_names(const R4CDist& dist)
{
    std::set<std::string> names;
    for (const auto& factor : dist.unique_factors())
    {
        names.insert(factor.name());
    }
    return names;
}

}  // namespace

TEST(T4CGeom11HrrElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const T4CGeom11HrrElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.is_electron_repulsion(geom_term(Px, Px)));
    EXPECT_EQ(geom_term(Px, Px).prefixes_order(), std::vector<int>({1, 1, 0, 0}));

    // Only one derivative -> wrong order.
    EXPECT_FALSE(drv.is_electron_repulsion(geom_term(Px, S)));
}

TEST(T4CGeom11HrrElectronRepulsionDriverTest, BraAuxHrrTransfersAndStepsBothPrefixes)
{
    const T4CGeom11HrrElectronRepulsionDriver drv;

    // BA transfer term, the B raise term, and the second-prefix down-step term.
    const auto rec = drv.bra_aux_hrr(geom_term(Px, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 3u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"BA"}));
}

TEST(T4CGeom11HrrElectronRepulsionDriverTest, BraAuxHrrWrongAxisIsNullopt)
{
    const T4CGeom11HrrElectronRepulsionDriver drv;

    // The A derivative prefix has no y component to lower.
    EXPECT_FALSE(drv.bra_aux_hrr(geom_term(Px, Px), 'y').has_value());
}
