// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>
#include <vector>

#include "t4c_geom_20_hrr_eri_driver.hpp"
#include "t4c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);
const TensorComponent Dx(2, 0, 0);

// Geometric-derivative (2,0,0,0) four-center ERI term: a second-order
// derivative prefix on center A, scalar prefixes elsewhere.
R4CTerm geom_term(const TensorComponent& a, const TensorComponent& pA)
{
    const VOperatorComponents prefixes({OperatorComponent("d/dA", pA, "bra", 0),
                                        OperatorComponent("d/dB", S, "bra", 1),
                                        OperatorComponent("d/dC", S, "ket", 2),
                                        OperatorComponent("d/dD", S, "ket", 3)});

    return R4CTerm(T4CIntegral(TwoCenterPairComponent({"a", "b"}, {a, S}),
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

TEST(T4CGeom20HrrElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const T4CGeom20HrrElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.is_electron_repulsion(geom_term(Px, Dx)));
    EXPECT_EQ(geom_term(Px, Dx).prefixes_order(), std::vector<int>({2, 0, 0, 0}));

    // First-order derivative -> wrong order for the (2,0,0,0) driver.
    EXPECT_FALSE(drv.is_electron_repulsion(geom_term(Px, Px)));
}

TEST(T4CGeom20HrrElectronRepulsionDriverTest, BraHrrTransfersAndStepsPrefix)
{
    const T4CGeom20HrrElectronRepulsionDriver drv;

    // (A=P) with second-order derivative: BA transfer, prefix down-step, B raise.
    const auto rec = drv.bra_hrr(geom_term(Px, Dx), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 3u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"BA"}));
}

TEST(T4CGeom20HrrElectronRepulsionDriverTest, BraHrrRequiresBraMomentum)
{
    const T4CGeom20HrrElectronRepulsionDriver drv;

    EXPECT_FALSE(drv.bra_hrr(geom_term(S, Dx), 'x').has_value());
}

TEST(T4CGeom20HrrElectronRepulsionDriverTest, BraHrrRejectsNonEriTerm)
{
    const T4CGeom20HrrElectronRepulsionDriver drv;

    EXPECT_FALSE(drv.bra_hrr(geom_term(Px, Px), 'x').has_value());  // (1,0,0,0) order
}
