// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>
#include <vector>

#include "t4c_geom_10_hrr_eri_driver.hpp"
#include "t4c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

// Geometric-derivative four-center ERI term with a derivative prefix of the
// given shape on center A (the (1,0,0,0) family) and scalar prefixes elsewhere.
R4CTerm geom_term(const TensorComponent& a, const TensorComponent& b, const TensorComponent& pA)
{
    const VOperatorComponents prefixes({OperatorComponent("d/dA", pA, "bra", 0),
                                        OperatorComponent("d/dB", S, "bra", 1),
                                        OperatorComponent("d/dC", S, "ket", 2),
                                        OperatorComponent("d/dD", S, "ket", 3)});

    return R4CTerm(T4CIntegral(TwoCenterPairComponent({"a", "b"}, {a, b}),
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

TEST(T4CGeom10HrrElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const T4CGeom10HrrElectronRepulsionDriver drv;

    // Accepted derivative orders include (1,0,0,0).
    EXPECT_TRUE(drv.is_electron_repulsion(geom_term(Px, S, Px)));
    EXPECT_EQ(geom_term(Px, S, Px).prefixes_order(), std::vector<int>({1, 0, 0, 0}));

    // No derivative prefix -> not in the accepted set.
    EXPECT_FALSE(drv.is_electron_repulsion(geom_term(Px, S, S)));
}

TEST(T4CGeom10HrrElectronRepulsionDriverTest, BraHrrTransfersAndStepsPrefix)
{
    const T4CGeom10HrrElectronRepulsionDriver drv;

    // (A=P, B=S | S, S) with derivative on A: BA transfer term, the prefix
    // down-step term, and the B raise term.
    const auto rec = drv.bra_hrr(geom_term(Px, S, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 3u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"BA"}));
}

TEST(T4CGeom10HrrElectronRepulsionDriverTest, BraHrrRequiresBraMomentum)
{
    const T4CGeom10HrrElectronRepulsionDriver drv;

    // bra_hrr lowers center A first; a scalar A has nothing to lower.
    EXPECT_FALSE(drv.bra_hrr(geom_term(S, S, Px), 'x').has_value());
}

TEST(T4CGeom10HrrElectronRepulsionDriverTest, BraHrrRejectsNonEriTerm)
{
    const T4CGeom10HrrElectronRepulsionDriver drv;

    EXPECT_FALSE(drv.bra_hrr(geom_term(Px, S, S), 'x').has_value());  // no derivative prefix
}
