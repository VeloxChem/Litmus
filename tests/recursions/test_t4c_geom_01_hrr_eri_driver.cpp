// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>
#include <vector>

#include "t4c_geom_01_hrr_eri_driver.hpp"
#include "t4c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

// Four-center geometric-derivative ERI term of the (0,1,0,0) family: a
// derivative prefix of shape pB on center B (prefix index 1), scalar prefixes
// elsewhere. The bra pair holds (A, B), the ket pair (C, D).
R4CTerm geom_term(const TensorComponent& a, const TensorComponent& b,
                  const TensorComponent& pB)
{
    const VOperatorComponents prefixes({OperatorComponent("d/dA", S, "bra", 0),
                                        OperatorComponent("d/dB", pB, "bra", 1),
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

TEST(T4CGeom01HrrElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const T4CGeom01HrrElectronRepulsionDriver drv;

    // Accepted derivative order is (0,1,0,0) with the ERI operator.
    EXPECT_TRUE(drv.is_electron_repulsion(geom_term(Px, S, Px)));
    EXPECT_EQ(geom_term(Px, S, Px).prefixes_order(), std::vector<int>({0, 1, 0, 0}));

    // No derivative prefix on B -> wrong order.
    EXPECT_FALSE(drv.is_electron_repulsion(geom_term(Px, S, S)));

    // Right derivative order but wrong operator.
    const VOperatorComponents prefixes({OperatorComponent("d/dA", S, "bra", 0),
                                        OperatorComponent("d/dB", Px, "bra", 1),
                                        OperatorComponent("d/dC", S, "ket", 2),
                                        OperatorComponent("d/dD", S, "ket", 3)});
    EXPECT_FALSE(drv.is_electron_repulsion(
        R4CTerm(T4CIntegral(TwoCenterPairComponent({"a", "b"}, {Px, S}),
                            TwoCenterPairComponent({"c", "d"}, {S, S}),
                            OperatorComponent("1"), 0, prefixes))));
}

TEST(T4CGeom01HrrElectronRepulsionDriverTest, BraHrrTransfersAndStepsPrefix)
{
    const T4CGeom01HrrElectronRepulsionDriver drv;

    // (A=P, B=S | S, S) with derivative on B: BA transfer term, the prefix
    // down-step term, and the B raise term.
    const auto rec = drv.bra_hrr(geom_term(Px, S, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 3u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"BA"}));
}

TEST(T4CGeom01HrrElectronRepulsionDriverTest, BraHrrRequiresBraMomentum)
{
    const T4CGeom01HrrElectronRepulsionDriver drv;

    // bra_hrr lowers center A first; a scalar A has nothing to lower.
    EXPECT_FALSE(drv.bra_hrr(geom_term(S, S, Px), 'x').has_value());
}

TEST(T4CGeom01HrrElectronRepulsionDriverTest, BraHrrRejectsNonEriTerm)
{
    const T4CGeom01HrrElectronRepulsionDriver drv;

    // No derivative prefix on B -> not in the accepted set.
    EXPECT_FALSE(drv.bra_hrr(geom_term(Px, S, S), 'x').has_value());
}

TEST(T4CGeom01HrrElectronRepulsionDriverTest, BraAuxHrrStepsPrefixOnScalarB)
{
    const T4CGeom01HrrElectronRepulsionDriver drv;

    // Lowering the B prefix; B center scalar -> only the raise-B term survives.
    const auto rec = drv.bra_aux_hrr(geom_term(S, S, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
}

TEST(T4CGeom01HrrElectronRepulsionDriverTest, BraAuxHrrOnPCenterB)
{
    const T4CGeom01HrrElectronRepulsionDriver drv;

    // B center = P: raise-B term plus the (scaled) lower-B term.
    const auto rec = drv.bra_aux_hrr(geom_term(S, Px, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
}
