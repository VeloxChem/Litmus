// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t3c_vrr_eri_driver.hpp"
#include "t3c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);
const TensorComponent Dx(2, 0, 0);

// Three-center ERI term: one-center bra (A) and ket pair (C, D).
// Centers index as 0 = A, 1 = C, 2 = D.
R3CTerm eri_term(const TensorComponent& a, const TensorComponent& c, const TensorComponent& d)
{
    return R3CTerm(T3CIntegral(OneCenterComponent("a", a),
                               TwoCenterPairComponent({"c", "d"}, {c, d}),
                               OperatorComponent("1/|r-r'|")));
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

TEST(T3CVrrElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const T3CVrrElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.is_electron_repulsion(eri_term(Px, S, S)));

    const auto overlap = R3CTerm(T3CIntegral(OneCenterComponent("a", Px),
                                             TwoCenterPairComponent({"c", "d"}, {S, S}),
                                             OperatorComponent("1")));
    EXPECT_FALSE(drv.is_electron_repulsion(overlap));
}

// bra_vrr lowers center A (index 0); only the WA term survives on a P shell.
TEST(T3CVrrElectronRepulsionDriverTest, BraVrrOnPShell)
{
    const T3CVrrElectronRepulsionDriver drv;

    const auto rec = drv.bra_vrr(eri_term(Px, S, S), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"WA"}));
}

TEST(T3CVrrElectronRepulsionDriverTest, BraVrrOnDShell)
{
    const T3CVrrElectronRepulsionDriver drv;

    const auto rec = drv.bra_vrr(eri_term(Dx, S, S), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 3u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"WA", "1/b_e", "zeta/b_e^2"}));
}

// ket_vrr lowers center D (index 2).
TEST(T3CVrrElectronRepulsionDriverTest, KetVrrOnDShell)
{
    const T3CVrrElectronRepulsionDriver drv;

    const auto rec = drv.ket_vrr(eri_term(S, S, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"QD", "WQ"}));
}

TEST(T3CVrrElectronRepulsionDriverTest, BraVrrWrongAxisIsNullopt)
{
    const T3CVrrElectronRepulsionDriver drv;

    EXPECT_FALSE(drv.bra_vrr(eri_term(Px, S, S), 'y').has_value());
}

TEST(T3CVrrElectronRepulsionDriverTest, VrrRejectsNonEriTerm)
{
    const T3CVrrElectronRepulsionDriver drv;

    const auto overlap = R3CTerm(T3CIntegral(OneCenterComponent("a", Px),
                                             TwoCenterPairComponent({"c", "d"}, {S, S}),
                                             OperatorComponent("1")));

    EXPECT_FALSE(drv.bra_vrr(overlap, 'x').has_value());
}

TEST(T3CVrrElectronRepulsionDriverTest, CreateRecursionBuildsOneExpansionPerIntegral)
{
    const T3CVrrElectronRepulsionDriver drv;

    const VT3CIntegrals vints({T3CIntegral(OneCenterComponent("a", Px),
                                           TwoCenterPairComponent({"c", "d"}, {S, Px}),
                                           OperatorComponent("1/|r-r'|"))});

    EXPECT_EQ(drv.create_recursion(vints).expansions(), 1u);
}
