// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t4c_vrr_eri_driver.hpp"
#include "t4c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);
const TensorComponent Dx(2, 0, 0);

// Four-center ERI term over bra pair (A, B) and ket pair (C, D).
R4CTerm eri_term(const TensorComponent& a, const TensorComponent& b,
                 const TensorComponent& c, const TensorComponent& d)
{
    return R4CTerm(T4CIntegral(TwoCenterPairComponent({"a", "b"}, {a, b}),
                               TwoCenterPairComponent({"c", "d"}, {c, d}),
                               OperatorComponent("1/|r-r'|")));
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

TEST(T4CVrrElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const T4CVrrElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.is_electron_repulsion(eri_term(S, Px, S, S)));

    const auto overlap = R4CTerm(T4CIntegral(TwoCenterPairComponent({"a", "b"}, {S, Px}),
                                             TwoCenterPairComponent({"c", "d"}, {S, S}),
                                             OperatorComponent("1")));
    EXPECT_FALSE(drv.is_electron_repulsion(overlap));
}

// bra_vrr lowers center B (index 1), not A.
TEST(T4CVrrElectronRepulsionDriverTest, BraVrrOnBShellPShell)
{
    const T4CVrrElectronRepulsionDriver drv;

    const auto rec = drv.bra_vrr(eri_term(S, Px, S, S), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PB", "WP"}));
}

TEST(T4CVrrElectronRepulsionDriverTest, BraVrrOnBShellDShell)
{
    const T4CVrrElectronRepulsionDriver drv;

    const auto rec = drv.bra_vrr(eri_term(S, Dx, S, S), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 4u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PB", "WP", "1/eta", "rho/eta"}));
}

TEST(T4CVrrElectronRepulsionDriverTest, BraVrrRequiresBMomentum)
{
    const T4CVrrElectronRepulsionDriver drv;

    // A has momentum but B is scalar: bra_vrr lowers B, so there is nothing to do.
    EXPECT_FALSE(drv.bra_vrr(eri_term(Px, S, S, S), 'x').has_value());
}

// ket_vrr lowers center D (index 3).
TEST(T4CVrrElectronRepulsionDriverTest, KetVrrOnDShellPShell)
{
    const T4CVrrElectronRepulsionDriver drv;

    const auto rec = drv.ket_vrr(eri_term(S, S, S, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"QD", "WQ"}));
}

TEST(T4CVrrElectronRepulsionDriverTest, KetVrrOnDShellDShell)
{
    const T4CVrrElectronRepulsionDriver drv;

    const auto rec = drv.ket_vrr(eri_term(S, S, S, Dx), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 4u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"QD", "WQ", "1/nu", "rho/nu"}));
}

TEST(T4CVrrElectronRepulsionDriverTest, KetVrrRequiresDMomentum)
{
    const T4CVrrElectronRepulsionDriver drv;

    // C has momentum but D is scalar: ket_vrr lowers D.
    EXPECT_FALSE(drv.ket_vrr(eri_term(S, S, Px, S), 'x').has_value());
}

TEST(T4CVrrElectronRepulsionDriverTest, VrrRejectsNonEriTerm)
{
    const T4CVrrElectronRepulsionDriver drv;

    const auto overlap = R4CTerm(T4CIntegral(TwoCenterPairComponent({"a", "b"}, {S, Px}),
                                             TwoCenterPairComponent({"c", "d"}, {S, S}),
                                             OperatorComponent("1")));

    EXPECT_FALSE(drv.bra_vrr(overlap, 'x').has_value());
}

TEST(T4CVrrElectronRepulsionDriverTest, CreateRecursionBuildsOneExpansionPerIntegral)
{
    const T4CVrrElectronRepulsionDriver drv;

    const VT4CIntegrals vints({T4CIntegral(TwoCenterPairComponent({"a", "b"}, {S, Px}),
                                           TwoCenterPairComponent({"c", "d"}, {S, Px}),
                                           OperatorComponent("1/|r-r'|"))});

    EXPECT_EQ(drv.create_recursion(vints).expansions(), 1u);
}
