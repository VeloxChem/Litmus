// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t4c_hrr_eri_driver.hpp"
#include "t4c_defs.hpp"

namespace {

// Four-center electron-repulsion term: bra pair (A, B) and ket pair (C, D).
// Centers index as 0 = A, 1 = B, 2 = C, 3 = D.
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

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

}  // namespace

TEST(T4CHrrElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const T4CHrrElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.is_electron_repulsion(eri_term(Px, S, S, S)));

    const auto overlap = R4CTerm(T4CIntegral(TwoCenterPairComponent({"a", "b"}, {Px, S}),
                                             TwoCenterPairComponent({"c", "d"}, {S, S}),
                                             OperatorComponent("1")));
    EXPECT_FALSE(drv.is_electron_repulsion(overlap));
}

TEST(T4CHrrElectronRepulsionDriverTest, BraHrrMovesMomentumAToB)
{
    const T4CHrrElectronRepulsionDriver drv;

    // (A=P, B=S | S, S): HRR lowers A and raises B, with the BA distance factor.
    const auto rec = drv.bra_hrr(eri_term(Px, S, S, S), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"BA"}));
}

TEST(T4CHrrElectronRepulsionDriverTest, BraHrrWrongAxisIsNullopt)
{
    const T4CHrrElectronRepulsionDriver drv;

    EXPECT_FALSE(drv.bra_hrr(eri_term(Px, S, S, S), 'y').has_value());
}

TEST(T4CHrrElectronRepulsionDriverTest, KetHrrMovesMomentumCToD)
{
    const T4CHrrElectronRepulsionDriver drv;

    // (S, S | C=P, D=S): HRR lowers C and raises D, with the DC distance factor.
    const auto rec = drv.ket_hrr(eri_term(S, S, Px, S), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"DC"}));
}

TEST(T4CHrrElectronRepulsionDriverTest, BraHrrRejectsNonEriTerm)
{
    const T4CHrrElectronRepulsionDriver drv;

    const auto overlap = R4CTerm(T4CIntegral(TwoCenterPairComponent({"a", "b"}, {Px, S}),
                                             TwoCenterPairComponent({"c", "d"}, {S, S}),
                                             OperatorComponent("1")));

    EXPECT_FALSE(drv.bra_hrr(overlap, 'x').has_value());
}

TEST(T4CHrrElectronRepulsionDriverTest, CreateRecursionBuildsOneExpansionPerIntegral)
{
    const T4CHrrElectronRepulsionDriver drv;

    const VT4CIntegrals vints({T4CIntegral(TwoCenterPairComponent({"a", "b"}, {Px, S}),
                                           TwoCenterPairComponent({"c", "d"}, {S, S}),
                                           OperatorComponent("1/|r-r'|"))});

    EXPECT_EQ(drv.create_recursion(vints).expansions(), 1u);
}
