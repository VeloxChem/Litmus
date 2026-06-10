// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t3c_hrr_eri_driver.hpp"
#include "t3c_defs.hpp"

namespace {

// Three-center electron-repulsion term: one-center bra (A) and a two-center ket
// pair (C, D). Centers index as 0 = A, 1 = C, 2 = D.
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

TEST(T3CHrrElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const T3CHrrElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.is_electron_repulsion(
        eri_term(TensorComponent(0, 0, 0), TensorComponent(1, 0, 0), TensorComponent(0, 0, 0))));

    const auto overlap = R3CTerm(T3CIntegral(OneCenterComponent("a", TensorComponent(0, 0, 0)),
                                             TwoCenterPairComponent({"c", "d"},
                                                                    {TensorComponent(1, 0, 0),
                                                                     TensorComponent(0, 0, 0)}),
                                             OperatorComponent("1")));
    EXPECT_FALSE(drv.is_electron_repulsion(overlap));
}

TEST(T3CHrrElectronRepulsionDriverTest, KetHrrMovesMomentumCToD)
{
    const T3CHrrElectronRepulsionDriver drv;

    // (A=S | C=P, D=S): HRR lowers C and raises D, with the DC distance factor.
    const auto rec =
        drv.ket_hrr(eri_term(TensorComponent(0, 0, 0), TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"DC"}));
}

TEST(T3CHrrElectronRepulsionDriverTest, KetHrrWrongAxisIsNullopt)
{
    const T3CHrrElectronRepulsionDriver drv;

    // C has no y momentum to lower.
    EXPECT_FALSE(
        drv.ket_hrr(eri_term(TensorComponent(0, 0, 0), TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)), 'y')
            .has_value());
}

TEST(T3CHrrElectronRepulsionDriverTest, KetHrrOnScalarCIsNullopt)
{
    const T3CHrrElectronRepulsionDriver drv;

    EXPECT_FALSE(
        drv.ket_hrr(eri_term(TensorComponent(0, 0, 0), TensorComponent(0, 0, 0), TensorComponent(1, 0, 0)), 'x')
            .has_value());
}

TEST(T3CHrrElectronRepulsionDriverTest, CreateRecursionBuildsOneExpansionPerIntegral)
{
    const T3CHrrElectronRepulsionDriver drv;

    const VT3CIntegrals vints({T3CIntegral(OneCenterComponent("a", TensorComponent(0, 0, 0)),
                                           TwoCenterPairComponent({"c", "d"},
                                                                  {TensorComponent(1, 0, 0),
                                                                   TensorComponent(0, 0, 0)}),
                                           OperatorComponent("1/|r-r'|"))});

    EXPECT_EQ(drv.create_recursion(vints).expansions(), 1u);
}
