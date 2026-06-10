// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t3c_ovl_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Three-center overlap recursion term (operator "G(r)"). The driver reuses the
// two-center recursion types, with the third (G) center carried by the operator.
R2CTerm ovl_term(const TensorComponent& bra, const TensorComponent& ket)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("G(r)")));
}

std::set<std::string> factor_names(const R2CDist& dist)
{
    std::set<std::string> names;
    for (const auto& factor : dist.unique_factors())
    {
        names.insert(factor.name());
    }
    return names;
}

}  // namespace

TEST(T3COverlapDriverTest, IsOverlap)
{
    const T3COverlapDriver drv;

    EXPECT_TRUE(drv.is_overlap(ovl_term(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0))));

    // The plain two-center overlap operator "1" is not a three-center G(r) term.
    const auto plain = R2CTerm(T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                           OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                           OperatorComponent("1")));
    EXPECT_FALSE(drv.is_overlap(plain));
}

TEST(T3COverlapDriverTest, BraVrrSingleTermOnPShell)
{
    const T3COverlapDriver drv;

    const auto rec = drv.bra_vrr(ovl_term(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"GA"}));
}

TEST(T3COverlapDriverTest, BraVrrTwoTermsOnDShell)
{
    const T3COverlapDriver drv;

    const auto rec = drv.bra_vrr(ovl_term(TensorComponent(2, 0, 0), TensorComponent(0, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"GA", "1/geta"}));
}

TEST(T3COverlapDriverTest, BraVrrWrongAxisIsNullopt)
{
    const T3COverlapDriver drv;

    EXPECT_FALSE(
        drv.bra_vrr(ovl_term(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)), 'y').has_value());
}

TEST(T3COverlapDriverTest, KetVrrSingleTermOnPShell)
{
    const T3COverlapDriver drv;

    const auto rec = drv.ket_vrr(ovl_term(TensorComponent(0, 0, 0), TensorComponent(1, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"GB"}));
}

TEST(T3COverlapDriverTest, VrrRejectsNonOverlapTerm)
{
    const T3COverlapDriver drv;

    const auto plain = R2CTerm(T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                           OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                           OperatorComponent("1")));

    EXPECT_FALSE(drv.bra_vrr(plain, 'x').has_value());
}

TEST(T3COverlapDriverTest, CreateRecursionBuildsOneExpansionPerIntegral)
{
    const T3COverlapDriver drv;

    const VT2CIntegrals vints({T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                           OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                           OperatorComponent("G(r)"))});

    EXPECT_EQ(drv.create_recursion(vints).expansions(), 1u);
}
