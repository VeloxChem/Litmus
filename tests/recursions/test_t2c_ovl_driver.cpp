// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t2c_ovl_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Overlap recursion term with the given bra/ket Cartesian components.
R2CTerm ovl_term(const TensorComponent& bra, const TensorComponent& ket)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("1")));
}

// Names of the factors appearing across an expansion's terms.
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

TEST(T2COverlapDriverTest, IsOverlap)
{
    const T2COverlapDriver drv;

    EXPECT_TRUE(drv.is_overlap(ovl_term(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0))));

    // A non-overlap integrand is rejected.
    const auto kinetic = R2CTerm(T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                             OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                             OperatorComponent("T")));
    EXPECT_FALSE(drv.is_overlap(kinetic));
}

TEST(T2COverlapDriverTest, BraVrrSingleTermOnPShell)
{
    const T2COverlapDriver drv;

    const auto rec = drv.bra_vrr(ovl_term(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);                                  // only the PA term
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PA"}));
    // The root is the original (P|S) term.
    EXPECT_EQ(rec->root(), ovl_term(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)));
    EXPECT_EQ((*rec)[0].prefactor(), Fraction(1));
}

TEST(T2COverlapDriverTest, BraVrrTwoTermsOnDShell)
{
    const T2COverlapDriver drv;

    // (D|S) along x yields the PA term plus the 1/eta down-step term.
    const auto rec = drv.bra_vrr(ovl_term(TensorComponent(2, 0, 0), TensorComponent(0, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PA", "1/eta"}));
}

TEST(T2COverlapDriverTest, BraVrrWrongAxisIsNullopt)
{
    const T2COverlapDriver drv;

    // Bra has no y momentum to lower.
    EXPECT_FALSE(
        drv.bra_vrr(ovl_term(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)), 'y').has_value());
}

TEST(T2COverlapDriverTest, KetVrrSingleTermOnPShell)
{
    const T2COverlapDriver drv;

    const auto rec = drv.ket_vrr(ovl_term(TensorComponent(0, 0, 0), TensorComponent(1, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PB"}));
}

TEST(T2COverlapDriverTest, VrrRejectsNonOverlapTerm)
{
    const T2COverlapDriver drv;

    const auto kinetic = R2CTerm(T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                             OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                             OperatorComponent("T")));

    EXPECT_FALSE(drv.bra_vrr(kinetic, 'x').has_value());
}

TEST(T2COverlapDriverTest, CreateRecursionBuildsOneExpansionPerIntegral)
{
    const T2COverlapDriver drv;

    const VT2CIntegrals vints({T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                           OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                           OperatorComponent("1"))});

    const auto group = drv.create_recursion(vints);

    EXPECT_EQ(group.expansions(), 1u);
}
