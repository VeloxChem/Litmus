// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t2c_proj_ecp_driver.hpp"
#include "t2c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

// Projected-ECP recursion term (operator "U_l").
R2CTerm ecp_term(const TensorComponent& bra, const TensorComponent& ket)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("U_l")));
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

TEST(T2CProjectedECPDriverTest, IsProjectedEcp)
{
    const T2CProjectedECPDriver drv;

    EXPECT_TRUE(drv.is_projected_ecp(ecp_term(Px, S)));

    const auto overlap = R2CTerm(T2CIntegral(OneCenterComponent("a", Px),
                                             OneCenterComponent("b", S),
                                             OperatorComponent("1")));
    EXPECT_FALSE(drv.is_projected_ecp(overlap));
}

TEST(T2CProjectedECPDriverTest, BraVrrCarriesRaFactor)
{
    const T2CProjectedECPDriver drv;

    const auto rec = drv.bra_vrr(ecp_term(Px, S), 'x');

    ASSERT_TRUE(rec.has_value());
    ASSERT_GT(rec->terms(), 0u);
    // The leading recursion term carries the RA distance factor.
    EXPECT_EQ(factor_names(*rec).count("RA"), 1u);
}

TEST(T2CProjectedECPDriverTest, BraVrrWrongAxisIsNullopt)
{
    const T2CProjectedECPDriver drv;

    EXPECT_FALSE(drv.bra_vrr(ecp_term(Px, S), 'y').has_value());
}

TEST(T2CProjectedECPDriverTest, BraVrrRejectsNonEcpTerm)
{
    const T2CProjectedECPDriver drv;

    const auto overlap = R2CTerm(T2CIntegral(OneCenterComponent("a", Px),
                                             OneCenterComponent("b", S),
                                             OperatorComponent("1")));

    EXPECT_FALSE(drv.bra_vrr(overlap, 'x').has_value());
}
