// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t2c_loc_ecp_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Local-ECP recursion term (operator "U_L").
R2CTerm ecp_term(const TensorComponent& bra, const TensorComponent& ket)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("U_L")));
}

T2CIntegral ecp_int(const TensorComponent& bra, const TensorComponent& ket)
{
    return T2CIntegral(OneCenterComponent("a", bra),
                       OneCenterComponent("b", ket),
                       OperatorComponent("U_L"));
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

TEST(T2CLocalECPDriverTest, IsLocalEcp)
{
    const T2CLocalECPDriver drv;

    EXPECT_TRUE(drv.is_local_ecp(ecp_term(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0))));

    const auto overlap = R2CTerm(T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                             OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                             OperatorComponent("1")));
    EXPECT_FALSE(drv.is_local_ecp(overlap));
}

TEST(T2CLocalECPDriverTest, FullBraVrrOnPShell)
{
    const T2CLocalECPDriver drv;

    const auto rec = drv.full_bra_vrr(ecp_term(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"RA"}));
}

TEST(T2CLocalECPDriverTest, FullBraVrrOnDShell)
{
    const T2CLocalECPDriver drv;

    const auto rec = drv.full_bra_vrr(ecp_term(TensorComponent(2, 0, 0), TensorComponent(0, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"RA", "1/xi"}));
}

TEST(T2CLocalECPDriverTest, CreateFullRecursionBuildsOneExpansionPerIntegral)
{
    const T2CLocalECPDriver drv;

    const VT2CIntegrals vints({ecp_int(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0))});

    EXPECT_EQ(drv.create_full_recursion(vints).expansions(), 1u);
}

// Regression for the apply_full_recursion(R2Group&) fix: the group overload must
// route through the *full* per-expansion recursion (it previously called the
// plain apply_recursion, dropping a term). Its integral set must therefore match
// create_full_recursion, which uses the correct full path.
TEST(T2CLocalECPDriverTest, ApplyFullRecursionGroupMatchesCreateFull)
{
    const T2CLocalECPDriver drv;

    const VT2CIntegrals vints({ecp_int(TensorComponent(1, 0, 0), TensorComponent(1, 0, 0))});

    R2Group group;
    group.add(R2CDist(R2CTerm(vints[0])));
    drv.apply_full_recursion(group);

    const auto reference = drv.create_full_recursion(vints);

    EXPECT_FALSE(group.components().empty());
    EXPECT_EQ(group.components(), reference.components());
}
