// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <algorithm>
#include <set>
#include <string>

#include "t2c_npot_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Nuclear-potential recursion term (operator "A", scalar shape).
R2CTerm npot_term(const TensorComponent& bra, const TensorComponent& ket)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("A")));
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

// Highest Boys order across an expansion's terms (the nuclear recursion bumps
// the order on the PC / penetration terms).
int max_term_order(const R2CDist& dist)
{
    int order = dist.root().order();
    for (size_t i = 0; i < dist.terms(); i++)
    {
        order = std::max(order, dist[i].order());
    }
    return order;
}

}  // namespace

TEST(T2CNuclearPotentialDriverTest, IsNuclearPotential)
{
    const T2CNuclearPotentialDriver drv;

    EXPECT_TRUE(
        drv.is_nuclear_potential(npot_term(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0))));

    // Wrong operator name.
    const auto overlap = R2CTerm(T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                             OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                             OperatorComponent("1")));
    EXPECT_FALSE(drv.is_nuclear_potential(overlap));

    // Right name but non-scalar operator shape.
    const auto tensorA = R2CTerm(T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                             OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                             OperatorComponent("A", TensorComponent(1, 0, 0))));
    EXPECT_FALSE(drv.is_nuclear_potential(tensorA));
}

TEST(T2CNuclearPotentialDriverTest, BraVrrOnPShell)
{
    const T2CNuclearPotentialDriver drv;

    // (P|S): PA term at Boys order 0 and PC term at Boys order 1.
    const auto rec = drv.bra_vrr(npot_term(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PA", "PC"}));
    EXPECT_EQ(max_term_order(*rec), 1);  // the PC term lives one Boys order up
}

TEST(T2CNuclearPotentialDriverTest, BraVrrOnDShell)
{
    const T2CNuclearPotentialDriver drv;

    const auto rec = drv.bra_vrr(npot_term(TensorComponent(2, 0, 0), TensorComponent(0, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 4u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PA", "PC", "1/eta"}));
}

TEST(T2CNuclearPotentialDriverTest, KetVrrOnPShell)
{
    const T2CNuclearPotentialDriver drv;

    const auto rec = drv.ket_vrr(npot_term(TensorComponent(0, 0, 0), TensorComponent(1, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PB", "PC"}));
}

TEST(T2CNuclearPotentialDriverTest, VrrRejectsNonNuclearTerm)
{
    const T2CNuclearPotentialDriver drv;

    const auto overlap = R2CTerm(T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                             OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                             OperatorComponent("1")));

    EXPECT_FALSE(drv.bra_vrr(overlap, 'x').has_value());
}

TEST(T2CNuclearPotentialDriverTest, CreateRecursionBuildsOneExpansionPerIntegral)
{
    const T2CNuclearPotentialDriver drv;

    const VT2CIntegrals vints({T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                           OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                           OperatorComponent("A"))});

    EXPECT_EQ(drv.create_recursion(vints).expansions(), 1u);
}
