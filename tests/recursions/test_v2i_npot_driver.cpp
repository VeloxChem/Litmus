// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v2i_npot_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Two-center nuclear-potential integral (operator "A") at the given Boys order.
I2CIntegral npot(int a, int b, int order = 0)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("A"), order, {});
}

bool contains(const SI2CIntegrals& set, const I2CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V2INuclearPotentialDriverTest, IsNuclearPotential)
{
    const V2INuclearPotentialDriver drv;

    EXPECT_TRUE(drv.is_nuclear_potential(npot(1, 0)));
    EXPECT_FALSE(drv.is_nuclear_potential(
        I2CIntegral(OneCenter("a", 1), OneCenter("b", 0), Operator("1"), 0, {})));
}

TEST(V2INuclearPotentialDriverTest, BraVrrSpawnsBoysOrders)
{
    const V2INuclearPotentialDriver drv;

    // (P|S) lowers to (S|S) at Boys orders 0 and 1.
    const auto tints = drv.bra_vrr(npot(1, 0));

    EXPECT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, npot(0, 0, 0)));
    EXPECT_TRUE(contains(tints, npot(0, 0, 1)));
}

TEST(V2INuclearPotentialDriverTest, BraVrrOnScalarIsEmpty)
{
    const V2INuclearPotentialDriver drv;

    EXPECT_TRUE(drv.bra_vrr(npot(0, 0)).empty());
}

TEST(V2INuclearPotentialDriverTest, CreateRecursionContainsBase)
{
    const V2INuclearPotentialDriver drv;

    const auto tints = drv.create_recursion(SI2CIntegrals({npot(1, 0)}));

    ASSERT_FALSE(tints.empty());
    // The original (P|S) and the fully reduced (S|S) Boys-order-0 base appear.
    // (The recursion also emits a plain-overlap "1" auxiliary at the bottom.)
    EXPECT_TRUE(contains(tints, npot(1, 0)));
    EXPECT_TRUE(contains(tints, npot(0, 0, 0)));
}
