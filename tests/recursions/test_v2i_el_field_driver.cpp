// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>

#include "v2i_el_field_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Two-center electric-field integral (operator "AG" of the given tensor rank).
I2CIntegral efield(int a, int b, int rank, int order = 0)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("AG", Tensor(rank)), order, {});
}

std::set<int> orders_of(const SI2CIntegrals& tints)
{
    std::set<int> orders;
    for (const auto& tint : tints)
    {
        orders.insert(tint.order());
    }
    return orders;
}

}  // namespace

TEST(V2IElectricFieldDriverTest, IsElectricField)
{
    const V2IElectricFieldDriver drv;

    EXPECT_TRUE(drv.is_electric_field(efield(0, 0, 1)));
    EXPECT_FALSE(drv.is_electric_field(
        I2CIntegral(OneCenter("a", 0), OneCenter("b", 0), Operator("1"), 0, {})));
}

TEST(V2IElectricFieldDriverTest, AuxVrrRankOne)
{
    const V2IElectricFieldDriver drv;

    const auto tints = drv.aux_vrr(efield(0, 0, 1));

    EXPECT_EQ(orders_of(tints), std::set<int>({1}));
    for (const auto& tint : tints)
    {
        EXPECT_EQ(tint.integrand(), Operator("A"));  // collapsed to scalar A
    }
}

TEST(V2IElectricFieldDriverTest, AuxVrrRankTwo)
{
    const V2IElectricFieldDriver drv;

    EXPECT_EQ(orders_of(drv.aux_vrr(efield(0, 0, 2))), std::set<int>({1, 2}));
}

// Regression: rank >= 3 previously fell through the TODO branch and returned an
// empty set; aux_vrr must now emit the full order range 1..N.
TEST(V2IElectricFieldDriverTest, AuxVrrRankThreeIsNotEmpty)
{
    const V2IElectricFieldDriver drv;

    const auto tints = drv.aux_vrr(efield(0, 0, 3));

    EXPECT_EQ(tints.size(), 3u);
    EXPECT_EQ(orders_of(tints), std::set<int>({1, 2, 3}));
    for (const auto& tint : tints)
    {
        EXPECT_EQ(tint.integrand(), Operator("A"));
    }
}

TEST(V2IElectricFieldDriverTest, AuxVrrRequiresScalarBraKet)
{
    const V2IElectricFieldDriver drv;

    // Auxiliaries are only generated from (S|S); a non-scalar bra yields nothing.
    EXPECT_TRUE(drv.aux_vrr(efield(1, 0, 2)).empty());
}

TEST(V2IElectricFieldDriverTest, AuxVrrRespectsStartingOrder)
{
    const V2IElectricFieldDriver drv;

    // Starting from order 2, a rank-2 operator yields orders 3 and 4.
    EXPECT_EQ(orders_of(drv.aux_vrr(efield(0, 0, 2, /*order=*/2))), std::set<int>({3, 4}));
}
