// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <array>

#include "v2i_proj_ecp_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Projected-ECP integral: a Boys/projection order triple paired with a
// two-center integral over the projected-ECP operator "U_l".
M2Integral proj_ecp(int a, int b, const std::array<int, 3>& order = {0, 0, 0})
{
    return {order, I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("U_l"), 0, {})};
}

}  // namespace

TEST(V2IProjectedECPDriverTest, IsProjectedEcp)
{
    const V2IProjectedECPDriver drv;

    EXPECT_TRUE(drv.is_projected_ecp(proj_ecp(1, 0)));

    const M2Integral overlap = {{0, 0, 0},
                                I2CIntegral(OneCenter("a", 1), OneCenter("b", 0), Operator("1"), 0, {})};
    EXPECT_FALSE(drv.is_projected_ecp(overlap));
}

TEST(V2IProjectedECPDriverTest, BraVrrLowersBraAndSpansOrders)
{
    const V2IProjectedECPDriver drv;

    // (P|S) lowers the bra and emits at least the base and Boys-incremented terms.
    const auto tints = drv.bra_vrr(proj_ecp(1, 0));

    EXPECT_GE(tints.size(), 2u);
}

TEST(V2IProjectedECPDriverTest, BraVrrOnNonEcpIsEmpty)
{
    const V2IProjectedECPDriver drv;

    const M2Integral overlap = {{0, 0, 0},
                                I2CIntegral(OneCenter("a", 1), OneCenter("b", 0), Operator("1"), 0, {})};

    EXPECT_TRUE(drv.bra_vrr(overlap).empty());
}
