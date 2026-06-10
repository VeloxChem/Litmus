// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v2i_eri_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Two-center electron-repulsion integral with bra/ket momenta a, b and Boys order.
I2CIntegral eri(int a, int b, int order = 0)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b),
                       Operator("1/|r-r'|"), order, {});
}

bool contains(const SI2CIntegrals& set, const I2CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V2IElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const V2IElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.is_electron_repulsion(eri(1, 0)));
    // Wrong operator (overlap) is rejected.
    EXPECT_FALSE(drv.is_electron_repulsion(
        I2CIntegral(OneCenter("a", 1), OneCenter("b", 0), Operator("1"), 0, {})));
    // Prefixed integral is rejected.
    EXPECT_FALSE(drv.is_electron_repulsion(
        I2CIntegral(OneCenter("a", 1), OneCenter("b", 0), Operator("1/|r-r'|"), 0,
                    {Operator("d/dA", Tensor(1))})));
}

TEST(V2IElectronRepulsionDriverTest, BraVrrLowersBraAndRaisesOrder)
{
    const V2IElectronRepulsionDriver drv;

    // (P|S)^(0): only the order-shifted (S|S)^(1) term survives (others vanish).
    const auto tints = drv.bra_vrr(eri(1, 0));

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints, eri(0, 0, 1)));
}

TEST(V2IElectronRepulsionDriverTest, KetVrrLowersKetAndRaisesOrder)
{
    const V2IElectronRepulsionDriver drv;

    // (S|P)^(0): only the order-shifted (S|S)^(1) term survives.
    const auto tints = drv.ket_vrr(eri(0, 1));

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints, eri(0, 0, 1)));
}

TEST(V2IElectronRepulsionDriverTest, BraVrrOnWrongOperatorIsEmpty)
{
    const V2IElectronRepulsionDriver drv;

    // Overlap operator is not electron repulsion; nothing is produced.
    EXPECT_TRUE(drv.bra_vrr(
        I2CIntegral(OneCenter("a", 1), OneCenter("b", 0), Operator("1"), 0, {})).empty());
}

TEST(V2IElectronRepulsionDriverTest, CreateRecursionContainsOrderedBase)
{
    const V2IElectronRepulsionDriver drv;

    const auto tints = drv.create_recursion(SI2CIntegrals({eri(1, 0)}));

    EXPECT_FALSE(tints.empty());
    // Original integral is preserved.
    EXPECT_TRUE(contains(tints, eri(1, 0)));
    // The order-raised scalar base is produced.
    EXPECT_TRUE(contains(tints, eri(0, 0, 1)));
}

TEST(V2IElectronRepulsionDriverTest, CreateRecursionLeavesNonEriUntouched)
{
    const V2IElectronRepulsionDriver drv;

    const auto ovl = I2CIntegral(OneCenter("a", 1), OneCenter("b", 0), Operator("1"), 0, {});
    const auto tints = drv.create_recursion(SI2CIntegrals({ovl}));

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints, ovl));
}
