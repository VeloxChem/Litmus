// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v2i_kin_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Two-center kinetic-energy integral (operator "T") with momenta a and b.
I2CIntegral kin(int a, int b, int order = 0)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("T"), order, {});
}

// Two-center overlap integral (operator "1"), the base of the kinetic recursion.
I2CIntegral ovl(int a, int b, int order = 0)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("1"), order, {});
}

bool contains(const SI2CIntegrals& set, const I2CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V2IKineticEnergyDriverTest, IsKineticEnergy)
{
    const V2IKineticEnergyDriver drv;

    EXPECT_TRUE(drv.is_kinetic_energy(kin(1, 0)));
    // Wrong operator (overlap) is rejected.
    EXPECT_FALSE(drv.is_kinetic_energy(ovl(1, 0)));
    EXPECT_FALSE(drv.is_kinetic_energy(
        I2CIntegral(OneCenter("a", 1), OneCenter("b", 0), Operator("p", Tensor(1)), 0, {})));
}

// bra_vrr on (P|S): tval = (S|S,T). Of the five recursion terms only two survive
// (the shift(-1,*) lower terms would go negative): the (S|S) kinetic term and the
// (P|S) overlap replacement (integral.replace("1")).
TEST(V2IKineticEnergyDriverTest, BraVrrOnPShellReducesToOverlapBase)
{
    const V2IKineticEnergyDriver drv;

    const auto tints = drv.bra_vrr(kin(1, 0));

    ASSERT_FALSE(tints.empty());
    EXPECT_TRUE(contains(tints, kin(0, 0)));   // (S|S) kinetic, tval
    EXPECT_TRUE(contains(tints, ovl(1, 0)));   // (P|S) overlap replacement
    EXPECT_EQ(tints.size(), 2u);
}

// ket_vrr on (S|P): tval = (S|S,T). Only (S|S) kinetic and the (S|P) overlap
// replacement survive; the ket-lower terms vanish (negative ket).
TEST(V2IKineticEnergyDriverTest, KetVrrOnPShellReducesToOverlapBase)
{
    const V2IKineticEnergyDriver drv;

    const auto tints = drv.ket_vrr(kin(0, 1));

    ASSERT_FALSE(tints.empty());
    EXPECT_TRUE(contains(tints, kin(0, 0)));   // (S|S) kinetic, tval
    EXPECT_TRUE(contains(tints, ovl(0, 1)));   // (S|P) overlap replacement
    EXPECT_EQ(tints.size(), 2u);
}

// Guard: nothing to lower on a scalar (S|S) shell.
TEST(V2IKineticEnergyDriverTest, BraVrrOnScalarIsEmpty)
{
    const V2IKineticEnergyDriver drv;

    EXPECT_TRUE(drv.bra_vrr(kin(0, 0)).empty());
}

TEST(V2IKineticEnergyDriverTest, KetVrrOnScalarIsEmpty)
{
    const V2IKineticEnergyDriver drv;

    EXPECT_TRUE(drv.ket_vrr(kin(0, 0)).empty());
}

// Guard: vrr on a non-kinetic operator yields nothing.
TEST(V2IKineticEnergyDriverTest, VrrOnNonKineticIsEmpty)
{
    const V2IKineticEnergyDriver drv;

    EXPECT_TRUE(drv.bra_vrr(ovl(1, 0)).empty());
    EXPECT_TRUE(drv.ket_vrr(ovl(0, 1)).empty());
}

// create_recursion smoke test: expansion is non-empty, contains the original
// kinetic integral, and reduces down to the overlap "1" base.
TEST(V2IKineticEnergyDriverTest, CreateRecursionReducesToOverlapBase)
{
    const V2IKineticEnergyDriver drv;

    const auto tints = drv.create_recursion(SI2CIntegrals({kin(1, 0)}));

    ASSERT_FALSE(tints.empty());
    EXPECT_TRUE(contains(tints, kin(1, 0)));   // original retained
    EXPECT_TRUE(contains(tints, ovl(0, 0)));   // reduced overlap base
}

// A non-kinetic integral passes through create_recursion unchanged.
TEST(V2IKineticEnergyDriverTest, CreateRecursionPassesThroughNonKinetic)
{
    const V2IKineticEnergyDriver drv;

    const auto tints = drv.create_recursion(SI2CIntegrals({ovl(1, 0)}));

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints, ovl(1, 0)));
}
