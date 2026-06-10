// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v2i_center_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Two-center overlap integral carrying geometric-derivative prefixes. Each entry
// of `pref` is the tensorial order of the prefix operator acting on a center.
I2CIntegral geom(int a, int b, const std::vector<int>& pref)
{
    VOperators prefixes;
    for (const auto order : pref)
    {
        prefixes.push_back(Operator("d/dR", Tensor(order)));
    }

    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("1"), 0, prefixes);
}

// Plain (prefix-free) two-center overlap integral.
I2CIntegral ovl(int a, int b)
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("1"), 0, {});
}

bool contains(const SI2CIntegrals& set, const I2CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V2ICenterDriverTest, IsAuxiliary)
{
    const V2ICenterDriver drv;

    // A scalar (order-0) prefix is auxiliary at that index.
    EXPECT_TRUE(drv.is_auxiliary(geom(1, 1, {0, 1}), 0));
    EXPECT_FALSE(drv.is_auxiliary(geom(1, 1, {0, 1}), 1));

    // A rank-1 prefix is not auxiliary.
    EXPECT_TRUE(drv.is_auxiliary(geom(1, 1, {1, 0}), 1));
    EXPECT_FALSE(drv.is_auxiliary(geom(1, 1, {1, 0}), 0));
}

TEST(V2ICenterDriverTest, BraKetVrrLowersPrefixAndShiftsBra)
{
    const V2ICenterDriver drv;

    // Prefix order 1 on bra (index 0); lowering it clears all prefixes (both
    // become scalar) and yields the raised and lowered bra integrals.
    const auto tints = drv.bra_ket_vrr(geom(1, 1, {1, 0}), 0);

    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, ovl(2, 1)));
    EXPECT_TRUE(contains(tints, ovl(0, 1)));
}

TEST(V2ICenterDriverTest, BraKetVrrLowersPrefixAndShiftsKet)
{
    const V2ICenterDriver drv;

    // Prefix order 1 on ket (index 1); lowering it clears prefixes and yields
    // the raised and lowered ket integrals.
    const auto tints = drv.bra_ket_vrr(geom(1, 1, {0, 1}), 1);

    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, ovl(1, 2)));
    EXPECT_TRUE(contains(tints, ovl(1, 0)));
}

TEST(V2ICenterDriverTest, BraKetVrrLoweredBraGuard)
{
    const V2ICenterDriver drv;

    // With bra = 0, the lowering shift(-1, 0) vanishes; only the raised term
    // survives.
    const auto tints = drv.bra_ket_vrr(geom(0, 1, {1, 0}), 0);

    ASSERT_EQ(tints.size(), 1u);
    EXPECT_TRUE(contains(tints, ovl(1, 1)));
}

TEST(V2ICenterDriverTest, BraKetVrrOnAuxiliaryIsEmpty)
{
    const V2ICenterDriver drv;

    // A scalar prefix at the targeted index is auxiliary; nothing is produced.
    EXPECT_TRUE(drv.bra_ket_vrr(geom(1, 1, {0, 1}), 0).empty());
}

TEST(V2ICenterDriverTest, ApplyRecursionExpandsPrefixedIntegral)
{
    const V2ICenterDriver drv;

    // Prefix order 1 on the ket (index 1); the recursion must reduce it down to
    // prefix-free overlap base integrals.
    const auto tints = drv.apply_recursion(SI2CIntegrals({geom(1, 1, {0, 1})}));

    EXPECT_FALSE(tints.empty());
    // Reduced (prefix-free) integrals produced by lowering the ket prefix.
    EXPECT_TRUE(contains(tints, ovl(1, 2)));
    EXPECT_TRUE(contains(tints, ovl(1, 0)));
}
