// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "v2i_translation_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Two-center integral carrying a rank-`rank` operator "A" (the operator the
// translation recursion differentiates). rank 0 is auxiliary; rank >= 1 drives
// the recursion.
I2CIntegral op_int(int a, int b, int rank, const VOperators& prefixes = {})
{
    return I2CIntegral(OneCenter("a", a), OneCenter("b", b), Operator("A", Tensor(rank)), 0, prefixes);
}

bool contains(const SI2CIntegrals& set, const I2CIntegral& integral)
{
    return set.find(integral) != set.cend();
}

}  // namespace

TEST(V2ITranslationDriverTest, IsAuxiliaryFromOperatorRank)
{
    const V2ITranslationDriver drv;

    // Scalar (rank-0) operator is auxiliary; a rank-1 operator is not.
    EXPECT_TRUE(drv.is_auxiliary(op_int(0, 0, 0)));
    EXPECT_TRUE(drv.is_auxiliary(op_int(1, 0, 0)));
    EXPECT_FALSE(drv.is_auxiliary(op_int(0, 0, 1)));
    EXPECT_FALSE(drv.is_auxiliary(op_int(1, 0, 1)));
}

// Guard: an auxiliary (scalar-operator) integral produces no recursion terms.
TEST(V2ITranslationDriverTest, OperatorVrrOnAuxiliaryIsEmpty)
{
    const V2ITranslationDriver drv;

    EXPECT_TRUE(drv.operator_vrr(op_int(1, 0, 0)).empty());
}

// First-derivative branch: rank-1 operator, no prefixes, (P|S).
// rint = base with operator lowered to rank 0; then bra/ket are shifted +/-1.
// (P|S) -> {(D|S), (S|S), (P|P)}; the ket-lower term vanishes (negative ket).
// Every produced integral carries the reduced rank-0 operator.
TEST(V2ITranslationDriverTest, OperatorVrrFirstDerivativeOnPShell)
{
    const V2ITranslationDriver drv;

    const auto tints = drv.operator_vrr(op_int(1, 0, 1));

    ASSERT_EQ(tints.size(), 3u);
    EXPECT_TRUE(contains(tints, op_int(2, 0, 0)));   // bra raised:  (D|S)
    EXPECT_TRUE(contains(tints, op_int(0, 0, 0)));   // bra lowered: (S|S)
    EXPECT_TRUE(contains(tints, op_int(1, 1, 0)));   // ket raised:  (P|P)
    for (const auto& tint : tints)
    {
        EXPECT_EQ(tint.integrand().shape(), Tensor(0));
    }
}

// First-derivative branch on (S|S): only the raise terms survive (both lowers
// would go negative) -> {(P|S), (S|P)} with the reduced rank-0 operator.
TEST(V2ITranslationDriverTest, OperatorVrrFirstDerivativeOnScalarShell)
{
    const V2ITranslationDriver drv;

    const auto tints = drv.operator_vrr(op_int(0, 0, 1));

    ASSERT_EQ(tints.size(), 2u);
    EXPECT_TRUE(contains(tints, op_int(1, 0, 0)));   // bra raised: (P|S)
    EXPECT_TRUE(contains(tints, op_int(0, 1, 0)));   // ket raised: (S|P)
}

// Second-derivative branch: rank-2 operator with no prefixes. The branch relies
// on shift_prefix, which has no prefix to act on, so every term is rejected and
// the result is empty (the integral carries no prefixes to be lowered).
TEST(V2ITranslationDriverTest, OperatorVrrSecondDerivativeWithoutPrefixesIsEmpty)
{
    const V2ITranslationDriver drv;

    EXPECT_TRUE(drv.operator_vrr(op_int(0, 0, 2)).empty());
}
