// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "t2c_trans_gen_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Recursion term for the generalized translation driver: an operator of the
// given tensor shape carrying two geometric prefixes (one per center) whose
// shapes are p0 and p1.
R2CTerm trans_term(const TensorComponent& op,
                   const TensorComponent& p0, const TensorComponent& p1)
{
    const VOperatorComponents prefixes({OperatorComponent("d/dA", p0, "bra", 0),
                                        OperatorComponent("d/dB", p1, "ket", 1)});

    return R2CTerm(T2CIntegral(OneCenterComponent("a", TensorComponent(0, 0, 0)),
                               OneCenterComponent("b", TensorComponent(0, 0, 0)),
                               OperatorComponent("G", op), 0, prefixes));
}

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);
const TensorComponent Dxx(2, 0, 0);

}  // namespace

TEST(T2CTransGenDriverTest, IsAuxilary)
{
    const T2CTransGenDriver drv;

    // Scalar operator -> auxiliary.
    EXPECT_TRUE(drv.is_auxilary(trans_term(S, S, S)));

    // Operator order 1 with prefix orders {1,0} -> NOT auxiliary (mixed g110).
    EXPECT_FALSE(drv.is_auxilary(trans_term(Px, Px, S)));

    // Operator order 1 with prefix orders {0,0} -> auxiliary.
    EXPECT_TRUE(drv.is_auxilary(trans_term(Px, S, S)));

    // Operator order 2 -> NOT auxiliary (g020).
    EXPECT_FALSE(drv.is_auxilary(trans_term(Dxx, S, S)));
}

TEST(T2CTransGenDriverTest, OperatorVrrMixedG110)
{
    const T2CTransGenDriver drv;

    // order-1 operator, prefix orders {1,0}: raises prefix 0 and prefix 1.
    const auto rec = drv.operator_vrr(trans_term(Px, Px, S));

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
}

TEST(T2CTransGenDriverTest, OperatorVrrG020)
{
    const T2CTransGenDriver drv;

    // order-2 operator, scalar prefixes: four raise-prefix combinations.
    const auto rec = drv.operator_vrr(trans_term(Dxx, S, S));

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 4u);
}

TEST(T2CTransGenDriverTest, OperatorVrrRejectsAuxilary)
{
    const T2CTransGenDriver drv;

    // Scalar operator is auxiliary -> nullopt.
    EXPECT_FALSE(drv.operator_vrr(trans_term(S, S, S)).has_value());

    // order-1 operator with {0,0} prefixes is auxiliary -> nullopt.
    EXPECT_FALSE(drv.operator_vrr(trans_term(Px, S, S)).has_value());
}
