// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t2c_linmom_driver.hpp"
#include "t2c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

// Linear-momentum recursion term (operator "p"); the operator shape carries the
// Cartesian component of the momentum operator.
R2CTerm linmom_term(const TensorComponent& bra, const TensorComponent& ket, const TensorComponent& op)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("p", op)));
}

std::set<std::string> term_operators(const R2CDist& dist)
{
    std::set<std::string> names;
    for (size_t i = 0; i < dist.terms(); i++)
    {
        names.insert(dist[i].integral().integrand().name());
    }
    return names;
}

}  // namespace

TEST(T2CLinearMomentumDriverTest, IsLinearMomentum)
{
    const T2CLinearMomentumDriver drv;

    EXPECT_TRUE(drv.is_linear_momentum(linmom_term(S, S, Px)));
    // Scalar operator shape is not a (rank-1) momentum operator.
    EXPECT_FALSE(drv.is_linear_momentum(linmom_term(S, S, S)));
}

TEST(T2CLinearMomentumDriverTest, OpVrrReducesToOverlap)
{
    const T2CLinearMomentumDriver drv;

    // (S | p_x | S): only the ket-raise term survives, over the overlap operator.
    const auto rec = drv.op_vrr(linmom_term(S, S, Px));

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
    EXPECT_EQ(term_operators(*rec), std::set<std::string>({"1"}));
}

TEST(T2CLinearMomentumDriverTest, OpVrrOnKetPShellGivesTwoTerms)
{
    const T2CLinearMomentumDriver drv;

    // (S | p_x | P): both the ket-raise and ket-lower terms, over overlap.
    const auto rec = drv.op_vrr(linmom_term(S, Px, Px));

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(term_operators(*rec), std::set<std::string>({"1"}));
}

TEST(T2CLinearMomentumDriverTest, OpVrrRejectsNonMomentumTerm)
{
    const T2CLinearMomentumDriver drv;

    EXPECT_FALSE(drv.op_vrr(linmom_term(S, S, S)).has_value());  // scalar operator
}
