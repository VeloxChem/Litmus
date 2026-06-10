// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t2c_trans_one_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Translation recursion term: the operator shape carries the rank (the driver
// acts on rank-1 integrands and reduces to the "R" operator).
R2CTerm trans_term(const TensorComponent& bra, const TensorComponent& ket, const TensorComponent& op)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("R", op)));
}

T2CIntegral trans_int(const TensorComponent& bra, const TensorComponent& ket, const TensorComponent& op)
{
    return T2CIntegral(OneCenterComponent("a", bra), OneCenterComponent("b", ket),
                       OperatorComponent("R", op));
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

std::set<std::string> term_operators(const R2CDist& dist)
{
    std::set<std::string> names;
    for (size_t i = 0; i < dist.terms(); i++)
    {
        names.insert(dist[i].integral().integrand().name());
    }
    return names;
}

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

}  // namespace

TEST(T2CTransOneDriverTest, IsAuxilaryFromOperatorRank)
{
    const T2CTransOneDriver drv;

    // A rank-1 operand drives the recursion; a scalar one is auxiliary.
    EXPECT_FALSE(drv.is_auxilary(trans_term(S, S, Px)));
    EXPECT_TRUE(drv.is_auxilary(trans_term(S, S, S)));
}

TEST(T2CTransOneDriverTest, BraKetVrrOnScalarShells)
{
    const T2CTransOneDriver drv;

    // (S|S): the 2*b_e bra-raise and 2*k_e ket-raise terms, both over operator R.
    const auto rec = drv.bra_ket_vrr(trans_term(S, S, Px));

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"b_e", "k_e"}));
    EXPECT_EQ(term_operators(*rec), std::set<std::string>({"R"}));
}

TEST(T2CTransOneDriverTest, BraKetVrrOnPShells)
{
    const T2CTransOneDriver drv;

    // (P|P): adds the two lower (scaled) terms on top of the raise terms.
    const auto rec = drv.bra_ket_vrr(trans_term(Px, Px, Px));

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 4u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"b_e", "k_e"}));
}

TEST(T2CTransOneDriverTest, BraKetVrrOnAuxilaryIsNullopt)
{
    const T2CTransOneDriver drv;

    EXPECT_FALSE(drv.bra_ket_vrr(trans_term(S, S, S)).has_value());
}

// Regression for the create_recursion fix: the optional from bra_ket_vrr is now
// guarded, so an auxiliary (scalar-operand) integral is skipped rather than
// dereferenced. Previously this dereferenced an empty optional.
TEST(T2CTransOneDriverTest, CreateRecursionSkipsAuxilaryIntegrals)
{
    const T2CTransOneDriver drv;

    const VT2CIntegrals vints({trans_int(S, S, Px),   // rank-1: produces an expansion
                               trans_int(S, S, S)});  // scalar: auxiliary, skipped

    const auto group = drv.create_recursion(vints);

    EXPECT_EQ(group.expansions(), 1u);
}
