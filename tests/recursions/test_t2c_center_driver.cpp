// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t2c_center_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Geometric-derivative recursion term: an overlap integral carrying one prefix
// operator (a Cartesian derivative) on each center. p0/p1 are the prefix shapes
// on the bra (center 0) and ket (center 1).
R2CTerm center_term(const TensorComponent& bra, const TensorComponent& ket,
                    const TensorComponent& p0, const TensorComponent& p1)
{
    const VOperatorComponents prefixes({OperatorComponent("d/dA", p0, "bra", 0),
                                        OperatorComponent("d/dB", p1, "ket", 1)});

    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("1"), 0, prefixes));
}

T2CIntegral center_int(const TensorComponent& bra, const TensorComponent& ket,
                       const TensorComponent& p0, const TensorComponent& p1)
{
    const VOperatorComponents prefixes({OperatorComponent("d/dA", p0, "bra", 0),
                                        OperatorComponent("d/dB", p1, "ket", 1)});

    return T2CIntegral(OneCenterComponent("a", bra), OneCenterComponent("b", ket),
                       OperatorComponent("1"), 0, prefixes);
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

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

}  // namespace

TEST(T2CCenterDriverTest, IsAuxilaryFromPrefixOrder)
{
    const T2CCenterDriver drv;

    // A first-order prefix on the bra is not auxiliary.
    EXPECT_FALSE(drv.is_auxilary(center_term(S, S, Px, S), 0));
    // A scalar (already reduced) prefix is auxiliary.
    EXPECT_TRUE(drv.is_auxilary(center_term(S, S, S, S), 0));
}

TEST(T2CCenterDriverTest, BraKetVrrOnBraPrefixScalarShell)
{
    const T2CCenterDriver drv;

    // Bra prefix order 1, scalar bra: only the 2*b_e raise term survives.
    const auto rec = drv.bra_ket_vrr(center_term(S, S, Px, S), 'x', 0);

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"b_e"}));
}

TEST(T2CCenterDriverTest, BraKetVrrOnBraPrefixPShell)
{
    const T2CCenterDriver drv;

    // Bra prefix order 1 with a P bra: raise term plus the lower (scaled) term.
    const auto rec = drv.bra_ket_vrr(center_term(Px, S, Px, S), 'x', 0);

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"b_e"}));
}

TEST(T2CCenterDriverTest, BraKetVrrOnKetPrefixUsesKetFactor)
{
    const T2CCenterDriver drv;

    const auto rec = drv.bra_ket_vrr(center_term(S, S, S, Px), 'x', 1);

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"k_e"}));
}

TEST(T2CCenterDriverTest, BraKetVrrOnAuxilaryIsNullopt)
{
    const T2CCenterDriver drv;

    // Scalar prefix on the bra -> nothing to recurse.
    EXPECT_FALSE(drv.bra_ket_vrr(center_term(S, S, S, S), 'x', 0).has_value());
}

TEST(T2CCenterDriverTest, CreateRecursionBuildsOneExpansionPerIntegral)
{
    const T2CCenterDriver drv;

    const VT2CIntegrals vints({center_int(S, S, Px, Px)});

    EXPECT_EQ(drv.create_recursion(vints).expansions(), 1u);
}
