// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t4c_center_driver.hpp"
#include "t4c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

// Four-center geometric-derivative ERI term. The bra pair holds centers
// (A, B), the ket pair (C, D); the four prefixes pA..pD are the per-center
// derivative shapes (prefix index 0=A, 1=B, 2=C, 3=D).
R4CTerm center_term(const TensorComponent& a, const TensorComponent& b,
                    const TensorComponent& c, const TensorComponent& d,
                    const TensorComponent& pA, const TensorComponent& pB,
                    const TensorComponent& pC, const TensorComponent& pD)
{
    const VOperatorComponents prefixes({OperatorComponent("d/dA", pA, "bra", 0),
                                        OperatorComponent("d/dB", pB, "bra", 1),
                                        OperatorComponent("d/dC", pC, "ket", 2),
                                        OperatorComponent("d/dD", pD, "ket", 3)});

    return R4CTerm(T4CIntegral(TwoCenterPairComponent({"a", "b"}, {a, b}),
                               TwoCenterPairComponent({"c", "d"}, {c, d}),
                               OperatorComponent("1/|r-r'|"), 0, prefixes));
}

T4CIntegral center_int(const TensorComponent& a, const TensorComponent& b,
                       const TensorComponent& pA)
{
    const VOperatorComponents prefixes({OperatorComponent("d/dA", pA, "bra", 0),
                                        OperatorComponent("d/dB", S, "bra", 1),
                                        OperatorComponent("d/dC", S, "ket", 2),
                                        OperatorComponent("d/dD", S, "ket", 3)});

    return T4CIntegral(TwoCenterPairComponent({"a", "b"}, {a, b}),
                       TwoCenterPairComponent({"c", "d"}, {S, S}),
                       OperatorComponent("1/|r-r'|"), 0, prefixes);
}

std::set<std::string> factor_names(const R4CDist& dist)
{
    std::set<std::string> names;
    for (const auto& factor : dist.unique_factors())
    {
        names.insert(factor.name());
    }
    return names;
}

}  // namespace

TEST(T4CCenterDriverTest, IsAuxilary)
{
    const T4CCenterDriver drv;

    // First-order derivative prefix on A -> not auxiliary at index 0.
    EXPECT_FALSE(drv.is_auxilary(center_term(S, S, S, S, Px, S, S, S), 0));
    // Scalar prefix on A -> auxiliary at index 0.
    EXPECT_TRUE(drv.is_auxilary(center_term(S, S, S, S, S, S, S, S), 0));
    // Scalar prefix on D -> auxiliary at index 3 even though A carries one.
    EXPECT_TRUE(drv.is_auxilary(center_term(S, S, S, S, Px, S, S, S), 3));
}

TEST(T4CCenterDriverTest, BraKetVrrOnScalarCenter)
{
    const T4CCenterDriver drv;

    // Prefix on A order 1, A center scalar: only the raise term survives.
    const auto rec = drv.bra_ket_vrr(center_term(S, S, S, S, Px, S, S, S), 'x', 0);

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"ba_e"}));
}

TEST(T4CCenterDriverTest, BraKetVrrOnPCenter)
{
    const T4CCenterDriver drv;

    // Prefix on A order 1, A center = P: raise term plus the (scaled) lower term.
    const auto rec = drv.bra_ket_vrr(center_term(Px, S, S, S, Px, S, S, S), 'x', 0);

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"ba_e"}));
}

TEST(T4CCenterDriverTest, BraKetVrrOnKetDCenterUsesDFactor)
{
    const T4CCenterDriver drv;

    // Prefix on D order 1, D center scalar: raise term with the d-exponent factor.
    const auto rec = drv.bra_ket_vrr(center_term(S, S, S, S, S, S, S, Px), 'x', 3);

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"kd_e"}));
}

TEST(T4CCenterDriverTest, BraKetVrrOnAuxilaryIsNullopt)
{
    const T4CCenterDriver drv;

    // Scalar prefix on A -> nothing to recurse.
    EXPECT_FALSE(drv.bra_ket_vrr(center_term(S, S, S, S, S, S, S, S), 'x', 0).has_value());
}

TEST(T4CCenterDriverTest, CreateRecursionBuildsOneExpansionPerIntegral)
{
    const T4CCenterDriver drv;

    const VT4CIntegrals vints({center_int(Px, S, Px)});

    EXPECT_EQ(drv.create_recursion(vints).expansions(), 1u);
}
