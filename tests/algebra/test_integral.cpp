// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "integral.hpp"
#include "one_center.hpp"

namespace {

// Two-center integral over one-center expansions, the form used by t2c drivers.
using TwoCenter = Integral<OneCenter, OneCenter>;

Operator prefix_op(const int order)
{
    return Operator("d/dA", Tensor(order), "bra", 0);
}

TwoCenter make_integral(const int bramom = 1,
                        const int ketmom = 0,
                        const VOperators& prefixes = VOperators({}))
{
    return TwoCenter(OneCenter("a", bramom), OneCenter("b", ketmom),
                     Operator("Overlap", Tensor(0)), 0, prefixes);
}

}  // namespace

TEST(IntegralTest, ConstructionAndAccessors)
{
    const auto integral = make_integral();

    EXPECT_EQ(integral.order(), 0);
    EXPECT_EQ(integral.integrand().name(), "Overlap");
    EXPECT_TRUE(integral.is_simple());
    EXPECT_TRUE(integral.prefixes().empty());
    EXPECT_EQ(integral[0], 1);  // bra angular momentum
    EXPECT_EQ(integral[1], 0);  // ket angular momentum
}

TEST(IntegralTest, Equality)
{
    EXPECT_EQ(make_integral(1, 0), make_integral(1, 0));
    EXPECT_NE(make_integral(1, 0), make_integral(2, 0));

    auto raised = make_integral();
    raised.set_order(1);
    EXPECT_NE(raised, make_integral());
}

TEST(IntegralTest, OrderingFallsThroughToOrder)
{
    auto i0 = make_integral();
    auto i1 = make_integral();
    i1.set_order(1);

    EXPECT_TRUE(i0 < i1);
    EXPECT_FALSE(i1 < i0);
}

TEST(IntegralTest, Label)
{
    EXPECT_EQ(make_integral(1, 0).label(), "PS");
    EXPECT_EQ(make_integral(2, 1).label(), "DP");
    EXPECT_EQ(make_integral(1, 0).label(/*use_order=*/true), "PS_0");
}

TEST(IntegralTest, ShiftRaisesBraAndKet)
{
    const auto integral = make_integral(1, 0);

    const auto bra = integral.shift(1, 0);
    ASSERT_TRUE(bra.has_value());
    EXPECT_EQ((*bra)[0], 2);

    const auto ket = integral.shift(1, 1);
    ASSERT_TRUE(ket.has_value());
    EXPECT_EQ((*ket)[1], 1);

    EXPECT_FALSE(integral.shift(-1, 1).has_value());  // ket already S
}

TEST(IntegralTest, ShiftOrder)
{
    const auto integral = make_integral();

    const auto raised = integral.shift_order(2);
    ASSERT_TRUE(raised.has_value());
    EXPECT_EQ(raised->order(), 2);

    EXPECT_FALSE(integral.shift_order(-1).has_value());
}

TEST(IntegralTest, ShiftOperatorRaisesIntegrandShell)
{
    const auto integral = make_integral();

    const auto shifted = integral.shift_operator(1);
    ASSERT_TRUE(shifted.has_value());
    EXPECT_EQ(shifted->integrand().shape().order(), 1);
}

TEST(IntegralTest, PrefixesAndReduction)
{
    const auto integral = make_integral(1, 0, {prefix_op(1)});

    EXPECT_FALSE(integral.is_simple());
    EXPECT_EQ(integral.prefixes_order(), std::vector<int>({1}));
    EXPECT_EQ(integral.prefixes_sum_order(), 1);
    EXPECT_EQ(integral.prefix_label(), "g{1}");

    // base() strips prefixes.
    EXPECT_TRUE(integral.base().is_simple());

    // A scalar prefix is cleared by reduce_prefixes; a tensor prefix is kept.
    auto scalar = make_integral(1, 0, {prefix_op(0)});
    scalar.reduce_prefixes();
    EXPECT_TRUE(scalar.is_simple());

    auto tensor = make_integral(1, 0, {prefix_op(1)});
    tensor.reduce_prefixes();
    EXPECT_FALSE(tensor.is_simple());
}

TEST(IntegralTest, ShiftPrefix)
{
    const auto integral = make_integral(1, 0, {prefix_op(1)});

    const auto raised = integral.shift_prefix(1, 0, false);
    ASSERT_TRUE(raised.has_value());
    EXPECT_EQ(raised->prefixes_order(), std::vector<int>({2}));

    // Lowering to scalar with noscalar erases the prefix entirely.
    const auto erased = integral.shift_prefix(-1, 0, /*noscalar=*/true);
    ASSERT_TRUE(erased.has_value());
    EXPECT_TRUE(erased->is_simple());

    // Out-of-range index is rejected.
    EXPECT_FALSE(integral.shift_prefix(1, 5, false).has_value());
}

TEST(IntegralTest, ComponentsExpandOverShells)
{
    // bra P (3) x ket S (1) x scalar integrand (1) = 3 components.
    const auto integral = make_integral(1, 0);
    const auto comps = integral.components<OneCenterComponent, OneCenterComponent>();

    EXPECT_EQ(comps.size(), 3u);
}

TEST(IntegralTest, DiagComponentsPairEqualShells)
{
    // Diagonal pairing requires identical bra/ket component counts; P x P -> 3.
    const auto integral = make_integral(1, 1);
    const auto comps = integral.diag_components<OneCenterComponent, OneCenterComponent>();

    EXPECT_EQ(comps.size(), 3u);
}
