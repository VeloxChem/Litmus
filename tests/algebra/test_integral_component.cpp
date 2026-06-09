// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "integral_component.hpp"
#include "one_center_component.hpp"
#include "operator_component.hpp"
#include "tensor_component.hpp"

namespace {

// Two-center integral component, the form used by the t2c integral drivers.
using TwoCenterComp = IntegralComponent<OneCenterComponent, OneCenterComponent>;

OneCenterComponent center(const std::string& name, int ax, int ay, int az)
{
    return OneCenterComponent(name, TensorComponent(ax, ay, az));
}

OperatorComponent prefix_op(int ax, int ay, int az)
{
    return OperatorComponent("d/dA", TensorComponent(ax, ay, az), "bra", 0);
}

TwoCenterComp make_integral(const VOperatorComponents& prefixes = VOperatorComponents({}))
{
    return TwoCenterComp(center("a", 1, 0, 0),
                         center("b", 0, 0, 0),
                         OperatorComponent("Overlap", TensorComponent(0, 0, 0)),
                         0,
                         prefixes);
}

}  // namespace

TEST(IntegralComponentTest, ConstructionAndAccessors)
{
    const auto integral = make_integral();

    EXPECT_EQ(integral.order(), 0);
    EXPECT_EQ(integral.bra(), center("a", 1, 0, 0));
    EXPECT_EQ(integral.ket(), center("b", 0, 0, 0));
    EXPECT_EQ(integral.integrand().name(), "Overlap");
    EXPECT_TRUE(integral.prefixes().empty());
}

TEST(IntegralComponentTest, CenterAccess)
{
    const auto integral = make_integral();

    EXPECT_EQ(integral[0], TensorComponent(1, 0, 0));  // bra
    EXPECT_EQ(integral[1], TensorComponent(0, 0, 0));  // ket
}

// Regression test for the shift_prefix out-of-bounds fix: with no prefixes,
// index 0 is out of range and previously indexed an empty vector (UB). It must
// now return nullopt.
TEST(IntegralComponentTest, ShiftPrefixOnEmptyPrefixesReturnsNullopt)
{
    const auto integral = make_integral();

    EXPECT_FALSE(integral.shift_prefix('x', -1, 0).has_value());
}

TEST(IntegralComponentTest, ShiftPrefixOutOfRangeIndexReturnsNullopt)
{
    const auto integral = make_integral({prefix_op(1, 0, 0)});

    EXPECT_FALSE(integral.shift_prefix('x', 1, 5).has_value());
}

TEST(IntegralComponentTest, ShiftPrefixNegativeIndexReturnsNullopt)
{
    const auto integral = make_integral({prefix_op(1, 0, 0)});

    EXPECT_FALSE(integral.shift_prefix('x', 1, -1).has_value());
}

TEST(IntegralComponentTest, ShiftPrefixValidIndexRaisesOrder)
{
    const auto integral = make_integral({prefix_op(1, 0, 0)});

    const auto shifted = integral.shift_prefix('x', 1, 0);

    ASSERT_TRUE(shifted.has_value());
    ASSERT_EQ(shifted->prefixes().size(), 1u);
    EXPECT_EQ(shifted->prefixes()[0].shape().order(), 2);
}

TEST(IntegralComponentTest, ShiftPrefixToScalarErasedWhenNoscalar)
{
    const auto integral = make_integral({prefix_op(1, 0, 0)});

    const auto shifted = integral.shift_prefix('x', -1, 0, /*noscalar=*/true);

    ASSERT_TRUE(shifted.has_value());
    EXPECT_TRUE(shifted->prefixes().empty());
}

TEST(IntegralComponentTest, ShiftPrefixToScalarKeptWhenNotNoscalar)
{
    const auto integral = make_integral({prefix_op(1, 0, 0)});

    const auto shifted = integral.shift_prefix('x', -1, 0, /*noscalar=*/false);

    ASSERT_TRUE(shifted.has_value());
    ASSERT_EQ(shifted->prefixes().size(), 1u);
    EXPECT_EQ(shifted->prefixes()[0].shape().order(), 0);
}

TEST(IntegralComponentTest, BaseRemovesPrefixes)
{
    const auto integral = make_integral({prefix_op(1, 0, 0)});
    const auto base = integral.base();

    EXPECT_TRUE(base.prefixes().empty());
    EXPECT_EQ(base.bra(), integral.bra());
    EXPECT_EQ(base.ket(), integral.ket());
    EXPECT_EQ(base.order(), integral.order());
}

TEST(IntegralComponentTest, ShiftOrder)
{
    const auto integral = make_integral();

    const auto raised = integral.shift_order(2);
    ASSERT_TRUE(raised.has_value());
    EXPECT_EQ(raised->order(), 2);

    EXPECT_FALSE(integral.shift_order(-1).has_value());  // below zero
}

TEST(IntegralComponentTest, EqualityAndOrderingByOrder)
{
    const auto i0 = TwoCenterComp(center("a", 1, 0, 0), center("b", 0, 0, 0),
                                  OperatorComponent("Overlap", TensorComponent(0, 0, 0)), 0);
    const auto i1 = TwoCenterComp(center("a", 1, 0, 0), center("b", 0, 0, 0),
                                  OperatorComponent("Overlap", TensorComponent(0, 0, 0)), 1);

    EXPECT_EQ(i0, make_integral());
    EXPECT_NE(i0, i1);

    // bra/ket/integrand equal, so ordering falls through to the order field.
    EXPECT_TRUE(i0 < i1);
    EXPECT_FALSE(i1 < i0);
}

TEST(IntegralComponentTest, Similarity)
{
    const auto integral = make_integral();  // bra a(1,0,0)

    // Same name and order on bra (different axial split) -> similar.
    const auto same_shell = TwoCenterComp(center("a", 0, 1, 0), center("b", 0, 0, 0),
                                          OperatorComponent("Overlap", TensorComponent(0, 0, 0)), 0);
    EXPECT_TRUE(integral.similar(same_shell));

    // Higher bra order -> not similar.
    const auto higher = TwoCenterComp(center("a", 2, 0, 0), center("b", 0, 0, 0),
                                      OperatorComponent("Overlap", TensorComponent(0, 0, 0)), 0);
    EXPECT_FALSE(integral.similar(higher));
}
