// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "recursion_term.hpp"
#include "integral_component.hpp"
#include "one_center_component.hpp"

namespace {

using Comp = IntegralComponent<OneCenterComponent, OneCenterComponent>;
using Term = RecursionTerm<Comp>;

OneCenterComponent oc(const std::string& name, int ax, int ay, int az)
{
    return OneCenterComponent(name, TensorComponent(ax, ay, az));
}

Comp make_comp(const TensorComponent& bra, const TensorComponent& ket, int order = 0)
{
    return Comp(OneCenterComponent("a", bra), OneCenterComponent("b", ket),
                OperatorComponent("Overlap"), order);
}

Factor make_factor()
{
    return Factor("X", "X", TensorComponent(1, 0, 0));
}

}  // namespace

TEST(RecursionTermTest, DefaultPrefactorIsUnity)
{
    const Term term(make_comp(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)));

    EXPECT_EQ(term.prefactor(), Fraction(1));
    EXPECT_EQ(term.order(), 0);
    EXPECT_TRUE(term.factors().empty());
}

TEST(RecursionTermTest, Accessors)
{
    const Term term(make_comp(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)), {}, Fraction(1, 2));

    EXPECT_EQ(term.prefactor(), Fraction(1, 2));
    EXPECT_EQ(term[0], TensorComponent(1, 0, 0));            // bra shape
    EXPECT_EQ(term[1], TensorComponent(0, 0, 0));            // ket shape
    EXPECT_EQ(term.bra<OneCenterComponent>(), oc("a", 1, 0, 0));
    EXPECT_EQ(term.ket<OneCenterComponent>(), oc("b", 0, 0, 0));
    EXPECT_EQ(term.integrand().name(), "Overlap");
    EXPECT_EQ(term.label(), "x_0");  // component label: bra axial "x", ket scalar "0"
}

TEST(RecursionTermTest, EqualityAndOrdering)
{
    const Term a(make_comp(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)));
    const Term b(make_comp(TensorComponent(2, 0, 0), TensorComponent(0, 0, 0)));

    EXPECT_EQ(a, Term(make_comp(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0))));
    EXPECT_NE(a, b);
    EXPECT_TRUE(a < b);
}

TEST(RecursionTermTest, SimilarityAndSameBase)
{
    const Term a(make_comp(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)));

    // Same shell, different axial split -> similar but not equal.
    const Term split(make_comp(TensorComponent(0, 1, 0), TensorComponent(0, 0, 0)));
    EXPECT_TRUE(a.similar(split));
    EXPECT_FALSE(a.same_base(split));

    // same_base ignores prefactor.
    const Term scaled(make_comp(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)), {}, Fraction(3));
    EXPECT_TRUE(a.same_base(scaled));
    EXPECT_NE(a, scaled);
}

TEST(RecursionTermTest, AddFactorTracksOrderAndScalesPrefactor)
{
    Term term(make_comp(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)));
    const auto factor = make_factor();

    term.add(factor, Fraction(2));
    EXPECT_EQ(term.factor_order(factor), 1);
    EXPECT_EQ(term.prefactor(), Fraction(2));

    term.add(factor, Fraction(3));
    EXPECT_EQ(term.factor_order(factor), 2);
    EXPECT_EQ(term.prefactor(), Fraction(6));

    EXPECT_EQ(term.factors().size(), 1u);
}

TEST(RecursionTermTest, ScaleAndAddFraction)
{
    Term term(make_comp(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)));

    term.scale(Fraction(1, 2));
    EXPECT_EQ(term.prefactor(), Fraction(1, 2));

    term.add(Fraction(4));  // fraction-only overload only touches the prefactor
    EXPECT_EQ(term.prefactor(), Fraction(2));
}

TEST(RecursionTermTest, Remove)
{
    const Factor keep("Y", "Y");
    Term term(make_comp(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)));
    term.add(make_factor());
    term.add(keep);

    const auto pruned = term.remove("X");
    EXPECT_EQ(pruned.factor_order(make_factor()), 0);
    EXPECT_EQ(pruned.factor_order(keep), 1);
}

TEST(RecursionTermTest, Replace)
{
    const Term term(make_comp(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)));

    const auto replaced = term.replace(OperatorComponent("Kinetic"));
    EXPECT_EQ(replaced.integrand().name(), "Kinetic");
}

TEST(RecursionTermTest, ShiftRaisesCenter)
{
    const Term term(make_comp(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)));

    const auto raised = term.shift('x', 1, 0);
    ASSERT_TRUE(raised.has_value());
    EXPECT_EQ((*raised)[0], TensorComponent(2, 0, 0));

    EXPECT_FALSE(term.shift('x', -2, 0).has_value());  // below zero
}

TEST(RecursionTermTest, ShiftOrder)
{
    const Term term(make_comp(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0), 1));

    const auto lowered = term.shift_order(-1);
    ASSERT_TRUE(lowered.has_value());
    EXPECT_EQ(lowered->order(), 0);

    EXPECT_FALSE(lowered->shift_order(-1).has_value());
}

TEST(RecursionTermTest, AuxilaryAtScalarCenter)
{
    const Term term(make_comp(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)));

    EXPECT_FALSE(term.auxilary(0));  // bra is P
    EXPECT_TRUE(term.auxilary(1));   // ket is S
}
