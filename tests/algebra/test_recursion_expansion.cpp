// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "recursion_expansion.hpp"
#include "integral_component.hpp"
#include "one_center_component.hpp"

namespace {

using Comp = IntegralComponent<OneCenterComponent, OneCenterComponent>;
using Term = RecursionTerm<Comp>;
using Expansion = RecursionExpansion<Comp>;

Comp make_comp(const TensorComponent& bra, int order = 0)
{
    return Comp(OneCenterComponent("a", bra), OneCenterComponent("b", TensorComponent(0, 0, 0)),
                OperatorComponent("Overlap"), order);
}

Term make_term(const TensorComponent& bra, int order = 0, const Fraction& prefactor = Fraction(1))
{
    return Term(make_comp(bra, order), {}, prefactor);
}

}  // namespace

TEST(RecursionExpansionTest, RootAndTerms)
{
    const Expansion exp(make_term(TensorComponent(2, 0, 0)),
                        {make_term(TensorComponent(1, 0, 0)), make_term(TensorComponent(0, 0, 0))});

    EXPECT_EQ(exp.root(), make_term(TensorComponent(2, 0, 0)));
    EXPECT_EQ(exp.terms(), 2u);
    EXPECT_EQ(exp[0], make_term(TensorComponent(1, 0, 0)));
}

TEST(RecursionExpansionTest, AddAppendsTerm)
{
    Expansion exp(make_term(TensorComponent(1, 0, 0)));
    EXPECT_EQ(exp.terms(), 0u);

    exp.add(make_term(TensorComponent(0, 0, 0)));
    EXPECT_EQ(exp.terms(), 1u);
}

TEST(RecursionExpansionTest, Equality)
{
    const Expansion a(make_term(TensorComponent(1, 0, 0)), {make_term(TensorComponent(0, 0, 0))});
    const Expansion b(make_term(TensorComponent(1, 0, 0)), {make_term(TensorComponent(0, 0, 0))});
    const Expansion c(make_term(TensorComponent(1, 0, 0)));

    EXPECT_EQ(a, b);
    EXPECT_NE(a, c);
}

TEST(RecursionExpansionTest, UniqueIntegralsAndNewCount)
{
    const Expansion exp(make_term(TensorComponent(2, 0, 0)),
                        {make_term(TensorComponent(1, 0, 0)),
                         make_term(TensorComponent(1, 0, 0)),   // duplicate integral
                         make_term(TensorComponent(0, 0, 0))});

    const auto integrals = exp.unique_integrals();
    EXPECT_EQ(integrals.size(), 2u);  // (1,0,0) and (0,0,0)

    // One of the two unique integrals is already known.
    const std::set<Comp> known({make_comp(TensorComponent(1, 0, 0))});
    EXPECT_EQ(exp.count_new_integrals(known), 1u);
}

TEST(RecursionExpansionTest, Split)
{
    const Expansion exp(make_term(TensorComponent(2, 0, 0)),
                        {make_term(TensorComponent(1, 0, 0)),
                         make_term(TensorComponent(1, 0, 0)),
                         make_term(TensorComponent(0, 0, 0))});

    const auto split = exp.split(make_comp(TensorComponent(1, 0, 0)));
    EXPECT_EQ(split.terms(), 2u);
    EXPECT_EQ(split.root(), exp.root());
}

TEST(RecursionExpansionTest, MinOrderSpansRootAndTerms)
{
    const Expansion exp(make_term(TensorComponent(1, 0, 0), /*order=*/2),
                        {make_term(TensorComponent(0, 0, 0), /*order=*/1),
                         make_term(TensorComponent(0, 0, 0), /*order=*/0)});

    EXPECT_EQ(exp.min_order(), 0);
}

TEST(RecursionExpansionTest, ReduceLowersOrders)
{
    Expansion exp(make_term(TensorComponent(1, 0, 0), /*order=*/2),
                  {make_term(TensorComponent(0, 0, 0), /*order=*/2)});

    exp.reduce(1);
    EXPECT_EQ(exp.root().order(), 1);
    EXPECT_EQ(exp[0].order(), 1);
}

TEST(RecursionExpansionTest, Similarity)
{
    const Expansion exp(make_term(TensorComponent(1, 0, 0)));

    EXPECT_TRUE(exp.similar(Expansion(make_term(TensorComponent(0, 1, 0)))));   // same root shell
    EXPECT_FALSE(exp.similar(Expansion(make_term(TensorComponent(2, 0, 0)))));  // higher shell
}

TEST(RecursionExpansionTest, SimplifyCombinesSameBaseTerms)
{
    Expansion exp(make_term(TensorComponent(1, 0, 0)),
                  {make_term(TensorComponent(0, 0, 0), 0, Fraction(1)),
                   make_term(TensorComponent(0, 0, 0), 0, Fraction(2))});

    exp.simplify();

    ASSERT_EQ(exp.terms(), 1u);
    EXPECT_EQ(exp[0].prefactor(), Fraction(3));
}

TEST(RecursionExpansionTest, PrefactorsAreAbsoluteAndUnique)
{
    const Expansion exp(make_term(TensorComponent(1, 0, 0)),
                        {make_term(TensorComponent(0, 0, 0), 0, Fraction(1)),
                         make_term(TensorComponent(1, 0, 0), 0, Fraction(-1))});

    const auto prefactors = exp.prefactors();
    EXPECT_EQ(prefactors.size(), 1u);  // |1| and |-1| collapse to one entry
    EXPECT_EQ(*prefactors.begin(), Fraction(1));
}
