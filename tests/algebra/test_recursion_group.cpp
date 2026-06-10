// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "recursion_group.hpp"
#include "integral.hpp"
#include "integral_component.hpp"
#include "one_center.hpp"
#include "one_center_component.hpp"

namespace {

using Comp = IntegralComponent<OneCenterComponent, OneCenterComponent>;
using Intg = Integral<OneCenter, OneCenter>;
using Term = RecursionTerm<Comp>;
using Expansion = RecursionExpansion<Comp>;
using Group = RecursionGroup<Comp>;

Comp make_comp(const TensorComponent& bra, int order = 0)
{
    return Comp(OneCenterComponent("a", bra), OneCenterComponent("b", TensorComponent(0, 0, 0)),
                OperatorComponent("Overlap"), order);
}

Term make_term(const TensorComponent& bra, int order = 0)
{
    return Term(make_comp(bra, order));
}

Expansion make_expansion(const TensorComponent& root_bra, const VRecursionTerms<Comp>& terms = {})
{
    return Expansion(make_term(root_bra), terms);
}

}  // namespace

TEST(RecursionGroupTest, ExpansionsCountAndIndexing)
{
    const Group group({make_expansion(TensorComponent(1, 0, 0)),
                       make_expansion(TensorComponent(0, 0, 0))});

    EXPECT_EQ(group.expansions(), 2u);
    EXPECT_EQ(group[0], make_expansion(TensorComponent(1, 0, 0)));
}

TEST(RecursionGroupTest, AddAndEquality)
{
    Group group;
    EXPECT_EQ(group.expansions(), 0u);

    group.add(make_expansion(TensorComponent(1, 0, 0)));
    EXPECT_EQ(group.expansions(), 1u);

    Group other({make_expansion(TensorComponent(1, 0, 0))});
    EXPECT_EQ(group, other);
    EXPECT_NE(group, Group());
}

TEST(RecursionGroupTest, EmptyDependsOnExpansionTerms)
{
    // An expansion with no terms makes the group empty.
    EXPECT_TRUE(Group({make_expansion(TensorComponent(1, 0, 0))}).empty());

    // An expansion carrying a term does not.
    const Group nonempty(
        {make_expansion(TensorComponent(1, 0, 0), {make_term(TensorComponent(0, 0, 0))})});
    EXPECT_FALSE(nonempty.empty());
}

TEST(RecursionGroupTest, ContainsMatchesByRoot)
{
    const Group group({make_expansion(TensorComponent(1, 0, 0))});

    EXPECT_TRUE(group.contains(make_expansion(TensorComponent(1, 0, 0))));
    EXPECT_FALSE(group.contains(make_expansion(TensorComponent(2, 0, 0))));
}

TEST(RecursionGroupTest, MergeSkipsDuplicates)
{
    Group group({make_expansion(TensorComponent(1, 0, 0))});
    group.merge(Group({make_expansion(TensorComponent(1, 0, 0)),   // duplicate root
                       make_expansion(TensorComponent(0, 0, 0))}));

    EXPECT_EQ(group.expansions(), 2u);
}

TEST(RecursionGroupTest, RootsAreStrippedTerms)
{
    const Group group({make_expansion(TensorComponent(1, 0, 0)),
                       make_expansion(TensorComponent(0, 0, 0))});

    const auto roots = group.roots();
    ASSERT_EQ(roots.size(), 2u);
    EXPECT_EQ(roots[0], make_term(TensorComponent(1, 0, 0)));
}

TEST(RecursionGroupTest, MinOrderAcrossExpansions)
{
    const Group group(
        {make_expansion(TensorComponent(1, 0, 0), {make_term(TensorComponent(0, 0, 0), 2)}),
         make_expansion(TensorComponent(0, 0, 0), {make_term(TensorComponent(0, 0, 0), 1)})});

    EXPECT_EQ(group.min_order(), 0);  // roots are order 0
}

TEST(RecursionGroupTest, ComponentsCollectRootsAndTerms)
{
    const Group group({make_expansion(TensorComponent(1, 0, 0),
                                      {make_term(TensorComponent(0, 0, 0))})});

    // root (1,0,0) plus term (0,0,0).
    EXPECT_EQ(group.components().size(), 2u);
}

TEST(RecursionGroupTest, IntegralsReduceComponentsToShells)
{
    // Two roots with the same shell but different axial split collapse to one integral.
    const Group group({make_expansion(TensorComponent(1, 0, 0)),
                       make_expansion(TensorComponent(0, 1, 0))});

    const auto integrals = group.integrals<Intg>();
    EXPECT_EQ(integrals.size(), 1u);
}

TEST(RecursionGroupTest, BaseIsDefinedOnlyForUniformRoots)
{
    // Roots differ in axial split but share a shell -> single base integral.
    const Group uniform({make_expansion(TensorComponent(1, 0, 0)),
                         make_expansion(TensorComponent(0, 1, 0))});
    EXPECT_TRUE(uniform.base<Intg>().has_value());

    // Roots of different shells -> no common base.
    const Group mixed({make_expansion(TensorComponent(1, 0, 0)),
                       make_expansion(TensorComponent(2, 0, 0))});
    EXPECT_FALSE(mixed.base<Intg>().has_value());
}

TEST(RecursionGroupTest, Similarity)
{
    const Group group({make_expansion(TensorComponent(1, 0, 0))});

    EXPECT_TRUE(group.similar(Group({make_expansion(TensorComponent(0, 1, 0))})));
    EXPECT_FALSE(group.similar(Group({make_expansion(TensorComponent(1, 0, 0)),
                                      make_expansion(TensorComponent(0, 0, 0))})));  // size differs
}
