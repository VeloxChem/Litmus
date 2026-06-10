// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t2c_dip_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Multipole recursion term (operator "r"); the operator shape carries the
// multipole component (e.g. (1,0,0) is the x dipole).
R2CTerm mpole_term(const TensorComponent& bra, const TensorComponent& op, const TensorComponent& ket)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("r", op)));
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

TEST(T2CMultipoleDriverTest, IsMultipole)
{
    const T2CMultipoleDriver drv;

    // Dipole: operator "r" with a non-scalar shape.
    EXPECT_TRUE(drv.is_multipole(mpole_term(Px, Px, S)));

    // Scalar "r" shape is not a multipole.
    EXPECT_FALSE(drv.is_multipole(mpole_term(Px, S, S)));

    // Wrong operator name.
    const auto overlap = R2CTerm(T2CIntegral(OneCenterComponent("a", Px),
                                             OneCenterComponent("b", S),
                                             OperatorComponent("1")));
    EXPECT_FALSE(drv.is_multipole(overlap));
}

TEST(T2CMultipoleDriverTest, BraVrrStepsBraAndOperator)
{
    const T2CMultipoleDriver drv;

    // (Px | r_x | S): the PA bra step and the operator down-step (r -> 1).
    const auto rec = drv.bra_vrr(mpole_term(Px, Px, S), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PA", "1/eta"}));
    // The operator down-step collapses the dipole "r" to the overlap "1".
    EXPECT_EQ(term_operators(*rec), std::set<std::string>({"r", "1"}));
}

TEST(T2CMultipoleDriverTest, BraVrrRequiresBraMomentum)
{
    const T2CMultipoleDriver drv;

    // No bra momentum to lower (the operator alone cannot drive the bra VRR).
    EXPECT_FALSE(drv.bra_vrr(mpole_term(S, Px, S), 'x').has_value());
}

TEST(T2CMultipoleDriverTest, VrrRejectsNonMultipoleTerm)
{
    const T2CMultipoleDriver drv;

    const auto overlap = R2CTerm(T2CIntegral(OneCenterComponent("a", Px),
                                             OneCenterComponent("b", S),
                                             OperatorComponent("1")));

    EXPECT_FALSE(drv.bra_vrr(overlap, 'x').has_value());
}

TEST(T2CMultipoleDriverTest, CreateRecursionBuildsOneExpansionPerIntegral)
{
    const T2CMultipoleDriver drv;

    const VT2CIntegrals vints({T2CIntegral(OneCenterComponent("a", Px),
                                           OneCenterComponent("b", S),
                                           OperatorComponent("r", Px))});

    EXPECT_EQ(drv.create_recursion(vints).expansions(), 1u);
}
