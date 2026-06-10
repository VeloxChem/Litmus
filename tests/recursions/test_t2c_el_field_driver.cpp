// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t2c_el_field_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Electric-field recursion term: operator "AG" with the given (non-scalar)
// tensor shape, carrying no geometric prefixes.
R2CTerm field_term(const TensorComponent& bra, const TensorComponent& ket,
                   const TensorComponent& op)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("AG", op)));
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

TEST(T2CElectricFieldDriverTest, IsElectricField)
{
    const T2CElectricFieldDriver drv;

    // Operator "AG" with a non-scalar shape -> electric field term.
    EXPECT_TRUE(drv.is_electric_field(field_term(Px, S, Px)));

    // Wrong operator name.
    EXPECT_FALSE(drv.is_electric_field(
        R2CTerm(T2CIntegral(OneCenterComponent("a", Px), OneCenterComponent("b", S),
                            OperatorComponent("A", Px)))));

    // Right name but scalar operator shape -> not an electric field term.
    EXPECT_FALSE(drv.is_electric_field(field_term(Px, S, S)));

    // Right name and shape but carries a geometric prefix -> rejected.
    const VOperatorComponents prefixes({OperatorComponent("d/dA", Px, "bra", 0)});
    EXPECT_FALSE(drv.is_electric_field(
        R2CTerm(T2CIntegral(OneCenterComponent("a", Px), OneCenterComponent("b", S),
                            OperatorComponent("AG", Px), 0, prefixes))));
}

TEST(T2CElectricFieldDriverTest, BraVrrOnPShell)
{
    const T2CElectricFieldDriver drv;

    // (P|S) with AG(1,0,0): PA term, PC term, and the operator-lowering term.
    const auto rec = drv.bra_vrr(field_term(Px, S, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 3u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PA", "PC"}));
}

TEST(T2CElectricFieldDriverTest, KetVrrOnPShell)
{
    const T2CElectricFieldDriver drv;

    // (S|P) with AG(1,0,0): PB term, PC term, and the operator-lowering term.
    const auto rec = drv.ket_vrr(field_term(S, Px, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 3u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PB", "PC"}));
}

TEST(T2CElectricFieldDriverTest, BraVrrRequiresBraMomentum)
{
    const T2CElectricFieldDriver drv;

    // Scalar bra: nothing to lower on the bra side -> nullopt.
    EXPECT_FALSE(drv.bra_vrr(field_term(S, S, Px), 'x').has_value());
}

TEST(T2CElectricFieldDriverTest, BraVrrRejectsNonFieldTerm)
{
    const T2CElectricFieldDriver drv;

    // Scalar AG operator is not an electric field term.
    EXPECT_FALSE(drv.bra_vrr(field_term(Px, S, S), 'x').has_value());
}

TEST(T2CElectricFieldDriverTest, CreateRecursionBuildsOneExpansionPerIntegral)
{
    const T2CElectricFieldDriver drv;

    const VT2CIntegrals vints({T2CIntegral(OneCenterComponent("a", Px),
                                           OneCenterComponent("b", S),
                                           OperatorComponent("AG", Px))});

    EXPECT_EQ(drv.create_recursion(vints).expansions(), 1u);
}
