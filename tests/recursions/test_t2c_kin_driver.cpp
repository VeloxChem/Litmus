// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t2c_kin_driver.hpp"
#include "t2c_defs.hpp"

namespace {

// Kinetic-energy recursion term (operator "T") with the given bra/ket shapes.
R2CTerm kin_term(const TensorComponent& bra, const TensorComponent& ket)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("T")));
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

// Integrand operator names across an expansion's terms.
std::set<std::string> term_operators(const R2CDist& dist)
{
    std::set<std::string> names;
    for (size_t i = 0; i < dist.terms(); i++)
    {
        names.insert(dist[i].integral().integrand().name());
    }
    return names;
}

}  // namespace

TEST(T2CKineticEnergyDriverTest, IsKineticEnergy)
{
    const T2CKineticEnergyDriver drv;

    EXPECT_TRUE(drv.is_kinetic_energy(kin_term(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0))));

    const auto overlap = R2CTerm(T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                             OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                             OperatorComponent("1")));
    EXPECT_FALSE(drv.is_kinetic_energy(overlap));
}

TEST(T2CKineticEnergyDriverTest, BraVrrOnPShell)
{
    const T2CKineticEnergyDriver drv;

    // (P|S): the 2*zeta overlap term and the PA kinetic term.
    const auto rec = drv.bra_vrr(kin_term(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"zeta", "PA"}));
    // The kinetic recursion spawns an overlap ("1") term alongside the "T" term.
    EXPECT_EQ(term_operators(*rec), std::set<std::string>({"1", "T"}));
}

TEST(T2CKineticEnergyDriverTest, BraVrrOnDShell)
{
    const T2CKineticEnergyDriver drv;

    const auto rec = drv.bra_vrr(kin_term(TensorComponent(2, 0, 0), TensorComponent(0, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 4u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"zeta", "PA", "1/eta", "1/b_e"}));
}

TEST(T2CKineticEnergyDriverTest, KetVrrOnPShell)
{
    const T2CKineticEnergyDriver drv;

    const auto rec = drv.ket_vrr(kin_term(TensorComponent(0, 0, 0), TensorComponent(1, 0, 0)), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"zeta", "PB"}));
}

TEST(T2CKineticEnergyDriverTest, VrrRejectsNonKineticTerm)
{
    const T2CKineticEnergyDriver drv;

    const auto overlap = R2CTerm(T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                             OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                             OperatorComponent("1")));

    EXPECT_FALSE(drv.bra_vrr(overlap, 'x').has_value());
}

TEST(T2CKineticEnergyDriverTest, CreateRecursionBuildsOneExpansionPerIntegral)
{
    const T2CKineticEnergyDriver drv;

    const VT2CIntegrals vints({T2CIntegral(OneCenterComponent("a", TensorComponent(1, 0, 0)),
                                           OneCenterComponent("b", TensorComponent(0, 0, 0)),
                                           OperatorComponent("T"))});

    EXPECT_EQ(drv.create_recursion(vints).expansions(), 1u);
}
