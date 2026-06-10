// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t2c_eri_driver.hpp"
#include "t2c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);
const TensorComponent Dx(2, 0, 0);

R2CTerm eri_term(const TensorComponent& bra, const TensorComponent& ket)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("1/|r-r'|")));
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

}  // namespace

TEST(T2CElectronRepulsionDriverTest, IsElectronRepulsion)
{
    const T2CElectronRepulsionDriver drv;

    EXPECT_TRUE(drv.is_electron_repulsion(eri_term(Px, S)));

    const auto overlap = R2CTerm(T2CIntegral(OneCenterComponent("a", Px),
                                             OneCenterComponent("b", S),
                                             OperatorComponent("1")));
    EXPECT_FALSE(drv.is_electron_repulsion(overlap));
}

TEST(T2CElectronRepulsionDriverTest, BraVrrSingleTermOnPShell)
{
    const T2CElectronRepulsionDriver drv;

    // (P|S): only the PA / WP-style leading term (Boys-order incremented).
    const auto rec = drv.bra_vrr(eri_term(Px, S), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PA"}));
}

TEST(T2CElectronRepulsionDriverTest, BraVrrOnDShell)
{
    const T2CElectronRepulsionDriver drv;

    const auto rec = drv.bra_vrr(eri_term(Dx, S), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 3u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PA", "1/b_e", "zeta/b_e^2"}));
}

TEST(T2CElectronRepulsionDriverTest, BraVrrWithKetMomentumAddsCrossTerm)
{
    const T2CElectronRepulsionDriver drv;

    // (P|P): the PA term plus the 1/eta bra-ket cross term.
    const auto rec = drv.bra_vrr(eri_term(Px, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"PA", "1/eta"}));
}

TEST(T2CElectronRepulsionDriverTest, BraVrrRejectsNonEriTerm)
{
    const T2CElectronRepulsionDriver drv;

    const auto overlap = R2CTerm(T2CIntegral(OneCenterComponent("a", Px),
                                             OneCenterComponent("b", S),
                                             OperatorComponent("1")));

    EXPECT_FALSE(drv.bra_vrr(overlap, 'x').has_value());
}

TEST(T2CElectronRepulsionDriverTest, CreateRecursionBuildsOneExpansionPerIntegral)
{
    const T2CElectronRepulsionDriver drv;

    const VT2CIntegrals vints({T2CIntegral(OneCenterComponent("a", Px),
                                           OneCenterComponent("b", S),
                                           OperatorComponent("1/|r-r'|"))});

    EXPECT_EQ(drv.create_recursion(vints).expansions(), 1u);
}
