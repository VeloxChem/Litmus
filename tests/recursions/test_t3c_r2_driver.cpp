// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t3c_r2_driver.hpp"
#include "t2c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

// Three-center r^2 term (operator "GR2(r)").
R2CTerm r2_term(const TensorComponent& bra, const TensorComponent& ket)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("GR2(r)")));
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

TEST(T3CR2DriverTest, IsR2)
{
    const T3CR2Driver drv;

    EXPECT_TRUE(drv.is_r2(r2_term(S, S)));

    const auto plain = R2CTerm(T2CIntegral(OneCenterComponent("a", S),
                                           OneCenterComponent("b", S),
                                           OperatorComponent("G(r)")));
    EXPECT_FALSE(drv.is_r2(plain));
}

TEST(T3CR2DriverTest, AuxVrrOnScalarShells)
{
    const T3CR2Driver drv;

    // (S|S): the r2gc contribution plus the operator-delta (1/geta) term; no
    // bra/ket down-steps are available.
    const auto rec = drv.aux_vrr(r2_term(S, S));

    EXPECT_EQ(rec.terms(), 2u);
    EXPECT_EQ(factor_names(rec), std::set<std::string>({"r2gc", "1/geta"}));
}

TEST(T3CR2DriverTest, AuxVrrOnPBraAddsDownStep)
{
    const T3CR2Driver drv;

    // Adds the bra-down GC/1/geta contribution on top of the scalar terms.
    const auto rec = drv.aux_vrr(r2_term(Px, S));

    EXPECT_EQ(rec.terms(), 3u);
    EXPECT_EQ(factor_names(rec), std::set<std::string>({"r2gc", "GC", "1/geta"}));
}

TEST(T3CR2DriverTest, AuxVrrOnNonR2HasNoTerms)
{
    const T3CR2Driver drv;

    const auto plain = R2CTerm(T2CIntegral(OneCenterComponent("a", Px),
                                           OneCenterComponent("b", S),
                                           OperatorComponent("G(r)")));

    EXPECT_EQ(drv.aux_vrr(plain).terms(), 0u);
}
