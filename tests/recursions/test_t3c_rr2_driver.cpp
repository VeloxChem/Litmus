// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t3c_rr2_driver.hpp"
#include "t2c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

// Three-center r.r^2 term (operator "GR.R2(r)"); the operator shape carries the
// gradient component stepped down by aux_vrr.
R2CTerm rr2_term(const TensorComponent& bra, const TensorComponent& ket, const TensorComponent& op)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("GR.R2(r)", op)));
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

TEST(T3CRR2DriverTest, IsRR2)
{
    const T3CRR2Driver drv;

    EXPECT_TRUE(drv.is_rr2(rr2_term(S, S, Px)));

    const auto plain = R2CTerm(T2CIntegral(OneCenterComponent("a", S),
                                           OneCenterComponent("b", S),
                                           OperatorComponent("G(r)")));
    EXPECT_FALSE(drv.is_rr2(plain));
}

TEST(T3CRR2DriverTest, AuxVrrOnScalarShells)
{
    const T3CRR2Driver drv;

    // (S|S): the GR2(r) GC term and the G(r) GC+1/geta term.
    const auto rec = drv.aux_vrr(rr2_term(S, S, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"GC", "1/geta"}));
}

TEST(T3CRR2DriverTest, AuxVrrOnPBraAddsDownSteps)
{
    const T3CRR2Driver drv;

    const auto rec = drv.aux_vrr(rr2_term(Px, S, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 4u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"GC", "1/geta", "1/geta2"}));
}

TEST(T3CRR2DriverTest, AuxVrrRejectsNonRR2Term)
{
    const T3CRR2Driver drv;

    const auto plain = R2CTerm(T2CIntegral(OneCenterComponent("a", Px),
                                           OneCenterComponent("b", S),
                                           OperatorComponent("G(r)", Px)));

    EXPECT_FALSE(drv.aux_vrr(plain, 'x').has_value());
}
