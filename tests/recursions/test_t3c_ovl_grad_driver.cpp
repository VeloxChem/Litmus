// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t3c_ovl_grad_driver.hpp"
#include "t2c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

// Three-center overlap-gradient term (operator "GX(r)"); the operator shape
// carries the gradient component.
R2CTerm grad_term(const TensorComponent& bra, const TensorComponent& ket, const TensorComponent& op)
{
    return R2CTerm(T2CIntegral(OneCenterComponent("a", bra),
                               OneCenterComponent("b", ket),
                               OperatorComponent("GX(r)", op)));
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

TEST(T3COverlapGradientDriverTest, IsOverlapGradient)
{
    const T3COverlapGradientDriver drv;

    EXPECT_TRUE(drv.is_overlap_gradient(grad_term(S, S, Px)));

    const auto plain = R2CTerm(T2CIntegral(OneCenterComponent("a", S),
                                           OneCenterComponent("b", S),
                                           OperatorComponent("G(r)")));
    EXPECT_FALSE(drv.is_overlap_gradient(plain));
}

TEST(T3COverlapGradientDriverTest, AuxVrrOnScalarShells)
{
    const T3COverlapGradientDriver drv;

    // (S|S): the operator down-step (GX(r) -> G(r)) GC term only.
    const auto rec = drv.aux_vrr(grad_term(S, S, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 1u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"GC", "c_e"}));
}

TEST(T3COverlapGradientDriverTest, AuxVrrOnPBraAddsDownStep)
{
    const T3COverlapGradientDriver drv;

    const auto rec = drv.aux_vrr(grad_term(Px, S, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"GC", "c_e", "1/geta"}));
}

TEST(T3COverlapGradientDriverTest, AuxVrrRejectsNonGradientTerm)
{
    const T3COverlapGradientDriver drv;

    const auto plain = R2CTerm(T2CIntegral(OneCenterComponent("a", Px),
                                           OneCenterComponent("b", S),
                                           OperatorComponent("G(r)", Px)));

    EXPECT_FALSE(drv.aux_vrr(plain, 'x').has_value());
}
