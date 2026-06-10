// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <set>
#include <string>

#include "t2c_hrr_driver.hpp"
#include "t2c_defs.hpp"

namespace {

const TensorComponent S(0, 0, 0);
const TensorComponent Px(1, 0, 0);

R2CTerm hrr_term(const TensorComponent& bra, const TensorComponent& ket)
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

TEST(T2CHRRDriverTest, BraVrrTransfersBraToKet)
{
    const T2CHRRDriver drv;

    // (P|S): lower the bra, raise the ket, with the AB distance factor.
    const auto rec = drv.bra_vrr(hrr_term(Px, S), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"AB"}));
}

TEST(T2CHRRDriverTest, KetVrrTransfersKetToBra)
{
    const T2CHRRDriver drv;

    const auto rec = drv.ket_vrr(hrr_term(S, Px), 'x');

    ASSERT_TRUE(rec.has_value());
    EXPECT_EQ(rec->terms(), 2u);
    EXPECT_EQ(factor_names(*rec), std::set<std::string>({"AB"}));
}

TEST(T2CHRRDriverTest, BraVrrRequiresBraMomentum)
{
    const T2CHRRDriver drv;

    EXPECT_FALSE(drv.bra_vrr(hrr_term(S, S), 'x').has_value());
}
