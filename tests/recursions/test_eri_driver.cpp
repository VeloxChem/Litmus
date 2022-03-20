// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.
// E-mail: rinkevic@kth.se
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "test_eri_driver.hpp"

#include "eri_driver.hpp"
#include "recursion_expansion.hpp"
#include "recursion_term.hpp"
#include "integral_component.hpp"
#include "two_center_pair_component.hpp"

using T2CPair = TwoCenterPairComponent;

using T4CIntegral = IntegralComponent<T2CPair, T2CPair>;

using R4CTerm = RecursionTerm<T4CIntegral>;

using R4CDist = RecursionExpansion<T4CIntegral>;

TEST_F(EriDriverTest, BraHrr)
{
    EriDriver eri_drv;
    
    // recursion data
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_zz = TensorComponent(0, 0, 2);
    
    const auto f_xzz = TensorComponent(1, 0, 2);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto b_xy_zz = T2CPair({"GA", "GB"}, {d_xy, d_zz});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    const auto tint = T4CIntegral(b_xy_zz, k_0_0, operi);
    
    const auto t4crec = R4CTerm(tint);
    
    // check recursion along x axis
    
    const auto b_y_xzz = T2CPair({"GA", "GB"}, {p_y, f_xzz});
    
    const auto r1aint = T4CIntegral(b_y_xzz, k_0_0, operi);
    
    const auto t1arec = R4CTerm(r1aint);
    
    const auto b_y_zz = T2CPair({"GA", "GB"}, {p_y, d_zz});
    
    const auto r2aint = T4CIntegral(b_y_zz, k_0_0, operi);
    
    const auto abx = Factor("AB", "rab", TensorComponent(1, 0, 0));
    
    const auto t2arec = R4CTerm(r2aint, {{abx, 1}, }, Fraction(-1));
    
    EXPECT_EQ(eri_drv.bra_hrr(t4crec, 'x'), R4CDist(t4crec, {t1arec, t2arec}));
    
    // check recursion along y axis
    
    const auto b_x_yzz = T2CPair({"GA", "GB"}, {p_x, f_yzz});

    const auto r1bint = T4CIntegral(b_x_yzz, k_0_0, operi);
    
    const auto t1brec = R4CTerm(r1bint);

    const auto b_x_zz = T2CPair({"GA", "GB"}, {p_x, d_zz});

    const auto r2bint = T4CIntegral(b_x_zz, k_0_0, operi);
    
    const auto aby = Factor("AB", "rab", TensorComponent(0, 1, 0));
    
    const auto t2brec = R4CTerm(r2bint, {{aby, 1}, }, Fraction(-1));
    
    EXPECT_EQ(eri_drv.bra_hrr(t4crec, 'y'), R4CDist(t4crec, {t1brec, t2brec}));
    
    // check recursion along z axis
    
    EXPECT_FALSE(eri_drv.bra_hrr(t4crec, 'z'));
}

TEST_F(EriDriverTest, ApplyBraHrr)
{
    EriDriver eri_drv;
    
    // recursion data
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_zz = TensorComponent(0, 0, 2);
    
    const auto f_xzz = TensorComponent(1, 0, 2);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto b_xy_zz = T2CPair({"GA", "GB"}, {d_xy, d_zz});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    const auto tint = T4CIntegral(b_xy_zz, k_0_0, operi);
    
    const auto t4crec = R4CTerm(tint);
    
    // witout initial set of integrals
    
    const auto b_y_xzz = T2CPair({"GA", "GB"}, {p_y, f_xzz});
    
    const auto r1aint = T4CIntegral(b_y_xzz, k_0_0, operi);
    
    const auto t1arec = R4CTerm(r1aint);
    
    const auto b_y_zz = T2CPair({"GA", "GB"}, {p_y, d_zz});
    
    const auto r2aint = T4CIntegral(b_y_zz, k_0_0, operi);
    
    const auto abx = Factor("AB", "rab", TensorComponent(1, 0, 0));
    
    const auto t2arec = R4CTerm(r2aint, {{abx, 1}, }, Fraction(-1));
    
    std::set<T4CIntegral> sints;
    
    auto r4cdist = eri_drv.apply_bra_hrr(t4crec, sints);
    
    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1arec, t2arec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1aint, r2aint}));
    
    // with initial set of integrals
    
    const auto b_x_yzz = T2CPair({"GA", "GB"}, {p_x, f_yzz});

    const auto r1bint = T4CIntegral(b_x_yzz, k_0_0, operi);
    
    const auto t1brec = R4CTerm(r1bint);

    const auto b_x_zz = T2CPair({"GA", "GB"}, {p_x, d_zz});

    const auto r2bint = T4CIntegral(b_x_zz, k_0_0, operi);
    
    const auto aby = Factor("AB", "rab", TensorComponent(0, 1, 0));
    
    const auto t2brec = R4CTerm(r2bint, {{aby, 1}, }, Fraction(-1));
    
    sints = {r2bint,};
    
    r4cdist = eri_drv.apply_bra_hrr(t4crec, sints);

    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1brec, t2brec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1bint, r2bint}));
}
