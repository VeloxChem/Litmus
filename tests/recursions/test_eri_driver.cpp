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

TEST_F(EriDriverTest, KetHrr)
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
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_xy_zz = T2CPair({"GC", "GD"}, {d_xy, d_zz});
    
    const auto tint = T4CIntegral(b_0_0, k_xy_zz, operi);
    
    const auto t4crec = R4CTerm(tint);
    
    // check recursion along x axis
    
    const auto k_y_xzz = T2CPair({"GC", "GD"}, {p_y, f_xzz});
    
    const auto r1aint = T4CIntegral(b_0_0, k_y_xzz, operi);
    
    const auto t1arec = R4CTerm(r1aint);
    
    const auto k_y_zz = T2CPair({"GC", "GD"}, {p_y, d_zz});
    
    const auto r2aint = T4CIntegral(b_0_0, k_y_zz, operi);
    
    const auto cdx = Factor("CD", "rcd", TensorComponent(1, 0, 0));
    
    const auto t2arec = R4CTerm(r2aint, {{cdx, 1}, }, Fraction(-1));
    
    EXPECT_EQ(eri_drv.ket_hrr(t4crec, 'x'), R4CDist(t4crec, {t1arec, t2arec}));
    
    // check recursion along y axis
    
    const auto k_x_yzz = T2CPair({"GC", "GD"}, {p_x, f_yzz});

    const auto r1bint = T4CIntegral(b_0_0, k_x_yzz, operi);
    
    const auto t1brec = R4CTerm(r1bint);

    const auto k_x_zz = T2CPair({"GC", "GD"}, {p_x, d_zz});

    const auto r2bint = T4CIntegral(b_0_0, k_x_zz, operi);
    
    const auto cdy = Factor("CD", "rcd", TensorComponent(0, 1, 0));
    
    const auto t2brec = R4CTerm(r2bint, {{cdy, 1}, }, Fraction(-1));
    
    EXPECT_EQ(eri_drv.ket_hrr(t4crec, 'y'), R4CDist(t4crec, {t1brec, t2brec}));
    
    // check recursion along z axis
    
    EXPECT_FALSE(eri_drv.ket_hrr(t4crec, 'z'));
}

TEST_F(EriDriverTest, BraVrr)
{
    EriDriver eri_drv;
    
    // recursion data
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);

    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);
    
    const auto f_xxx = TensorComponent(3, 0, 0);
    
    const auto f_xyy = TensorComponent(1, 2, 0);
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto b_0_xyy = T2CPair({"GA", "GB"}, {s_0, f_xyy});
    
    const auto k_0_xxx = T2CPair({"GC", "GD"}, {s_0, f_xxx});
    
    const auto tint = T4CIntegral(b_0_xyy, k_0_xxx, operi);
    
    const auto t4crec = R4CTerm(tint);
    
    // check recursion along x axis
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto r1aint = T4CIntegral(b_0_yy, k_0_xxx, operi);
    
    const auto pbx = Factor("PB", "rpb", TensorComponent(1, 0, 0));
    
    const auto t1arec = R4CTerm(r1aint, {{pbx, 1},}, Fraction(1));
    
    const auto r2aint = T4CIntegral(b_0_yy, k_0_xxx, operi, 1);
    
    const auto wpx = Factor("WP", "rwp", TensorComponent(1, 0, 0));
    
    const auto t2arec = R4CTerm(r2aint, {{wpx, 1},}, Fraction(1));
    
    const auto k_0_xx = T2CPair({"GC", "GD"}, {s_0, d_xx});
    
    const auto r3aint = T4CIntegral(b_0_yy, k_0_xx, operi, 1);
    
    const auto fze = Factor("1/(zeta+eta)", "fze", TensorComponent(0, 0, 0));
    
    const auto t3arec = R4CTerm(r3aint, {{fze, 1},}, Fraction(3, 2));
    
    EXPECT_EQ(eri_drv.bra_vrr(t4crec, 'x'), R4CDist(t4crec, {t1arec, t2arec, t3arec}));
    
    // check recursion along y axis
    
    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto r1bint = T4CIntegral(b_0_xy, k_0_xxx, operi);
    
    const auto pby = Factor("PB", "rpb", TensorComponent(0, 1, 0));
    
    const auto t1brec = R4CTerm(r1bint, {{pby, 1},}, Fraction(1));
    
    const auto r2bint = T4CIntegral(b_0_xy, k_0_xxx, operi, 1);
    
    const auto wpy = Factor("WP", "rwp", TensorComponent(0, 1, 0));
    
    const auto t2brec = R4CTerm(r2bint, {{wpy, 1},}, Fraction(1));
    
    const auto b_0_x = T2CPair({"GA", "GB"}, {s_0, p_x});
    
    const auto r3bint = T4CIntegral(b_0_x, k_0_xxx, operi);
    
    const auto fz = Factor("1/zeta", "fz", TensorComponent(0, 0, 0));
    
    const auto t3brec = R4CTerm(r3bint, {{fz, 1},}, Fraction(1, 2));
    
    const auto r4bint = T4CIntegral(b_0_x, k_0_xxx, operi, 1);
    
    const auto frz2 = Factor("rho/zeta^2", "frz2", TensorComponent(0, 0, 0));
    
    const auto t4brec = R4CTerm(r4bint, {{frz2, 1},}, Fraction(-1, 2));
    
    EXPECT_EQ(eri_drv.bra_vrr(t4crec, 'y'), R4CDist(t4crec, {t1brec, t2brec, t3brec, t4brec}));
    
    EXPECT_FALSE(eri_drv.bra_vrr(t4crec, 'z'));
}

TEST_F(EriDriverTest, KetVrr)
{
    EriDriver eri_drv;
    
    // recursion data
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);
    
    const auto f_xxx = TensorComponent(3, 0, 0);
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_0_xxx = T2CPair({"GC", "GD"}, {s_0, f_xxx});
    
    const auto tint = T4CIntegral(b_0_0, k_0_xxx, operi);
    
    const auto t4crec = R4CTerm(tint);
    
    // check recursion along x axis
    
    const auto k_0_xx = T2CPair({"GC", "GD"}, {s_0, d_xx});
    
    const auto r1aint = T4CIntegral(b_0_0, k_0_xx, operi);
    
    const auto qdx = Factor("QD", "rqd", TensorComponent(1, 0, 0));
    
    const auto t1arec = R4CTerm(r1aint, {{qdx, 1},}, Fraction(1));
    
    const auto r2aint = T4CIntegral(b_0_0, k_0_xx, operi, 1);
    
    const auto wqx = Factor("WQ", "rwq", TensorComponent(1, 0, 0));
    
    const auto t2arec = R4CTerm(r2aint, {{wqx, 1},}, Fraction(1));
    
    const auto k_0_x = T2CPair({"GC", "GD"}, {s_0, p_x});
    
    const auto r3aint = T4CIntegral(b_0_0, k_0_x, operi);
    
    const auto fe = Factor("1/eta", "fe", TensorComponent(0, 0, 0));
    
    const auto t3arec = R4CTerm(r3aint, {{fe, 1},}, Fraction(1));
    
    const auto r4aint = T4CIntegral(b_0_0, k_0_x, operi, 1);
    
    const auto fre2 = Factor("rho/eta^2", "fre2", TensorComponent(0, 0, 0));
    
    const auto t4arec = R4CTerm(r4aint, {{fre2, 1},}, Fraction(-1));
    
    EXPECT_EQ(eri_drv.ket_vrr(t4crec, 'x'), R4CDist(t4crec, {t1arec, t2arec, t3arec, t4arec}));
    
    EXPECT_FALSE(eri_drv.ket_vrr(t4crec, 'y'));
    
    EXPECT_FALSE(eri_drv.ket_vrr(t4crec, 'z'));
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

TEST_F(EriDriverTest, ApplyKetHrr)
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
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_xy_zz = T2CPair({"GC", "GD"}, {d_xy, d_zz});
    
    const auto tint = T4CIntegral(b_0_0, k_xy_zz, operi);
    
    const auto t4crec = R4CTerm(tint);
    
    // witout initial set of integrals
    
    const auto k_y_xzz = T2CPair({"GC", "GD"}, {p_y, f_xzz});
    
    const auto r1aint = T4CIntegral(b_0_0, k_y_xzz, operi);
    
    const auto t1arec = R4CTerm(r1aint);
    
    const auto k_y_zz = T2CPair({"GC", "GD"}, {p_y, d_zz});
    
    const auto r2aint = T4CIntegral(b_0_0, k_y_zz, operi);
    
    const auto cdx = Factor("CD", "rcd", TensorComponent(1, 0, 0));
    
    const auto t2arec = R4CTerm(r2aint, {{cdx, 1}, }, Fraction(-1));
    
    std::set<T4CIntegral> sints;
    
    auto r4cdist = eri_drv.apply_ket_hrr(t4crec, sints);
    
    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1arec, t2arec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1aint, r2aint}));
    
    // with initial set of integrals
    
    const auto k_x_yzz = T2CPair({"GC", "GD"}, {p_x, f_yzz});

    const auto r1bint = T4CIntegral(b_0_0, k_x_yzz, operi);
    
    const auto t1brec = R4CTerm(r1bint);

    const auto k_x_zz = T2CPair({"GC", "GD"}, {p_x, d_zz});

    const auto r2bint = T4CIntegral(b_0_0, k_x_zz, operi);
    
    const auto cdy = Factor("CD", "rcd", TensorComponent(0, 1, 0));
    
    const auto t2brec = R4CTerm(r2bint, {{cdy, 1}, }, Fraction(-1));
    
    sints = {r2bint,};
    
    r4cdist = eri_drv.apply_ket_hrr(t4crec, sints);

    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1brec, t2brec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1bint, r2bint}));
}

TEST_F(EriDriverTest, ApplyBraVrr)
{
    EriDriver eri_drv;
    
    // recursion data
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);

    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);
    
    const auto f_xxx = TensorComponent(3, 0, 0);
    
    const auto f_xyy = TensorComponent(1, 2, 0);
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto b_0_xyy = T2CPair({"GA", "GB"}, {s_0, f_xyy});
    
    const auto k_0_xxx = T2CPair({"GC", "GD"}, {s_0, f_xxx});
    
    const auto tint = T4CIntegral(b_0_xyy, k_0_xxx, operi);
    
    const auto t4crec = R4CTerm(tint);
    
    // witout initial set of integrals
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto r1aint = T4CIntegral(b_0_yy, k_0_xxx, operi);
    
    const auto pbx = Factor("PB", "rpb", TensorComponent(1, 0, 0));
    
    const auto t1arec = R4CTerm(r1aint, {{pbx, 1},}, Fraction(1));
    
    const auto r2aint = T4CIntegral(b_0_yy, k_0_xxx, operi, 1);
    
    const auto wpx = Factor("WP", "rwp", TensorComponent(1, 0, 0));
    
    const auto t2arec = R4CTerm(r2aint, {{wpx, 1},}, Fraction(1));
    
    const auto k_0_xx = T2CPair({"GC", "GD"}, {s_0, d_xx});
    
    const auto r3aint = T4CIntegral(b_0_yy, k_0_xx, operi, 1);
    
    const auto fze = Factor("1/(zeta+eta)", "fze", TensorComponent(0, 0, 0));
    
    const auto t3arec = R4CTerm(r3aint, {{fze, 1},}, Fraction(3, 2));
    
    std::set<T4CIntegral> sints;
    
    auto r4cdist = eri_drv.apply_bra_vrr(t4crec, sints);
    
    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1arec, t2arec, t3arec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1aint, r2aint, r3aint}));
    
    // with initial set of integrals
    
    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto r1bint = T4CIntegral(b_0_xy, k_0_xxx, operi);
    
    const auto pby = Factor("PB", "rpb", TensorComponent(0, 1, 0));
    
    const auto t1brec = R4CTerm(r1bint, {{pby, 1},}, Fraction(1));
    
    const auto r2bint = T4CIntegral(b_0_xy, k_0_xxx, operi, 1);
    
    const auto wpy = Factor("WP", "rwp", TensorComponent(0, 1, 0));
    
    const auto t2brec = R4CTerm(r2bint, {{wpy, 1},}, Fraction(1));
    
    const auto b_0_x = T2CPair({"GA", "GB"}, {s_0, p_x});
    
    const auto r3bint = T4CIntegral(b_0_x, k_0_xxx, operi);
    
    const auto fz = Factor("1/zeta", "fz", TensorComponent(0, 0, 0));
    
    const auto t3brec = R4CTerm(r3bint, {{fz, 1},}, Fraction(1, 2));
    
    const auto r4bint = T4CIntegral(b_0_x, k_0_xxx, operi, 1);
    
    const auto frz2 = Factor("rho/zeta^2", "frz2", TensorComponent(0, 0, 0));
    
    const auto t4brec = R4CTerm(r4bint, {{frz2, 1},}, Fraction(-1, 2));
    
    sints = {r1bint, r2bint, r4bint};
    
    r4cdist = eri_drv.apply_bra_vrr(t4crec, sints);

    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1brec, t2brec, t3brec, t4brec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1bint, r2bint, r3bint, r4bint}));
}

TEST_F(EriDriverTest, AppyKetVrr)
{
    EriDriver eri_drv;
    
    // recursion data
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);
    
    const auto f_xxx = TensorComponent(3, 0, 0);
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_0_xxx = T2CPair({"GC", "GD"}, {s_0, f_xxx});
    
    const auto tint = T4CIntegral(b_0_0, k_0_xxx, operi);
    
    const auto t4crec = R4CTerm(tint);
    
    // witout initial set of integrals
    
    const auto k_0_xx = T2CPair({"GC", "GD"}, {s_0, d_xx});
    
    const auto r1aint = T4CIntegral(b_0_0, k_0_xx, operi);
    
    const auto qdx = Factor("QD", "rqd", TensorComponent(1, 0, 0));
    
    const auto t1arec = R4CTerm(r1aint, {{qdx, 1},}, Fraction(1));
    
    const auto r2aint = T4CIntegral(b_0_0, k_0_xx, operi, 1);
    
    const auto wqx = Factor("WQ", "rwq", TensorComponent(1, 0, 0));
    
    const auto t2arec = R4CTerm(r2aint, {{wqx, 1},}, Fraction(1));
    
    const auto k_0_x = T2CPair({"GC", "GD"}, {s_0, p_x});
    
    const auto r3aint = T4CIntegral(b_0_0, k_0_x, operi);
    
    const auto fe = Factor("1/eta", "fe", TensorComponent(0, 0, 0));
    
    const auto t3arec = R4CTerm(r3aint, {{fe, 1},}, Fraction(1));
    
    const auto r4aint = T4CIntegral(b_0_0, k_0_x, operi, 1);
    
    const auto fre2 = Factor("rho/eta^2", "fre2", TensorComponent(0, 0, 0));
    
    const auto t4arec = R4CTerm(r4aint, {{fre2, 1},}, Fraction(-1));
    
    std::set<T4CIntegral> sints;
    
    auto r4cdist = eri_drv.apply_ket_vrr(t4crec, sints);
    
    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1arec, t2arec, t3arec, t4arec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1aint, r2aint, r3aint, r4aint}));
    
    // with initial set of integrals
    
    sints = {r3aint, r4aint};
    
    r4cdist = eri_drv.apply_ket_vrr(t4crec, sints);

    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1arec, t2arec, t3arec, t4arec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1aint, r2aint, r3aint, r4aint}));
}

TEST_F(EriDriverTest, ApplyBraHrrForGroup)
{
    EriDriver eri_drv;
    
    // recursion data
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto b_x_x = T2CPair({"GA", "GB"}, {p_x, p_x});
    
    const auto b_y_x = T2CPair({"GA", "GB"}, {p_y, p_x});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    const auto taint = T4CIntegral(b_x_x, k_0_0, operi);
    
    const auto tbint = T4CIntegral(b_y_x, k_0_0, operi);
    
    const auto t4arec = R4CTerm(taint);
    
    const auto t4brec = R4CTerm(tbint);
    
    // generate recursion group
    
    std::set<T4CIntegral> sints;
    
    const auto t4g = eri_drv.apply_bra_hrr({t4arec, t4brec}, sints);
    
    // witout initial set of integrals
    
    std::set<T4CIntegral> rints;
    
    const auto r4adist = eri_drv.apply_bra_hrr(t4arec, rints);
    
    const auto r4bdist = eri_drv.apply_bra_hrr(t4brec, rints);
    
    EXPECT_EQ(t4g, R4Group({r4adist, r4bdist}));
    
    EXPECT_EQ(sints, rints);
}

TEST_F(EriDriverTest, ApplyKetHrrForGroup)
{
    EriDriver eri_drv;
    
    // recursion data
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto operi = OperatorComponent("1/|r-r'|");
   
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_x_x = T2CPair({"GC", "GD"}, {p_x, p_x});
    
    const auto k_y_x = T2CPair({"GC", "GD"}, {p_y, p_x});
    
    const auto taint = T4CIntegral(b_0_0, k_x_x, operi);
    
    const auto tbint = T4CIntegral(b_0_0, k_y_x, operi);
    
    const auto t4arec = R4CTerm(taint);
    
    const auto t4brec = R4CTerm(tbint);
    
    // generate recursion group
    
    std::set<T4CIntegral> sints;
    
    const auto t4g = eri_drv.apply_ket_hrr({t4arec, t4brec}, sints);
    
    // witout initial set of integrals
    
    std::set<T4CIntegral> rints;
    
    const auto r4adist = eri_drv.apply_ket_hrr(t4arec, rints);
    
    const auto r4bdist = eri_drv.apply_ket_hrr(t4brec, rints);
    
    EXPECT_EQ(t4g, R4Group({r4adist, r4bdist}));
    
    EXPECT_EQ(sints, rints);
}

TEST_F(EriDriverTest, ApplyBraVrrForGroup)
{
    EriDriver eri_drv;
    
    // recursion data
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);
    
    const auto f_xyy = TensorComponent(1, 3, 0);
    
    const auto f_xyz = TensorComponent(1, 1, 1);
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto k_0_xyy = T2CPair({"GC", "GD"}, {s_0, f_xyy});
                                       
    const auto k_0_xyz = T2CPair({"GC", "GD"}, {s_0, f_xyz});
                                       
    const auto taint = T4CIntegral(b_0_xy, k_0_xyy, operi);
    
    const auto tbint = T4CIntegral(b_0_yy, k_0_xyz, operi);
    
    const auto t4arec = R4CTerm(taint);
    
    const auto t4brec = R4CTerm(tbint);
    
    // generate recursion group
    
    std::set<T4CIntegral> sints;
    
    const auto t4g = eri_drv.apply_bra_vrr({t4arec, t4brec}, sints);
    
    // witout initial set of integrals
    
    std::set<T4CIntegral> rints;
    
    const auto r4adist = eri_drv.apply_bra_vrr(t4arec, rints);
    
    const auto r4bdist = eri_drv.apply_bra_vrr(t4brec, rints);
    
    EXPECT_EQ(t4g, R4Group({r4adist, r4bdist}));
    
    EXPECT_EQ(sints, rints);
}

TEST_F(EriDriverTest, ApplyKetVrrForGroup)
{
    EriDriver eri_drv;
    
    // recursion data
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto f_xyy = TensorComponent(1, 3, 0);
    
    const auto f_xyz = TensorComponent(1, 1, 1);
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_0_xyy = T2CPair({"GC", "GD"}, {s_0, f_xyy});
                                       
    const auto k_0_xyz = T2CPair({"GC", "GD"}, {s_0, f_xyz});
                                       
    const auto taint = T4CIntegral(b_0_0, k_0_xyy, operi);
    
    const auto tbint = T4CIntegral(b_0_0, k_0_xyz, operi);
    
    const auto t4arec = R4CTerm(taint);
    
    const auto t4brec = R4CTerm(tbint);
    
    // generate recursion group
    
    std::set<T4CIntegral> sints;
    
    const auto t4g = eri_drv.apply_ket_vrr({t4arec, t4brec}, sints);
    
    // witout initial set of integrals
    
    std::set<T4CIntegral> rints;
    
    const auto r4adist = eri_drv.apply_ket_vrr(t4arec, rints);
    
    const auto r4bdist = eri_drv.apply_ket_vrr(t4brec, rints);
    
    EXPECT_EQ(t4g, R4Group({r4adist, r4bdist}));
    
    EXPECT_EQ(sints, rints);
}
