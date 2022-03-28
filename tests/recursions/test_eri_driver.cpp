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

using R4Group = RecursionGroup<T4CIntegral>;

using R4Graph = Graph<R4Group>;

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

TEST_F(EriDriverTest, ApplyBraHrrWithGraphForPP)
{
    EriDriver eri_drv;
    
    // recursion data
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto b_x_x = T2CPair({"GA", "GB"}, {p_x, p_x});
    
    const auto b_x_y = T2CPair({"GA", "GB"}, {p_x, p_y});
    
    const auto b_y_y = T2CPair({"GA", "GB"}, {p_y, p_y});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    const auto t_x_x = T4CIntegral(b_x_x, k_0_0, operi);
    
    const auto t_x_y = T4CIntegral(b_x_y, k_0_0, operi);
    
    const auto t_y_y = T4CIntegral(b_y_y, k_0_0, operi);
    
    const auto rd_x_x =  R4CDist(R4CTerm(t_x_x));
    
    const auto rd_x_y =  R4CDist(R4CTerm(t_x_y));
    
    const auto rd_y_y =  R4CDist(R4CTerm(t_y_y));
    
    // initialize graph
    
    R4Graph rgraph(R4Group({rd_x_x, rd_x_y, rd_y_y}));
    
    std::set<T4CIntegral> sints;
    
    // apply horizontal bra recursion
    
    eri_drv.apply_bra_hrr(rgraph, sints);
    
    // set up reference data:
    
    const auto b_0_xx = T2CPair({"GA", "GB"}, {s_0, d_xx});
    
    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto t_0_xx = T4CIntegral(b_0_xx, k_0_0, operi);
    
    const auto t_0_xy = T4CIntegral(b_0_xy, k_0_0, operi);
    
    const auto t_0_yy = T4CIntegral(b_0_yy, k_0_0, operi);
    
    const auto rd_0_xx = R4CDist(R4CTerm(t_0_xx));
    
    const auto rd_0_xy = R4CDist(R4CTerm(t_0_xy));
    
    const auto rd_0_yy = R4CDist(R4CTerm(t_0_yy));
    
    const auto b_0_x = T2CPair({"GA", "GB"}, {s_0, p_x});
    
    const auto b_0_y = T2CPair({"GA", "GB"}, {s_0, p_y});
    
    const auto t_0_x = T4CIntegral(b_0_x, k_0_0, operi);
    
    const auto t_0_y = T4CIntegral(b_0_y, k_0_0, operi);
    
    const auto rd_0_x = R4CDist(R4CTerm(t_0_x));
    
    const auto rd_0_y = R4CDist(R4CTerm(t_0_y));
    
    std::set<T4CIntegral> rints;
    
    const auto hrr_x_x = eri_drv.apply_bra_hrr(R4CTerm(t_x_x), rints);
    
    const auto hrr_x_y = eri_drv.apply_bra_hrr(R4CTerm(t_x_y), rints);
    
    const auto hrr_y_y = eri_drv.apply_bra_hrr(R4CTerm(t_y_y), rints);
    
    // compare by terms
    
    EXPECT_EQ(rgraph.vertices(), 3);
    
    EXPECT_EQ(rgraph[0], R4Group({hrr_x_x, hrr_x_y, hrr_y_y}));
    
    EXPECT_EQ(rgraph[1], R4Group({rd_0_y, rd_0_x}));
    
    EXPECT_EQ(rgraph[2], R4Group({rd_0_yy, rd_0_xy, rd_0_xx}));
    
    EXPECT_EQ(sints, rints);
    
    EXPECT_EQ(rgraph.edge(0), std::set<int>({1,2}));
    
    EXPECT_EQ(rgraph.edge(1), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(2), std::set<int>({}));
}

#include <iostream>

TEST_F(EriDriverTest, ApplyBraHrrWithGraphForDD)
{
    EriDriver eri_drv;
    
    // recursion data
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);
    
    const auto f_xxx = TensorComponent(3, 0, 0);
    
    const auto f_xxy = TensorComponent(2, 1, 0);
    
    const auto f_xyy = TensorComponent(1, 2, 0);
    
    const auto f_yyy = TensorComponent(0, 3, 0);
    
    const auto g_xxxx = TensorComponent(4, 0, 0);
    
    const auto g_xxyy = TensorComponent(2, 2, 0);
    
    const auto g_yyyy = TensorComponent(0, 4, 0);
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    const auto b_xx_xx = T2CPair({"GA", "GB"}, {d_xx, d_xx});
    
    const auto b_xy_xy = T2CPair({"GA", "GB"}, {d_xy, d_xy});
    
    const auto b_yy_yy = T2CPair({"GA", "GB"}, {d_yy, d_yy});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    const auto t_xx_xx = T4CIntegral(b_xx_xx, k_0_0, operi);
    
    const auto t_xy_xy = T4CIntegral(b_xy_xy, k_0_0, operi);
    
    const auto t_yy_yy = T4CIntegral(b_yy_yy, k_0_0, operi);
    
    auto rd_xx_xx =  R4CDist(R4CTerm(t_xx_xx));
    
    auto rd_xy_xy =  R4CDist(R4CTerm(t_xy_xy));
    
    auto rd_yy_yy =  R4CDist(R4CTerm(t_yy_yy));
    
    // initialize graph
    
    R4Graph rgraph(R4Group({rd_xx_xx, rd_xy_xy, rd_yy_yy}));
    
    std::set<T4CIntegral> sints;
    
    // apply horizontal bra recursion
    
    eri_drv.apply_bra_hrr(rgraph, sints);
    
    // check number of vertices
    
    EXPECT_EQ(rgraph.vertices(), 6);
    
    // reference (pd|ss) integrals
    
    const auto b_x_xx = T2CPair({"GA", "GB"}, {p_x, d_xx});
    
    const auto b_y_xy = T2CPair({"GA", "GB"}, {p_y, d_xy});
    
    const auto b_y_yy = T2CPair({"GA", "GB"}, {p_y, d_yy});
    
    const auto t_x_xx = T4CIntegral(b_x_xx, k_0_0, operi);
    
    const auto t_y_xy = T4CIntegral(b_y_xy, k_0_0, operi);
    
    const auto t_y_yy = T4CIntegral(b_y_yy, k_0_0, operi);
    
    const auto rabx = Factor("AB", "rab", TensorComponent(1, 0, 0));
    
    const auto raby = Factor("AB", "rab", TensorComponent(0, 1, 0));
    
    const auto rt_x_xx = R4CTerm(t_x_xx, {{rabx, 1},}, Fraction(-1));
    
    const auto rt_y_xy = R4CTerm(t_y_xy, {{rabx, 1},}, Fraction(-1));
    
    const auto rt_y_yy = R4CTerm(t_y_yy, {{raby, 1},}, Fraction(-1));
    
    // reference (pf|ss) integrals
    
    const auto b_x_xxx = T2CPair({"GA", "GB"}, {p_x, f_xxx});
    
    const auto b_y_xxy = T2CPair({"GA", "GB"}, {p_y, f_xxy});
    
    const auto b_y_yyy = T2CPair({"GA", "GB"}, {p_y, f_yyy});
    
    const auto t_x_xxx = T4CIntegral(b_x_xxx, k_0_0, operi);
    
    const auto t_y_xxy = T4CIntegral(b_y_xxy, k_0_0, operi);
    
    const auto t_y_yyy = T4CIntegral(b_y_yyy, k_0_0, operi);
    
    const auto rt_x_xxx = R4CTerm(t_x_xxx);
    
    const auto rt_y_xxy = R4CTerm(t_y_xxy);
    
    const auto rt_y_yyy = R4CTerm(t_y_yyy);
    
    // check first recursion step
    
    rd_xx_xx =  R4CDist(R4CTerm(t_xx_xx), {rt_x_xxx, rt_x_xx,});
    
    rd_xy_xy =  R4CDist(R4CTerm(t_xy_xy), {rt_y_xxy, rt_y_xy});
    
    rd_yy_yy =  R4CDist(R4CTerm(t_yy_yy), {rt_y_yyy, rt_y_yy});
    
    EXPECT_EQ(rgraph[0], R4Group({rd_xx_xx, rd_xy_xy, rd_yy_yy,}));
    
    // reference (sd|ss) integrals
    
    const auto b_0_xx = T2CPair({"GA", "GB"}, {s_0, d_xx});
    
    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto t_0_xx = T4CIntegral(b_0_xx, k_0_0, operi);
    
    const auto t_0_xy = T4CIntegral(b_0_xy, k_0_0, operi);
    
    const auto t_0_yy = T4CIntegral(b_0_yy, k_0_0, operi);
    
    const auto rt_0_xx = R4CTerm(t_0_xx, {{rabx, 1},}, Fraction(-1));
    
    const auto rt_0_xy = R4CTerm(t_0_xy, {{raby, 1},}, Fraction(-1));
    
    const auto rt_0_yy = R4CTerm(t_0_yy, {{raby, 1},}, Fraction(-1));
    
    // reference (sf|ss) integrals
    
    const auto b_0_xxx = T2CPair({"GA", "GB"}, {s_0, f_xxx});
    
    const auto b_0_xxy = T2CPair({"GA", "GB"}, {s_0, f_xxy});
    
    const auto b_0_xyy = T2CPair({"GA", "GB"}, {s_0, f_xyy});
    
    const auto b_0_yyy = T2CPair({"GA", "GB"}, {s_0, f_yyy});
    
    const auto t_0_xxx = T4CIntegral(b_0_xxx, k_0_0, operi);
    
    const auto t_0_xxy = T4CIntegral(b_0_xxy, k_0_0, operi);
    
    const auto t_0_xyy = T4CIntegral(b_0_xyy, k_0_0, operi);
    
    const auto t_0_yyy = T4CIntegral(b_0_yyy, k_0_0, operi);
    
    auto rt_0_xxx = R4CTerm(t_0_xxx);
    
    auto rt_0_xxy = R4CTerm(t_0_xxy);
    
    auto rt_0_xyy = R4CTerm(t_0_xyy);
   
    auto rt_0_yyy = R4CTerm(t_0_yyy);
    
    // check second step in recursion
    
    auto rd_x_xx =  R4CDist(R4CTerm(t_x_xx), {rt_0_xxx, rt_0_xx,});
    
    auto rd_y_xy =  R4CDist(R4CTerm(t_y_xy), {rt_0_xyy, rt_0_xy});
    
    auto rd_y_yy =  R4CDist(R4CTerm(t_y_yy), {rt_0_yyy, rt_0_yy});
    
    EXPECT_EQ(rgraph[1], R4Group({rd_x_xx, rd_y_xy, rd_y_yy,}));
    
    // reference (sg|ss) integrals
    
    const auto b_0_xxxx = T2CPair({"GA", "GB"}, {s_0, g_xxxx});
    
    const auto b_0_xxyy = T2CPair({"GA", "GB"}, {s_0, g_xxyy});
    
    const auto b_0_yyyy = T2CPair({"GA", "GB"}, {s_0, g_yyyy});
    
    const auto t_0_xxxx = T4CIntegral(b_0_xxxx, k_0_0, operi);
    
    const auto t_0_xxyy = T4CIntegral(b_0_xxyy, k_0_0, operi);
    
    const auto t_0_yyyy = T4CIntegral(b_0_yyyy, k_0_0, operi);
    
    auto rt_0_xxxx = R4CTerm(t_0_xxxx);
    
    auto rt_0_xxyy = R4CTerm(t_0_xxyy);
   
    auto rt_0_yyyy = R4CTerm(t_0_yyyy);
    
    // update recursion (sf|ss) terms
    
    rt_0_yyy.add(raby); rt_0_yyy.scale(Fraction(-1));
    
    rt_0_xxy.add(raby); rt_0_xxy.scale(Fraction(-1));
    
    rt_0_xxx.add(rabx); rt_0_xxx.scale(Fraction(-1));
    
    // check third recursion term
    
    auto rd_x_xxx =  R4CDist(R4CTerm(t_x_xxx), {rt_0_xxxx, rt_0_xxx,});
    
    auto rd_y_xxy =  R4CDist(R4CTerm(t_y_xxy), {rt_0_xxyy, rt_0_xxy});
    
    auto rd_y_yyy =  R4CDist(R4CTerm(t_y_yyy), {rt_0_yyyy, rt_0_yyy});
    
    EXPECT_EQ(rgraph[2], R4Group({rd_x_xxx, rd_y_xxy, rd_y_yyy,}));
    
    // check fourth recursion term
    
    auto rd_0_xx =  R4CDist(R4CTerm(t_0_xx));
    
    auto rd_0_xy =  R4CDist(R4CTerm(t_0_xy));
    
    auto rd_0_yy =  R4CDist(R4CTerm(t_0_yy));
    
    EXPECT_EQ(rgraph[3], R4Group({rd_0_xx, rd_0_xy, rd_0_yy,}));
    
    // check fifth recursion term
    
    auto rd_0_xxx =  R4CDist(R4CTerm(t_0_xxx));
    
    auto rd_0_xxy =  R4CDist(R4CTerm(t_0_xxy));
    
    auto rd_0_xyy =  R4CDist(R4CTerm(t_0_xyy));
    
    auto rd_0_yyy =  R4CDist(R4CTerm(t_0_yyy));
    
    EXPECT_EQ(rgraph[4], R4Group({rd_0_xxx, rd_0_xxy, rd_0_xyy, rd_0_yyy,}));
    
    // check sixth recursion term
    
    auto rd_0_xxxx =  R4CDist(R4CTerm(t_0_xxxx));
    
    auto rd_0_xxyy =  R4CDist(R4CTerm(t_0_xxyy));
    
    auto rd_0_yyyy =  R4CDist(R4CTerm(t_0_yyyy));
    
    EXPECT_EQ(rgraph[5], R4Group({rd_0_xxxx, rd_0_xxyy, rd_0_yyyy,}));
    
    // check common integrals set
    
    std::set<T4CIntegral> rints = {t_0_xx, t_0_xy, t_0_yy,
                                   t_0_xxx, t_0_xxy, t_0_xyy, t_0_yyy,
                                   t_0_xxxx, t_0_xxyy, t_0_yyyy,
                                   t_x_xx, t_y_xy, t_y_yy,
                                   t_x_xxx, t_y_xxy, t_y_yyy};
    
    EXPECT_EQ(sints, rints);
    
    // check edges
    
    EXPECT_EQ(rgraph.edge(0), std::set<int>({1, 2}));
    
    EXPECT_EQ(rgraph.edge(1), std::set<int>({3, 4}));
    
    EXPECT_EQ(rgraph.edge(2), std::set<int>({4, 5}));
    
    EXPECT_EQ(rgraph.edge(3), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(4), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(5), std::set<int>({}));
}


//const auto nverts = rgraph.vertices();
//
//for (int i = 0; i < nverts; i++)
//{
//    std::cout << "Vertice: " << i << std::endl;
//
//    for (const auto& tval : rgraph[i].roots())
//    {
//        std::cout << tval.integral().label() << " ";
//    }
//
//    std::cout << std::endl;
//
//    for (const auto& tval : rgraph.edge(i))
//    {
//        std::cout << tval << " ";
//    }
//
//    std::cout << std::endl;
//}
