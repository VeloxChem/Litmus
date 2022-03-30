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

#include "repository.hpp"

TEST_F(EriDriverTest, BraHrr)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_zz = TensorComponent(0, 0, 2);
    
    const auto f_xzz = TensorComponent(1, 0, 2);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    // bra and ket pairs
    
    const auto b_xy_zz = T2CPair({"GA", "GB"}, {d_xy, d_zz});
    
    const auto b_y_xzz = T2CPair({"GA", "GB"}, {p_y, f_xzz});
    
    const auto b_x_yzz = T2CPair({"GA", "GB"}, {p_x, f_yzz});
    
    const auto b_y_zz = T2CPair({"GA", "GB"}, {p_y, d_zz});
    
    const auto b_x_zz = T2CPair({"GA", "GB"}, {p_x, d_zz});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto tint = T4CIntegral(b_xy_zz, k_0_0, operi);
    
    const auto r1aint = T4CIntegral(b_y_xzz, k_0_0, operi);
    
    const auto r2aint = T4CIntegral(b_y_zz, k_0_0, operi);
    
    const auto r1bint = T4CIntegral(b_x_yzz, k_0_0, operi);
    
    const auto r2bint = T4CIntegral(b_x_zz, k_0_0, operi);
    
    // recursion factors
    
    const auto abx = Factor("AB", "rab", TensorComponent(1, 0, 0));
    
    const auto aby = Factor("AB", "rab", TensorComponent(0, 1, 0));
    
    // reference recursion term
    
    const auto t4crec = R4CTerm(tint);
    
    // check recursion along x axis
    
    const auto t1arec = R4CTerm(r1aint);
    
    const auto t2arec = R4CTerm(r2aint, {{abx, 1}, }, Fraction(-1));
    
    EXPECT_EQ(eri_drv.bra_hrr(t4crec, 'x'), R4CDist(t4crec, {t1arec, t2arec}));
    
    // check recursion along y axis
    
    const auto t1brec = R4CTerm(r1bint);

    const auto t2brec = R4CTerm(r2bint, {{aby, 1}, }, Fraction(-1));
    
    EXPECT_EQ(eri_drv.bra_hrr(t4crec, 'y'), R4CDist(t4crec, {t1brec, t2brec}));
    
    // check recursion along z axis
    
    EXPECT_FALSE(eri_drv.bra_hrr(t4crec, 'z'));
}

TEST_F(EriDriverTest, KetHrr)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_zz = TensorComponent(0, 0, 2);
    
    const auto f_xzz = TensorComponent(1, 0, 2);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    // bra and ket pairs
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_xy_zz = T2CPair({"GC", "GD"}, {d_xy, d_zz});
    
    const auto k_y_xzz = T2CPair({"GC", "GD"}, {p_y, f_xzz});
    
    const auto k_x_yzz = T2CPair({"GC", "GD"}, {p_x, f_yzz});
    
    const auto k_y_zz = T2CPair({"GC", "GD"}, {p_y, d_zz});
    
    const auto k_x_zz = T2CPair({"GC", "GD"}, {p_x, d_zz});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto tint = T4CIntegral(b_0_0, k_xy_zz, operi);
    
    const auto r1aint = T4CIntegral(b_0_0, k_y_xzz, operi);
    
    const auto r2aint = T4CIntegral(b_0_0, k_y_zz, operi);
    
    const auto r1bint = T4CIntegral(b_0_0, k_x_yzz, operi);
    
    const auto r2bint = T4CIntegral(b_0_0, k_x_zz, operi);
    
    // recursion factors
    
    const auto cdx = Factor("CD", "rcd", TensorComponent(1, 0, 0));
    
    const auto cdy = Factor("CD", "rcd", TensorComponent(0, 1, 0));
    
    // reference recursion term
    
    const auto t4crec = R4CTerm(tint);
    
    // check recursion along x axis
    
    const auto t1arec = R4CTerm(r1aint);
    
    const auto t2arec = R4CTerm(r2aint, {{cdx, 1}, }, Fraction(-1));
    
    EXPECT_EQ(eri_drv.ket_hrr(t4crec, 'x'), R4CDist(t4crec, {t1arec, t2arec}));
    
    // check recursion along y axis
    
    const auto t1brec = R4CTerm(r1bint);

    const auto t2brec = R4CTerm(r2bint, {{cdy, 1}, }, Fraction(-1));
    
    EXPECT_EQ(eri_drv.ket_hrr(t4crec, 'y'), R4CDist(t4crec, {t1brec, t2brec}));
    
    // check recursion along z axis
    
    EXPECT_FALSE(eri_drv.ket_hrr(t4crec, 'z'));
}

TEST_F(EriDriverTest, BraVrr)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);

    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);
    
    const auto f_xxx = TensorComponent(3, 0, 0);
    
    const auto f_xyy = TensorComponent(1, 2, 0);
    
    // bra and ket pairs
    
    const auto b_0_xyy = T2CPair({"GA", "GB"}, {s_0, f_xyy});
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto b_0_x = T2CPair({"GA", "GB"}, {s_0, p_x});
    
    const auto k_0_xxx = T2CPair({"GC", "GD"}, {s_0, f_xxx});
    
    const auto k_0_xx = T2CPair({"GC", "GD"}, {s_0, d_xx});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto tint = T4CIntegral(b_0_xyy, k_0_xxx, operi);
    
    const auto r1aint = T4CIntegral(b_0_yy, k_0_xxx, operi);
    
    const auto r2aint = T4CIntegral(b_0_yy, k_0_xxx, operi, 1);
    
    const auto r3aint = T4CIntegral(b_0_yy, k_0_xx, operi, 1);
    
    const auto r1bint = T4CIntegral(b_0_xy, k_0_xxx, operi);
    
    const auto r2bint = T4CIntegral(b_0_xy, k_0_xxx, operi, 1);
    
    const auto r3bint = T4CIntegral(b_0_x, k_0_xxx, operi);
    
    const auto r4bint = T4CIntegral(b_0_x, k_0_xxx, operi, 1);
    
    // recursion factors
    
    const auto pbx = Factor("PB", "rpb", TensorComponent(1, 0, 0));
    
    const auto pby = Factor("PB", "rpb", TensorComponent(0, 1, 0));
    
    const auto wpx = Factor("WP", "rwp", TensorComponent(1, 0, 0));
    
    const auto wpy = Factor("WP", "rwp", TensorComponent(0, 1, 0));
    
    const auto fze = Factor("1/(zeta+eta)", "fze");
    
    const auto fz = Factor("1/zeta", "fz");
    
    const auto frz2 = Factor("rho/zeta^2", "frz2");
    
    // reference recursion term
    
    const auto t4crec = R4CTerm(tint);
    
    // check recursion along x axis
    
    const auto t1arec = R4CTerm(r1aint, {{pbx, 1},}, Fraction(1));
    
    const auto t2arec = R4CTerm(r2aint, {{wpx, 1},}, Fraction(1));
    
    const auto t3arec = R4CTerm(r3aint, {{fze, 1},}, Fraction(3, 2));
    
    EXPECT_EQ(eri_drv.bra_vrr(t4crec, 'x'), R4CDist(t4crec, {t1arec, t2arec, t3arec}));
    
    // check recursion along y axis
    
    const auto t1brec = R4CTerm(r1bint, {{pby, 1},}, Fraction(1));
    
    const auto t2brec = R4CTerm(r2bint, {{wpy, 1},}, Fraction(1));
    
    const auto t3brec = R4CTerm(r3bint, {{fz, 1},}, Fraction(1, 2));
    
    const auto t4brec = R4CTerm(r4bint, {{frz2, 1},}, Fraction(-1, 2));
    
    EXPECT_EQ(eri_drv.bra_vrr(t4crec, 'y'), R4CDist(t4crec, {t1brec, t2brec, t3brec, t4brec}));
    
    // check recursion along z axis
    
    EXPECT_FALSE(eri_drv.bra_vrr(t4crec, 'z'));
}

TEST_F(EriDriverTest, KetVrr)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);
    
    const auto f_xxx = TensorComponent(3, 0, 0);
    
    // bra and ket pairs
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_0_xxx = T2CPair({"GC", "GD"}, {s_0, f_xxx});
    
    const auto k_0_xx = T2CPair({"GC", "GD"}, {s_0, d_xx});
    
    const auto k_0_x = T2CPair({"GC", "GD"}, {s_0, p_x});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto tint = T4CIntegral(b_0_0, k_0_xxx, operi);
    
    const auto r1aint = T4CIntegral(b_0_0, k_0_xx, operi);
    
    const auto r2aint = T4CIntegral(b_0_0, k_0_xx, operi, 1);
    
    const auto r3aint = T4CIntegral(b_0_0, k_0_x, operi);
    
    const auto r4aint = T4CIntegral(b_0_0, k_0_x, operi, 1);
    
    // recursion factors
    
    const auto qdx = Factor("QD", "rqd", TensorComponent(1, 0, 0));
    
    const auto wqx = Factor("WQ", "rwq", TensorComponent(1, 0, 0));
    
    const auto fe = Factor("1/eta", "fe");
    
    const auto fre2 = Factor("rho/eta^2", "fre2");
    
    // reference recursion term
    
    const auto t4crec = R4CTerm(tint);
    
    // check recursion along x axis
    
    const auto t1arec = R4CTerm(r1aint, {{qdx, 1},}, Fraction(1));
    
    const auto t2arec = R4CTerm(r2aint, {{wqx, 1},}, Fraction(1));
    
    const auto t3arec = R4CTerm(r3aint, {{fe, 1},}, Fraction(1));
    
    const auto t4arec = R4CTerm(r4aint, {{fre2, 1},}, Fraction(-1));
    
    EXPECT_EQ(eri_drv.ket_vrr(t4crec, 'x'), R4CDist(t4crec, {t1arec, t2arec, t3arec, t4arec}));
    
    // check recursion along y axis
    
    EXPECT_FALSE(eri_drv.ket_vrr(t4crec, 'y'));
    
    // check recursion along z axis
    
    EXPECT_FALSE(eri_drv.ket_vrr(t4crec, 'z'));
}

TEST_F(EriDriverTest, ApplyBraHrr)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_zz = TensorComponent(0, 0, 2);
    
    const auto f_xzz = TensorComponent(1, 0, 2);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    // bra and ket pairs
    
    const auto b_xy_zz = T2CPair({"GA", "GB"}, {d_xy, d_zz});
    
    const auto b_y_xzz = T2CPair({"GA", "GB"}, {p_y, f_xzz});
    
    const auto b_x_yzz = T2CPair({"GA", "GB"}, {p_x, f_yzz});
    
    const auto b_y_zz = T2CPair({"GA", "GB"}, {p_y, d_zz});
    
    const auto b_x_zz = T2CPair({"GA", "GB"}, {p_x, d_zz});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto tint = T4CIntegral(b_xy_zz, k_0_0, operi);
    
    const auto r1aint = T4CIntegral(b_y_xzz, k_0_0, operi);
    
    const auto r2aint = T4CIntegral(b_y_zz, k_0_0, operi);
    
    const auto r1bint = T4CIntegral(b_x_yzz, k_0_0, operi);
    
    const auto r2bint = T4CIntegral(b_x_zz, k_0_0, operi);
    
    // recursion factors
    
    const auto abx = Factor("AB", "rab", TensorComponent(1, 0, 0));
    
    const auto aby = Factor("AB", "rab", TensorComponent(0, 1, 0));
    
    // reference recursion term
    
    const auto t4crec = R4CTerm(tint);
    
    // witout initial set of integrals
    
    const auto t1arec = R4CTerm(r1aint);
    
    const auto t2arec = R4CTerm(r2aint, {{abx, 1}, }, Fraction(-1));
    
    std::set<T4CIntegral> sints;
    
    auto r4cdist = eri_drv.apply_bra_hrr(t4crec, sints);
    
    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1arec, t2arec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1aint, r2aint}));
    
    // with initial set of integrals
    
    const auto t1brec = R4CTerm(r1bint);

    const auto t2brec = R4CTerm(r2bint, {{aby, 1}, }, Fraction(-1));
    
    sints = {r2bint,};
    
    r4cdist = eri_drv.apply_bra_hrr(t4crec, sints);

    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1brec, t2brec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1bint, r2bint}));
}

TEST_F(EriDriverTest, ApplyKetHrr)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_zz = TensorComponent(0, 0, 2);
    
    const auto f_xzz = TensorComponent(1, 0, 2);
    
    const auto f_yzz = TensorComponent(0, 1, 2);
    
    // bra and ket pairs
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_xy_zz = T2CPair({"GC", "GD"}, {d_xy, d_zz});
    
    const auto k_y_xzz = T2CPair({"GC", "GD"}, {p_y, f_xzz});
    
    const auto k_x_yzz = T2CPair({"GC", "GD"}, {p_x, f_yzz});
    
    const auto k_y_zz = T2CPair({"GC", "GD"}, {p_y, d_zz});
    
    const auto k_x_zz = T2CPair({"GC", "GD"}, {p_x, d_zz});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto tint = T4CIntegral(b_0_0, k_xy_zz, operi);
    
    const auto r1aint = T4CIntegral(b_0_0, k_y_xzz, operi);
    
    const auto r2aint = T4CIntegral(b_0_0, k_y_zz, operi);
    
    const auto r1bint = T4CIntegral(b_0_0, k_x_yzz, operi);
    
    const auto r2bint = T4CIntegral(b_0_0, k_x_zz, operi);
    
    // recursion factors
    
    const auto cdx = Factor("CD", "rcd", TensorComponent(1, 0, 0));
    
    const auto cdy = Factor("CD", "rcd", TensorComponent(0, 1, 0));
    
    // reference recursion term
    
    const auto t4crec = R4CTerm(tint);
    
    // witout initial set of integrals
    
    const auto t1arec = R4CTerm(r1aint);
    
    const auto t2arec = R4CTerm(r2aint, {{cdx, 1}, }, Fraction(-1));
    
    std::set<T4CIntegral> sints;
    
    auto r4cdist = eri_drv.apply_ket_hrr(t4crec, sints);
    
    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1arec, t2arec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1aint, r2aint}));
    
    // with initial set of integrals
    
    const auto t1brec = R4CTerm(r1bint);
    
    const auto t2brec = R4CTerm(r2bint, {{cdy, 1}, }, Fraction(-1));
    
    sints = {r2bint,};
    
    r4cdist = eri_drv.apply_ket_hrr(t4crec, sints);

    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1brec, t2brec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1bint, r2bint}));
}

TEST_F(EriDriverTest, ApplyBraVrr)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);

    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);
    
    const auto f_xxx = TensorComponent(3, 0, 0);
    
    const auto f_xyy = TensorComponent(1, 2, 0);
    
    // bra and ket pairs
    
    const auto b_0_xyy = T2CPair({"GA", "GB"}, {s_0, f_xyy});
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto b_0_x = T2CPair({"GA", "GB"}, {s_0, p_x});
    
    const auto k_0_xxx = T2CPair({"GC", "GD"}, {s_0, f_xxx});
    
    const auto k_0_xx = T2CPair({"GC", "GD"}, {s_0, d_xx});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto tint = T4CIntegral(b_0_xyy, k_0_xxx, operi);
    
    const auto r1aint = T4CIntegral(b_0_yy, k_0_xxx, operi);
    
    const auto r2aint = T4CIntegral(b_0_yy, k_0_xxx, operi, 1);
    
    const auto r3aint = T4CIntegral(b_0_yy, k_0_xx, operi, 1);
    
    const auto r1bint = T4CIntegral(b_0_xy, k_0_xxx, operi);
    
    const auto r2bint = T4CIntegral(b_0_xy, k_0_xxx, operi, 1);
    
    const auto r3bint = T4CIntegral(b_0_x, k_0_xxx, operi);
    
    const auto r4bint = T4CIntegral(b_0_x, k_0_xxx, operi, 1);
    
    // recursion factors
    
    const auto pbx = Factor("PB", "rpb", TensorComponent(1, 0, 0));
    
    const auto pby = Factor("PB", "rpb", TensorComponent(0, 1, 0));
    
    const auto wpx = Factor("WP", "rwp", TensorComponent(1, 0, 0));
    
    const auto wpy = Factor("WP", "rwp", TensorComponent(0, 1, 0));
    
    const auto fze = Factor("1/(zeta+eta)", "fze");
    
    const auto fz = Factor("1/zeta", "fz");
    
    const auto frz2 = Factor("rho/zeta^2", "frz2");
    
    // reference recursion term
    
    const auto t4crec = R4CTerm(tint);
    
    // witout initial set of integrals
    
    const auto t1arec = R4CTerm(r1aint, {{pbx, 1},}, Fraction(1));
    
    const auto t2arec = R4CTerm(r2aint, {{wpx, 1},}, Fraction(1));
    
    const auto t3arec = R4CTerm(r3aint, {{fze, 1},}, Fraction(3, 2));
    
    std::set<T4CIntegral> sints;
    
    auto r4cdist = eri_drv.apply_bra_vrr(t4crec, sints);
    
    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1arec, t2arec, t3arec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1aint, r2aint, r3aint}));
    
    // with initial set of integrals
    
    const auto t1brec = R4CTerm(r1bint, {{pby, 1},}, Fraction(1));
    
    const auto t2brec = R4CTerm(r2bint, {{wpy, 1},}, Fraction(1));
    
    const auto t3brec = R4CTerm(r3bint, {{fz, 1},}, Fraction(1, 2));
    
    const auto t4brec = R4CTerm(r4bint, {{frz2, 1},}, Fraction(-1, 2));
    
    sints = {r1bint, r2bint, r4bint};
    
    r4cdist = eri_drv.apply_bra_vrr(t4crec, sints);

    EXPECT_EQ(r4cdist, R4CDist(t4crec, {t1brec, t2brec, t3brec, t4brec}));
    
    EXPECT_EQ(sints, std::set<T4CIntegral>({r1bint, r2bint, r3bint, r4bint}));
}

TEST_F(EriDriverTest, AppyKetVrr)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);
    
    const auto f_xxx = TensorComponent(3, 0, 0);
    
    // bra and ket pairs
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_0_xxx = T2CPair({"GC", "GD"}, {s_0, f_xxx});
    
    const auto k_0_xx = T2CPair({"GC", "GD"}, {s_0, d_xx});
    
    const auto k_0_x = T2CPair({"GC", "GD"}, {s_0, p_x});
    
    // operator

    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto tint = T4CIntegral(b_0_0, k_0_xxx, operi);
    
    const auto r1aint = T4CIntegral(b_0_0, k_0_xx, operi);
    
    const auto r2aint = T4CIntegral(b_0_0, k_0_xx, operi, 1);
    
    const auto r3aint = T4CIntegral(b_0_0, k_0_x, operi);
    
    const auto r4aint = T4CIntegral(b_0_0, k_0_x, operi, 1);
    
    // recursion factors
    
    const auto qdx = Factor("QD", "rqd", TensorComponent(1, 0, 0));
    
    const auto wqx = Factor("WQ", "rwq", TensorComponent(1, 0, 0));
    
    const auto fe = Factor("1/eta", "fe");
    
    const auto fre2 = Factor("rho/eta^2", "fre2");
    
    // reference recursion term
    
    const auto t4crec = R4CTerm(tint);
    
    // witout initial set of integrals
    
    const auto t1arec = R4CTerm(r1aint, {{qdx, 1},}, Fraction(1));
    
    const auto t2arec = R4CTerm(r2aint, {{wqx, 1},}, Fraction(1));
    
    const auto t3arec = R4CTerm(r3aint, {{fe, 1},}, Fraction(1));
    
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
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    // bra and ket pairs
    
    const auto b_x_x = T2CPair({"GA", "GB"}, {p_x, p_x});
    
    const auto b_y_x = T2CPair({"GA", "GB"}, {p_y, p_x});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto taint = T4CIntegral(b_x_x, k_0_0, operi);
    
    const auto tbint = T4CIntegral(b_y_x, k_0_0, operi);
    
    // generated recursion group
    
    const auto t4arec = R4CTerm(taint);
    
    const auto t4brec = R4CTerm(tbint);
    
    std::set<T4CIntegral> sints;
    
    const auto t4g = eri_drv.apply_bra_hrr({t4arec, t4brec}, sints);
    
    // reference recursion group
    
    std::set<T4CIntegral> rints;
    
    const auto r4adist = eri_drv.apply_bra_hrr(t4arec, rints);
    
    const auto r4bdist = eri_drv.apply_bra_hrr(t4brec, rints);
    
    EXPECT_EQ(t4g, R4Group({r4adist, r4bdist}));
    
    EXPECT_EQ(sints, rints);
}

TEST_F(EriDriverTest, ApplyKetHrrForGroup)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    // bra and ket pairs
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_x_x = T2CPair({"GC", "GD"}, {p_x, p_x});
    
    const auto k_y_x = T2CPair({"GC", "GD"}, {p_y, p_x});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
   
    // integral components
    
    const auto taint = T4CIntegral(b_0_0, k_x_x, operi);
    
    const auto tbint = T4CIntegral(b_0_0, k_y_x, operi);
    
    // generated recursion group
    
    const auto t4arec = R4CTerm(taint);
    
    const auto t4brec = R4CTerm(tbint);
    
    std::set<T4CIntegral> sints;
    
    const auto t4g = eri_drv.apply_ket_hrr({t4arec, t4brec}, sints);
    
    // reference recursion group
    
    std::set<T4CIntegral> rints;
    
    const auto r4adist = eri_drv.apply_ket_hrr(t4arec, rints);
    
    const auto r4bdist = eri_drv.apply_ket_hrr(t4brec, rints);
    
    EXPECT_EQ(t4g, R4Group({r4adist, r4bdist}));
    
    EXPECT_EQ(sints, rints);
}

TEST_F(EriDriverTest, ApplyBraVrrForGroup)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);
    
    const auto f_xyy = TensorComponent(1, 3, 0);
    
    const auto f_xyz = TensorComponent(1, 1, 1);
    
    // bra and ket pairs
    
    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto k_0_xyy = T2CPair({"GC", "GD"}, {s_0, f_xyy});
                                       
    const auto k_0_xyz = T2CPair({"GC", "GD"}, {s_0, f_xyz});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
                                       
    const auto taint = T4CIntegral(b_0_xy, k_0_xyy, operi);
    
    const auto tbint = T4CIntegral(b_0_yy, k_0_xyz, operi);
    
    // generated recursion group
    
    const auto t4arec = R4CTerm(taint);
    
    const auto t4brec = R4CTerm(tbint);
    
    std::set<T4CIntegral> sints;
    
    const auto t4g = eri_drv.apply_bra_vrr({t4arec, t4brec}, sints);
    
    // reference recursion group
    
    std::set<T4CIntegral> rints;
    
    const auto r4adist = eri_drv.apply_bra_vrr(t4arec, rints);
    
    const auto r4bdist = eri_drv.apply_bra_vrr(t4brec, rints);
    
    EXPECT_EQ(t4g, R4Group({r4adist, r4bdist}));
    
    EXPECT_EQ(sints, rints);
}

TEST_F(EriDriverTest, ApplyKetVrrForGroup)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto f_xyy = TensorComponent(1, 3, 0);
    
    const auto f_xyz = TensorComponent(1, 1, 1);
    
    // bra and ket pairs
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_0_xyy = T2CPair({"GC", "GD"}, {s_0, f_xyy});
                                       
    const auto k_0_xyz = T2CPair({"GC", "GD"}, {s_0, f_xyz});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
                                       
    const auto taint = T4CIntegral(b_0_0, k_0_xyy, operi);
    
    const auto tbint = T4CIntegral(b_0_0, k_0_xyz, operi);
    
    // generated recursion group
    
    const auto t4arec = R4CTerm(taint);
    
    const auto t4brec = R4CTerm(tbint);
    
    std::set<T4CIntegral> sints;
    
    const auto t4g = eri_drv.apply_ket_vrr({t4arec, t4brec}, sints);
    
    // reference recursion group
    
    std::set<T4CIntegral> rints;
    
    const auto r4adist = eri_drv.apply_ket_vrr(t4arec, rints);
    
    const auto r4bdist = eri_drv.apply_ket_vrr(t4brec, rints);
    
    EXPECT_EQ(t4g, R4Group({r4adist, r4bdist}));
    
    EXPECT_EQ(sints, rints);
}

TEST_F(EriDriverTest, ApplyBraHrrWithGraphForPP)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);
    
    // bra and ket pairs
    
    const auto b_x_x = T2CPair({"GA", "GB"}, {p_x, p_x});
    
    const auto b_x_y = T2CPair({"GA", "GB"}, {p_x, p_y});
    
    const auto b_y_y = T2CPair({"GA", "GB"}, {p_y, p_y});
    
    const auto b_0_xx = T2CPair({"GA", "GB"}, {s_0, d_xx});
    
    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto b_0_x = T2CPair({"GA", "GB"}, {s_0, p_x});
    
    const auto b_0_y = T2CPair({"GA", "GB"}, {s_0, p_y});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto t_x_x = T4CIntegral(b_x_x, k_0_0, operi);
    
    const auto t_x_y = T4CIntegral(b_x_y, k_0_0, operi);
    
    const auto t_y_y = T4CIntegral(b_y_y, k_0_0, operi);
    
    const auto t_0_xx = T4CIntegral(b_0_xx, k_0_0, operi);
    
    const auto t_0_xy = T4CIntegral(b_0_xy, k_0_0, operi);
    
    const auto t_0_yy = T4CIntegral(b_0_yy, k_0_0, operi);
    
    const auto t_0_x = T4CIntegral(b_0_x, k_0_0, operi);
    
    const auto t_0_y = T4CIntegral(b_0_y, k_0_0, operi);
    
    // generate graph
    
    const auto rd_x_x = R4CDist(R4CTerm(t_x_x));
    
    const auto rd_x_y = R4CDist(R4CTerm(t_x_y));
    
    const auto rd_y_y = R4CDist(R4CTerm(t_y_y));
    
    R4Graph rgraph(R4Group({rd_x_x, rd_x_y, rd_y_y}));
    
    std::set<T4CIntegral> sints;
    
    eri_drv.apply_bra_hrr(rgraph, sints);
    
    // set up reference data
    
    std::set<T4CIntegral> rints;
    
    const auto rr_x_x = eri_drv.apply_bra_hrr(R4CTerm(t_x_x), rints);
    
    const auto rr_x_y = eri_drv.apply_bra_hrr(R4CTerm(t_x_y), rints);
    
    const auto rr_y_y = eri_drv.apply_bra_hrr(R4CTerm(t_y_y), rints);
    
    const auto rd_0_xx = R4CDist(R4CTerm(t_0_xx));
    
    const auto rd_0_xy = R4CDist(R4CTerm(t_0_xy));
    
    const auto rd_0_yy = R4CDist(R4CTerm(t_0_yy));
    
    const auto rd_0_x = R4CDist(R4CTerm(t_0_x));
    
    const auto rd_0_y = R4CDist(R4CTerm(t_0_y));
    
    // compare vertices and edges of graph
    
    EXPECT_EQ(rgraph.vertices(), 3);
    
    EXPECT_EQ(rgraph[0], R4Group({rr_x_x, rr_x_y, rr_y_y}));
    
    EXPECT_EQ(rgraph[1], R4Group({rd_0_yy, rd_0_xy, rd_0_xx}));
    
    EXPECT_EQ(rgraph[2], R4Group({rd_0_y, rd_0_x}));
    
    EXPECT_EQ(sints, rints);
    
    EXPECT_EQ(rgraph.edge(0), std::set<int>({1, 2}));
    
    EXPECT_EQ(rgraph.edge(1), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(2), std::set<int>({}));
}

TEST_F(EriDriverTest, ApplyBraHrrWithGraphForDD)
{
    EriDriver eri_drv;
    
    // tensor components
    
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
    
    // bra and ket pairs
    
    const auto b_xx_xx = T2CPair({"GA", "GB"}, {d_xx, d_xx});
    
    const auto b_xy_xy = T2CPair({"GA", "GB"}, {d_xy, d_xy});
    
    const auto b_yy_yy = T2CPair({"GA", "GB"}, {d_yy, d_yy});
    
    const auto b_x_xxx = T2CPair({"GA", "GB"}, {p_x, f_xxx});
    
    const auto b_y_xxy = T2CPair({"GA", "GB"}, {p_y, f_xxy});
    
    const auto b_y_yyy = T2CPair({"GA", "GB"}, {p_y, f_yyy});
    
    const auto b_x_xx = T2CPair({"GA", "GB"}, {p_x, d_xx});
    
    const auto b_y_xy = T2CPair({"GA", "GB"}, {p_y, d_xy});
    
    const auto b_y_yy = T2CPair({"GA", "GB"}, {p_y, d_yy});
    
    const auto b_0_xxxx = T2CPair({"GA", "GB"}, {s_0, g_xxxx});
    
    const auto b_0_xxyy = T2CPair({"GA", "GB"}, {s_0, g_xxyy});
    
    const auto b_0_yyyy = T2CPair({"GA", "GB"}, {s_0, g_yyyy});
    
    const auto b_0_xxx = T2CPair({"GA", "GB"}, {s_0, f_xxx});
    
    const auto b_0_xxy = T2CPair({"GA", "GB"}, {s_0, f_xxy});
    
    const auto b_0_xyy = T2CPair({"GA", "GB"}, {s_0, f_xyy});
    
    const auto b_0_yyy = T2CPair({"GA", "GB"}, {s_0, f_yyy});
    
    const auto b_0_xx = T2CPair({"GA", "GB"}, {s_0, d_xx});
    
    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto t_xx_xx = T4CIntegral(b_xx_xx, k_0_0, operi);
    
    const auto t_xy_xy = T4CIntegral(b_xy_xy, k_0_0, operi);
    
    const auto t_yy_yy = T4CIntegral(b_yy_yy, k_0_0, operi);
    
    const auto t_x_xxx = T4CIntegral(b_x_xxx, k_0_0, operi);
    
    const auto t_y_xxy = T4CIntegral(b_y_xxy, k_0_0, operi);
    
    const auto t_y_yyy = T4CIntegral(b_y_yyy, k_0_0, operi);
    
    const auto t_x_xx = T4CIntegral(b_x_xx, k_0_0, operi);
    
    const auto t_y_xy = T4CIntegral(b_y_xy, k_0_0, operi);
    
    const auto t_y_yy = T4CIntegral(b_y_yy, k_0_0, operi);
    
    const auto t_0_xxxx = T4CIntegral(b_0_xxxx, k_0_0, operi);
    
    const auto t_0_xxyy = T4CIntegral(b_0_xxyy, k_0_0, operi);
    
    const auto t_0_yyyy = T4CIntegral(b_0_yyyy, k_0_0, operi);
    
    const auto t_0_xxx = T4CIntegral(b_0_xxx, k_0_0, operi);
    
    const auto t_0_xxy = T4CIntegral(b_0_xxy, k_0_0, operi);
    
    const auto t_0_xyy = T4CIntegral(b_0_xyy, k_0_0, operi);
    
    const auto t_0_yyy = T4CIntegral(b_0_yyy, k_0_0, operi);
    
    const auto t_0_xx = T4CIntegral(b_0_xx, k_0_0, operi);
    
    const auto t_0_xy = T4CIntegral(b_0_xy, k_0_0, operi);
    
    const auto t_0_yy = T4CIntegral(b_0_yy, k_0_0, operi);
    
    // generate graph
    
    const auto rd_xx_xx = R4CDist(R4CTerm(t_xx_xx));
    
    const auto rd_xy_xy = R4CDist(R4CTerm(t_xy_xy));
    
    const auto rd_yy_yy = R4CDist(R4CTerm(t_yy_yy));
    
    R4Graph rgraph(R4Group({rd_xx_xx, rd_xy_xy, rd_yy_yy}));
    
    std::set<T4CIntegral> sints;
    
    eri_drv.apply_bra_hrr(rgraph, sints);
    
    // set up reference data
    
    std::set<T4CIntegral> rints;
    
    const auto rr_yy_yy = eri_drv.apply_bra_hrr(R4CTerm(t_yy_yy), rints);
    
    const auto rr_xy_xy = eri_drv.apply_bra_hrr(R4CTerm(t_xy_xy), rints);
    
    const auto rr_xx_xx = eri_drv.apply_bra_hrr(R4CTerm(t_xx_xx), rints);
    
    const auto rr_y_yyy = eri_drv.apply_bra_hrr(R4CTerm(t_y_yyy), rints);
    
    const auto rr_y_xxy = eri_drv.apply_bra_hrr(R4CTerm(t_y_xxy), rints);
    
    const auto rr_x_xxx = eri_drv.apply_bra_hrr(R4CTerm(t_x_xxx), rints);

    const auto rr_y_yy = eri_drv.apply_bra_hrr(R4CTerm(t_y_yy), rints);
    
    const auto rr_y_xy = eri_drv.apply_bra_hrr(R4CTerm(t_y_xy), rints);
    
    const auto rr_x_xx = eri_drv.apply_bra_hrr(R4CTerm(t_x_xx), rints);
    
    const auto rd_0_xxxx = R4CDist(R4CTerm(t_0_xxxx));
    
    const auto rd_0_xxyy = R4CDist(R4CTerm(t_0_xxyy));
    
    const auto rd_0_yyyy = R4CDist(R4CTerm(t_0_yyyy));
    
    const auto rd_0_xxx = R4CDist(R4CTerm(t_0_xxx));
    
    const auto rd_0_xxy = R4CDist(R4CTerm(t_0_xxy));
    
    const auto rd_0_xyy = R4CDist(R4CTerm(t_0_xyy));
    
    const auto rd_0_yyy = R4CDist(R4CTerm(t_0_yyy));
    
    const auto rd_0_xx =  R4CDist(R4CTerm(t_0_xx));
    
    const auto rd_0_xy =  R4CDist(R4CTerm(t_0_xy));
    
    const auto rd_0_yy =  R4CDist(R4CTerm(t_0_yy));
    
    // compare vertices and edges of graph
    
    EXPECT_EQ(rgraph.vertices(), 6);
    
    EXPECT_EQ(rgraph[0], R4Group({rr_xx_xx, rr_xy_xy, rr_yy_yy,}));
    
    EXPECT_EQ(rgraph[1], R4Group({rr_x_xxx, rr_y_xxy, rr_y_yyy,}));
    
    EXPECT_EQ(rgraph[2], R4Group({rr_x_xx, rr_y_xy, rr_y_yy,}));
    
    EXPECT_EQ(rgraph[3], R4Group({rd_0_xxxx, rd_0_xxyy, rd_0_yyyy,}));
    
    EXPECT_EQ(rgraph[4], R4Group({rd_0_xxx, rd_0_xxy, rd_0_xyy, rd_0_yyy,}));
    
    EXPECT_EQ(rgraph[5], R4Group({rd_0_xx, rd_0_xy, rd_0_yy,}));
    
    EXPECT_EQ(sints, rints);
    
    EXPECT_EQ(rgraph.edge(0), std::set<int>({1, 2}));
    
    EXPECT_EQ(rgraph.edge(1), std::set<int>({3, 4}));
    
    EXPECT_EQ(rgraph.edge(2), std::set<int>({4, 5}));
    
    EXPECT_EQ(rgraph.edge(3), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(4), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(5), std::set<int>({}));
}

TEST_F(EriDriverTest, ApplyKetHrrWithGraphForDD)
{
    EriDriver eri_drv;
    
    // tensor components
    
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
    
    // bra and ket pairs
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_xx_xx = T2CPair({"GC", "GD"}, {d_xx, d_xx});
    
    const auto k_xy_xy = T2CPair({"GC", "GD"}, {d_xy, d_xy});
    
    const auto k_yy_yy = T2CPair({"GC", "GD"}, {d_yy, d_yy});
    
    const auto k_x_xxx = T2CPair({"GC", "GD"}, {p_x, f_xxx});
    
    const auto k_y_xxy = T2CPair({"GC", "GD"}, {p_y, f_xxy});
    
    const auto k_y_yyy = T2CPair({"GC", "GD"}, {p_y, f_yyy});
    
    const auto k_x_xx = T2CPair({"GC", "GD"}, {p_x, d_xx});
    
    const auto k_y_xy = T2CPair({"GC", "GD"}, {p_y, d_xy});
    
    const auto k_y_yy = T2CPair({"GC", "GD"}, {p_y, d_yy});
    
    const auto k_0_xxxx = T2CPair({"GC", "GD"}, {s_0, g_xxxx});
    
    const auto k_0_xxyy = T2CPair({"GC", "GD"}, {s_0, g_xxyy});
    
    const auto k_0_yyyy = T2CPair({"GC", "GD"}, {s_0, g_yyyy});
    
    const auto k_0_xxx = T2CPair({"GC", "GD"}, {s_0, f_xxx});
    
    const auto k_0_xxy = T2CPair({"GC", "GD"}, {s_0, f_xxy});
    
    const auto k_0_xyy = T2CPair({"GC", "GD"}, {s_0, f_xyy});
    
    const auto k_0_yyy = T2CPair({"GC", "GD"}, {s_0, f_yyy});
    
    const auto k_0_xx = T2CPair({"GC", "GD"}, {s_0, d_xx});
    
    const auto k_0_xy = T2CPair({"GC", "GD"}, {s_0, d_xy});
    
    const auto k_0_yy = T2CPair({"GC", "GD"}, {s_0, d_yy});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto t_xx_xx = T4CIntegral(b_0_0, k_xx_xx, operi);
    
    const auto t_xy_xy = T4CIntegral(b_0_0, k_xy_xy, operi);
    
    const auto t_yy_yy = T4CIntegral(b_0_0, k_yy_yy, operi);
    
    const auto t_x_xxx = T4CIntegral(b_0_0, k_x_xxx, operi);
    
    const auto t_y_xxy = T4CIntegral(b_0_0, k_y_xxy, operi);
    
    const auto t_y_yyy = T4CIntegral(b_0_0, k_y_yyy, operi);
    
    const auto t_x_xx = T4CIntegral(b_0_0, k_x_xx, operi);
    
    const auto t_y_xy = T4CIntegral(b_0_0, k_y_xy, operi);
    
    const auto t_y_yy = T4CIntegral(b_0_0, k_y_yy, operi);
    
    const auto t_0_xxxx = T4CIntegral(b_0_0, k_0_xxxx, operi);
    
    const auto t_0_xxyy = T4CIntegral(b_0_0, k_0_xxyy, operi);
    
    const auto t_0_yyyy = T4CIntegral(b_0_0, k_0_yyyy, operi);
    
    const auto t_0_xxx = T4CIntegral(b_0_0, k_0_xxx, operi);
    
    const auto t_0_xxy = T4CIntegral(b_0_0, k_0_xxy, operi);
    
    const auto t_0_xyy = T4CIntegral(b_0_0, k_0_xyy, operi);
    
    const auto t_0_yyy = T4CIntegral(b_0_0, k_0_yyy, operi);
    
    const auto t_0_xx = T4CIntegral(b_0_0, k_0_xx, operi);
    
    const auto t_0_xy = T4CIntegral(b_0_0, k_0_xy,  operi);
    
    const auto t_0_yy = T4CIntegral(b_0_0, k_0_yy, operi);
    
    // generate graph
    
    const auto rd_xx_xx =  R4CDist(R4CTerm(t_xx_xx));
    
    const auto rd_xy_xy =  R4CDist(R4CTerm(t_xy_xy));
    
    const auto rd_yy_yy =  R4CDist(R4CTerm(t_yy_yy));
    
    R4Graph rgraph(R4Group({rd_xx_xx, rd_xy_xy, rd_yy_yy}));
    
    std::set<T4CIntegral> sints;
    
    eri_drv.apply_ket_hrr(rgraph, sints);
    
    // set up reference data
    
    std::set<T4CIntegral> rints;
    
    const auto rr_yy_yy = eri_drv.apply_ket_hrr(R4CTerm(t_yy_yy), rints);
    
    const auto rr_xy_xy = eri_drv.apply_ket_hrr(R4CTerm(t_xy_xy), rints);
    
    const auto rr_xx_xx = eri_drv.apply_ket_hrr(R4CTerm(t_xx_xx), rints);
    
    const auto rr_y_yyy = eri_drv.apply_ket_hrr(R4CTerm(t_y_yyy), rints);
    
    const auto rr_y_xxy = eri_drv.apply_ket_hrr(R4CTerm(t_y_xxy), rints);
    
    const auto rr_x_xxx = eri_drv.apply_ket_hrr(R4CTerm(t_x_xxx), rints);

    const auto rr_y_yy = eri_drv.apply_ket_hrr(R4CTerm(t_y_yy), rints);
    
    const auto rr_y_xy = eri_drv.apply_ket_hrr(R4CTerm(t_y_xy), rints);
    
    const auto rr_x_xx = eri_drv.apply_ket_hrr(R4CTerm(t_x_xx), rints);
    
    const auto rd_0_xxxx = R4CDist(R4CTerm(t_0_xxxx));
    
    const auto rd_0_xxyy = R4CDist(R4CTerm(t_0_xxyy));
    
    const auto rd_0_yyyy = R4CDist(R4CTerm(t_0_yyyy));
    
    const auto rd_0_xxx = R4CDist(R4CTerm(t_0_xxx));
    
    const auto rd_0_xxy = R4CDist(R4CTerm(t_0_xxy));
    
    const auto rd_0_xyy = R4CDist(R4CTerm(t_0_xyy));
    
    const auto rd_0_yyy = R4CDist(R4CTerm(t_0_yyy));
    
    const auto rd_0_xx =  R4CDist(R4CTerm(t_0_xx));
    
    const auto rd_0_xy =  R4CDist(R4CTerm(t_0_xy));
    
    const auto rd_0_yy =  R4CDist(R4CTerm(t_0_yy));
    
    // compare vertices and edges of graph
    
    EXPECT_EQ(rgraph.vertices(), 6);
    
    EXPECT_EQ(rgraph[0], R4Group({rr_xx_xx, rr_xy_xy, rr_yy_yy,}));
    
    EXPECT_EQ(rgraph[1], R4Group({rr_x_xxx, rr_y_xxy, rr_y_yyy,}));
    
    EXPECT_EQ(rgraph[2], R4Group({rr_x_xx, rr_y_xy, rr_y_yy,}));
    
    EXPECT_EQ(rgraph[3], R4Group({rd_0_xxxx, rd_0_xxyy, rd_0_yyyy,}));
    
    EXPECT_EQ(rgraph[4], R4Group({rd_0_xxx, rd_0_xxy, rd_0_xyy, rd_0_yyy,}));
    
    EXPECT_EQ(rgraph[5], R4Group({rd_0_xx, rd_0_xy, rd_0_yy,}));
    
    EXPECT_EQ(sints, rints);
    
    EXPECT_EQ(rgraph.edge(0), std::set<int>({1, 2}));
    
    EXPECT_EQ(rgraph.edge(1), std::set<int>({3, 4}));
    
    EXPECT_EQ(rgraph.edge(2), std::set<int>({4, 5}));
    
    EXPECT_EQ(rgraph.edge(3), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(4), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(5), std::set<int>({}));
}

TEST_F(EriDriverTest, ApplyBraVrrWithGraphForDD)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);
    
    // bra and ket pairs
    
    const auto b_0_xx = T2CPair({"GA", "GB"}, {s_0, d_xx});

    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto b_0_x = T2CPair({"GA", "GB"}, {s_0, p_x});

    const auto b_0_y = T2CPair({"GA", "GB"}, {s_0, p_y});
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_0_xx = T2CPair({"GC", "GD"}, {s_0, d_xx});
    
    const auto k_0_xy = T2CPair({"GC", "GD"}, {s_0, d_xy});
    
    const auto k_0_yy = T2CPair({"GC", "GD"}, {s_0, d_yy});
    
    const auto k_0_x = T2CPair({"GC", "GD"}, {s_0, p_x});
    
    const auto k_0_y = T2CPair({"GC", "GD"}, {s_0, p_y});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto t_xx_xx = T4CIntegral(b_0_xx, k_0_xx, operi);
    
    const auto t_xy_xy = T4CIntegral(b_0_xy, k_0_xy, operi);
    
    const auto t_yy_yy = T4CIntegral(b_0_yy, k_0_yy, operi);
    
    const auto t_x_xx_0 = T4CIntegral(b_0_x, k_0_xx, operi);
    
    const auto t_y_xy_0 = T4CIntegral(b_0_y, k_0_xy, operi);
    
    const auto t_y_yy_0 = T4CIntegral(b_0_y, k_0_yy, operi);
    
    const auto t_x_xx_1 = T4CIntegral(b_0_x, k_0_xx, operi, 1);
    
    const auto t_y_xy_1 = T4CIntegral(b_0_y, k_0_xy, operi, 1);
    
    const auto t_y_yy_1 = T4CIntegral(b_0_y, k_0_yy, operi, 1);
    
    const auto t_0_xx_0 = T4CIntegral(b_0_0, k_0_xx, operi);
    
    const auto t_0_xy_0 = T4CIntegral(b_0_0, k_0_xy, operi);
    
    const auto t_0_yy_0 = T4CIntegral(b_0_0, k_0_yy, operi);
    
    const auto t_0_xx_1 = T4CIntegral(b_0_0, k_0_xx, operi, 1);
    
    const auto t_0_xy_1 = T4CIntegral(b_0_0, k_0_xy, operi, 1);
    
    const auto t_0_yy_1 = T4CIntegral(b_0_0, k_0_yy, operi, 1);
    
    const auto t_0_xx_2 = T4CIntegral(b_0_0, k_0_xx, operi, 2);
    
    const auto t_0_xy_2 = T4CIntegral(b_0_0, k_0_xy, operi, 2);
    
    const auto t_0_yy_2 = T4CIntegral(b_0_0, k_0_yy, operi, 2);
    
    const auto t_x_x_1 = T4CIntegral(b_0_x, k_0_x, operi, 1);
    
    const auto t_y_y_1 = T4CIntegral(b_0_y, k_0_y, operi, 1);
    
    const auto t_0_x_2 = T4CIntegral(b_0_0, k_0_x, operi, 2);
    
    const auto t_0_y_2 = T4CIntegral(b_0_0, k_0_y, operi, 2);
    
    const auto t_0_x_1 = T4CIntegral(b_0_0, k_0_x, operi, 1);
    
    const auto t_0_y_1 = T4CIntegral(b_0_0, k_0_y, operi, 1);
    
    const auto t_0_0_2 = T4CIntegral(b_0_0, k_0_0, operi, 2);
    
    // generate graph
    
    const auto rd_xx_xx =  R4CDist(R4CTerm(t_xx_xx));
    
    const auto rd_xy_xy =  R4CDist(R4CTerm(t_xy_xy));
    
    const auto rd_yy_yy =  R4CDist(R4CTerm(t_yy_yy));
    
    R4Graph rgraph(R4Group({rd_xx_xx, rd_xy_xy, rd_yy_yy}));
    
    std::set<T4CIntegral> sints;
    
    eri_drv.apply_bra_vrr(rgraph, sints);
    
    // set up reference data
    
    std::set<T4CIntegral> rints;
    
    const auto rr_yy_yy = eri_drv.apply_bra_vrr(R4CTerm(t_yy_yy), rints);
    
    const auto rr_xy_xy = eri_drv.apply_bra_vrr(R4CTerm(t_xy_xy), rints);
    
    const auto rr_xx_xx = eri_drv.apply_bra_vrr(R4CTerm(t_xx_xx), rints);
    
    const auto rr_y_yy_1 = eri_drv.apply_bra_vrr(R4CTerm(t_y_yy_1), rints);
    
    const auto rr_y_xy_1 = eri_drv.apply_bra_vrr(R4CTerm(t_y_xy_1), rints);
    
    const auto rr_x_xx_1 = eri_drv.apply_bra_vrr(R4CTerm(t_x_xx_1), rints);
    
    const auto rr_y_yy_0 = eri_drv.apply_bra_vrr(R4CTerm(t_y_yy_0), rints);
    
    const auto rr_y_xy_0 = eri_drv.apply_bra_vrr(R4CTerm(t_y_xy_0), rints);
    
    const auto rr_x_xx_0 = eri_drv.apply_bra_vrr(R4CTerm(t_x_xx_0), rints);
    
    const auto rr_y_y_1 = eri_drv.apply_bra_vrr(R4CTerm(t_y_y_1), rints);
    
    const auto rr_x_x_1 = eri_drv.apply_bra_vrr(R4CTerm(t_x_x_1), rints);
    
    const auto rd_0_yy_2 = R4CDist(R4CTerm(t_0_yy_2));
    
    const auto rd_0_xy_2 = R4CDist(R4CTerm(t_0_xy_2));
    
    const auto rd_0_xx_2 = R4CDist(R4CTerm(t_0_xx_2));
    
    const auto rd_0_yy_1 = R4CDist(R4CTerm(t_0_yy_1));
    
    const auto rd_0_xy_1 = R4CDist(R4CTerm(t_0_xy_1));
    
    const auto rd_0_xx_1 = R4CDist(R4CTerm(t_0_xx_1));
    
    const auto rd_0_yy_0 = R4CDist(R4CTerm(t_0_yy_0));
    
    const auto rd_0_xy_0 = R4CDist(R4CTerm(t_0_xy_0));
    
    const auto rd_0_xx_0 = R4CDist(R4CTerm(t_0_xx_0));
    
    const auto rd_0_y_2 = R4CDist(R4CTerm(t_0_y_2));
    
    const auto rd_0_x_2 = R4CDist(R4CTerm(t_0_x_2));
    
    const auto rd_0_y_1 = R4CDist(R4CTerm(t_0_y_1));
    
    const auto rd_0_x_1 = R4CDist(R4CTerm(t_0_x_1));
    
    const auto rd_0_0_2 = R4CDist(R4CTerm(t_0_0_2));
    
    // compare vertices and edges of graph
    
    EXPECT_EQ(rgraph.vertices(), 10);
    
    EXPECT_EQ(rgraph[0], R4Group({rr_xx_xx, rr_xy_xy, rr_yy_yy,}));
    
    EXPECT_EQ(rgraph[1], R4Group({rr_x_xx_1, rr_y_xy_1, rr_y_yy_1,}));
    
    EXPECT_EQ(rgraph[2], R4Group({rr_x_xx_0, rr_y_xy_0, rr_y_yy_0,}));
    
    EXPECT_EQ(rgraph[3], R4Group({rr_x_x_1, rr_y_y_1,}));
    
    EXPECT_EQ(rgraph[4], R4Group({rd_0_xx_2, rd_0_xy_2, rd_0_yy_2}));
    
    EXPECT_EQ(rgraph[5], R4Group({rd_0_xx_1, rd_0_xy_1, rd_0_yy_1}));
    
    EXPECT_EQ(rgraph[6], R4Group({rd_0_xx_0, rd_0_xy_0, rd_0_yy_0}));
    
    EXPECT_EQ(rgraph[7], R4Group({rd_0_x_2, rd_0_y_2,}));
    
    EXPECT_EQ(rgraph[8], R4Group({rd_0_x_1, rd_0_y_1,}));
    
    EXPECT_EQ(rgraph[9], R4Group({rd_0_0_2,}));
    
    EXPECT_EQ(rints, sints);
    
    EXPECT_EQ(rgraph.edge(0), std::set<int>({1, 2, 3, 5, 6}));
    
    EXPECT_EQ(rgraph.edge(1), std::set<int>({4, 7}));
    
    EXPECT_EQ(rgraph.edge(2), std::set<int>({8,}));
    
    EXPECT_EQ(rgraph.edge(3), std::set<int>({7, 8, 9}));
    
    EXPECT_EQ(rgraph.edge(4), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(5), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(6), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(7), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(8), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(9), std::set<int>({}));
}

TEST_F(EriDriverTest, ApplyKetVrrWithGraphForSD)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);
    
    // bra and ket pairs
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_0_xx = T2CPair({"GC", "GD"}, {s_0, d_xx});
    
    const auto k_0_xy = T2CPair({"GC", "GD"}, {s_0, d_xy});
    
    const auto k_0_yy = T2CPair({"GC", "GD"}, {s_0, d_yy});
    
    const auto k_0_x = T2CPair({"GC", "GD"}, {s_0, p_x});
    
    const auto k_0_y = T2CPair({"GC", "GD"}, {s_0, p_y});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto t_0_xx = T4CIntegral(b_0_0, k_0_xx, operi);
    
    const auto t_0_xy = T4CIntegral(b_0_0, k_0_xy, operi);
    
    const auto t_0_yy = T4CIntegral(b_0_0, k_0_yy, operi);
    
    const auto t_0_x_0 = T4CIntegral(b_0_0, k_0_x, operi);
    
    const auto t_0_y_0 = T4CIntegral(b_0_0, k_0_y, operi);
    
    const auto t_0_x_1 = T4CIntegral(b_0_0, k_0_x, operi, 1);
    
    const auto t_0_y_1 = T4CIntegral(b_0_0, k_0_y, operi, 1);
    
    const auto t_0_0_0 = T4CIntegral(b_0_0, k_0_0, operi);
    
    const auto t_0_0_1 = T4CIntegral(b_0_0, k_0_0, operi, 1);
    
    const auto t_0_0_2 = T4CIntegral(b_0_0, k_0_0, operi, 2);
    
    // generate graph
    
    const auto rd_0_xx =  R4CDist(R4CTerm(t_0_xx));
    
    const auto rd_0_xy =  R4CDist(R4CTerm(t_0_xy));
    
    const auto rd_0_yy =  R4CDist(R4CTerm(t_0_yy));
    
    R4Graph rgraph(R4Group({rd_0_xx, rd_0_xy, rd_0_yy}));
    
    std::set<T4CIntegral> sints;
    
    eri_drv.apply_ket_vrr(rgraph, sints);
    
    // set up reference data
    
    std::set<T4CIntegral> rints;
    
    const auto rr_0_yy = eri_drv.apply_ket_vrr(R4CTerm(t_0_yy), rints);
    
    const auto rr_0_xy = eri_drv.apply_ket_vrr(R4CTerm(t_0_xy), rints);
    
    const auto rr_0_xx = eri_drv.apply_ket_vrr(R4CTerm(t_0_xx), rints);
    
    const auto rr_0_y_1 = eri_drv.apply_ket_vrr(R4CTerm(t_0_y_1), rints);
    
    const auto rr_0_x_1 = eri_drv.apply_ket_vrr(R4CTerm(t_0_x_1), rints);
    
    const auto rr_0_y_0 = eri_drv.apply_ket_vrr(R4CTerm(t_0_y_0), rints);
    
    const auto rr_0_x_0 = eri_drv.apply_ket_vrr(R4CTerm(t_0_x_0), rints);

    // compare vertices and edges of graph
    
    EXPECT_EQ(rgraph.vertices(), 6);
    
    EXPECT_EQ(rgraph[0], R4Group({rr_0_xx, rr_0_xy, rr_0_yy,}));
    
    EXPECT_EQ(rgraph[1], R4Group({rr_0_x_1, rr_0_y_1}));
    
    EXPECT_EQ(rgraph[2], R4Group({rr_0_x_0, rr_0_y_0}));
    
    EXPECT_EQ(rgraph[3], R4Group({R4CTerm(t_0_0_2), }));
    
    EXPECT_EQ(rgraph[4], R4Group({R4CTerm(t_0_0_1), }));
    
    EXPECT_EQ(rgraph[5], R4Group({R4CTerm(t_0_0_0), }));
    
    EXPECT_EQ(sints, rints);
    
    EXPECT_EQ(rgraph.edge(0), std::set<int>({1, 2, 4, 5}));
    
    EXPECT_EQ(rgraph.edge(1), std::set<int>({3, 4}));
    
    EXPECT_EQ(rgraph.edge(2), std::set<int>({4, 5}));
    
    EXPECT_EQ(rgraph.edge(3), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(4), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(5), std::set<int>({}));
}

TEST_F(EriDriverTest, ApplyRecursionPPPP)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);

    // bra and ket pairs
    
    const auto b_x_x = T2CPair({"GA", "GB"}, {p_x, p_x});
    
    const auto b_x_y = T2CPair({"GA", "GB"}, {p_x, p_y});
    
    const auto b_y_y = T2CPair({"GA", "GB"}, {p_y, p_y});
    
    const auto b_0_xx = T2CPair({"GA", "GB"}, {s_0, d_xx});
    
    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto b_0_x = T2CPair({"GA", "GB"}, {s_0, p_x});
    
    const auto b_0_y = T2CPair({"GA", "GB"}, {s_0, p_y});
    
    const auto b_0_0 = T2CPair({"GA", "GB"}, {s_0, s_0});
    
    const auto k_x_x = T2CPair({"GC", "GD"}, {p_x, p_x});
    
    const auto k_x_y = T2CPair({"GC", "GD"}, {p_x, p_y});
    
    const auto k_y_y = T2CPair({"GC", "GD"}, {p_y, p_y});
    
    const auto k_0_xx = T2CPair({"GC", "GD"}, {s_0, d_xx});
    
    const auto k_0_xy = T2CPair({"GC", "GD"}, {s_0, d_xy});
    
    const auto k_0_yy = T2CPair({"GC", "GD"}, {s_0, d_yy});
    
    const auto k_0_x = T2CPair({"GC", "GD"}, {s_0, p_x});
    
    const auto k_0_y = T2CPair({"GC", "GD"}, {s_0, p_y});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto t_x_x_x_x = T4CIntegral(b_x_x, k_x_x, operi);
    
    const auto t_x_y_x_y = T4CIntegral(b_x_y, k_x_y, operi);
    
    const auto t_y_y_y_y = T4CIntegral(b_y_y, k_y_y, operi);
    
    const auto t_0_xx_x_x = T4CIntegral(b_0_xx, k_x_x, operi);
    
    const auto t_0_xy_x_y = T4CIntegral(b_0_xy, k_x_y, operi);
    
    const auto t_0_yy_y_y = T4CIntegral(b_0_yy, k_y_y, operi);
    
    const auto t_0_xx_0_xx = T4CIntegral(b_0_xx, k_0_xx, operi);
    
    const auto t_0_xy_0_xy = T4CIntegral(b_0_xy, k_0_xy, operi);
    
    const auto t_0_yy_0_yy = T4CIntegral(b_0_yy, k_0_yy, operi);
    
    const auto t_0_xx_0_x = T4CIntegral(b_0_xx, k_0_x, operi);
    
    const auto t_0_xy_0_y = T4CIntegral(b_0_xy, k_0_y, operi);
    
    const auto t_0_yy_0_y = T4CIntegral(b_0_yy, k_0_y, operi);
    
    const auto t_0_x_x_x = T4CIntegral(b_0_x, k_x_x, operi);
    
    const auto t_0_y_x_y = T4CIntegral(b_0_y, k_x_y, operi);
    
    const auto t_0_y_y_y = T4CIntegral(b_0_y, k_y_y, operi);
    
    const auto t_0_x_0_xx_1 = T4CIntegral(b_0_x, k_0_xx, operi, 1);
    
    const auto t_0_y_0_xy_1 = T4CIntegral(b_0_y, k_0_xy, operi, 1);
    
    const auto t_0_y_0_yy_1 = T4CIntegral(b_0_y, k_0_yy, operi, 1);
    
    const auto t_0_x_0_xx = T4CIntegral(b_0_x, k_0_xx, operi);
    
    const auto t_0_y_0_xy = T4CIntegral(b_0_y, k_0_xy, operi);
    
    const auto t_0_y_0_yy = T4CIntegral(b_0_y, k_0_yy, operi);
    
    const auto t_0_x_0_x_1 = T4CIntegral(b_0_x, k_0_x, operi, 1);
    
    const auto t_0_y_0_y_1 = T4CIntegral(b_0_y, k_0_y, operi, 1);
    
    const auto t_0_x_0_x = T4CIntegral(b_0_x, k_0_x, operi);
    
    const auto t_0_y_0_y = T4CIntegral(b_0_y, k_0_y, operi);
    
    const auto t_0_x_0_0_1 = T4CIntegral(b_0_x, k_0_0, operi, 1);
    
    const auto t_0_y_0_0_1 = T4CIntegral(b_0_y, k_0_0, operi, 1);
    
    const auto t_0_0_0_xx_2 = T4CIntegral(b_0_0, k_0_xx, operi, 2);
    
    const auto t_0_0_0_xy_2 = T4CIntegral(b_0_0, k_0_xy, operi, 2);
    
    const auto t_0_0_0_yy_2 = T4CIntegral(b_0_0, k_0_yy, operi, 2);
    
    const auto t_0_0_0_xx_1 = T4CIntegral(b_0_0, k_0_xx, operi, 1);
    
    const auto t_0_0_0_xy_1 = T4CIntegral(b_0_0, k_0_xy, operi, 1);
    
    const auto t_0_0_0_yy_1 = T4CIntegral(b_0_0, k_0_yy, operi, 1);
    
    const auto t_0_0_0_xx = T4CIntegral(b_0_0, k_0_xx, operi);
    
    const auto t_0_0_0_xy = T4CIntegral(b_0_0, k_0_xy, operi);
    
    const auto t_0_0_0_yy = T4CIntegral(b_0_0, k_0_yy, operi);
    
    const auto t_0_0_0_x_3 = T4CIntegral(b_0_0, k_0_x, operi, 3);
    
    const auto t_0_0_0_y_3 = T4CIntegral(b_0_0, k_0_y, operi, 3);
    
    const auto t_0_0_0_x_2 = T4CIntegral(b_0_0, k_0_x, operi, 2);
    
    const auto t_0_0_0_y_2 = T4CIntegral(b_0_0, k_0_y, operi, 2);
    
    const auto t_0_0_0_x_1 = T4CIntegral(b_0_0, k_0_x, operi, 1);
    
    const auto t_0_0_0_y_1 = T4CIntegral(b_0_0, k_0_y, operi, 1);
    
    const auto t_0_0_0_x = T4CIntegral(b_0_0, k_0_x, operi);
    
    const auto t_0_0_0_y = T4CIntegral(b_0_0, k_0_y, operi);
    
    const auto t_0_0_0_0_4 = T4CIntegral(b_0_0, k_0_0, operi, 4);
    
    const auto t_0_0_0_0_3 = T4CIntegral(b_0_0, k_0_0, operi, 3);
    
    const auto t_0_0_0_0_2 = T4CIntegral(b_0_0, k_0_0, operi, 2);
    
    const auto t_0_0_0_0_1 = T4CIntegral(b_0_0, k_0_0, operi, 1);
    
    const auto t_0_0_0_0_0 = T4CIntegral(b_0_0, k_0_0, operi);
    
    // generate graph
    
    const auto rd_x_x_x_x = R4CDist(R4CTerm(t_x_x_x_x));
    
    const auto rd_x_y_x_y = R4CDist(R4CTerm(t_x_y_x_y));
    
    const auto rd_y_y_y_y = R4CDist(R4CTerm(t_y_y_y_y));
    
    R4Graph rgraph(R4Group({rd_x_x_x_x, rd_x_y_x_y,
                            rd_y_y_y_y}));
    
    std::set<T4CIntegral> sints;
    
    eri_drv.apply_recursion(rgraph, sints);
    
    // set up reference data
    
    std::set<T4CIntegral> rints;
    
    const auto rr_y_y_y_y = eri_drv.apply_bra_hrr(R4CTerm(t_y_y_y_y), rints);
    
    const auto rr_x_y_x_y = eri_drv.apply_bra_hrr(R4CTerm(t_x_y_x_y), rints);
    
    const auto rr_x_x_x_x = eri_drv.apply_bra_hrr(R4CTerm(t_x_x_x_x), rints);
    
    const auto rr_0_yy_y_y = eri_drv.apply_ket_hrr(R4CTerm(t_0_yy_y_y), rints);
    
    const auto rr_0_xy_x_y = eri_drv.apply_ket_hrr(R4CTerm(t_0_xy_x_y), rints);
    
    const auto rr_0_xx_x_x = eri_drv.apply_ket_hrr(R4CTerm(t_0_xx_x_x), rints);
    
    const auto rr_0_yy_0_yy = eri_drv.apply_bra_vrr(R4CTerm(t_0_yy_0_yy), rints);
    
    const auto rr_0_xy_0_xy = eri_drv.apply_bra_vrr(R4CTerm(t_0_xy_0_xy), rints);
    
    const auto rr_0_xx_0_xx = eri_drv.apply_bra_vrr(R4CTerm(t_0_xx_0_xx), rints);
    
    const auto rr_0_yy_0_y = eri_drv.apply_bra_vrr(R4CTerm(t_0_yy_0_y), rints);
    
    const auto rr_0_xy_0_y = eri_drv.apply_bra_vrr(R4CTerm(t_0_xy_0_y), rints);
    
    const auto rr_0_xx_0_x = eri_drv.apply_bra_vrr(R4CTerm(t_0_xx_0_x), rints);
    
    const auto rr_0_y_y_y = eri_drv.apply_ket_hrr(R4CTerm(t_0_y_y_y), rints);
    
    const auto rr_0_y_x_y = eri_drv.apply_ket_hrr(R4CTerm(t_0_y_x_y), rints);
    
    const auto rr_0_x_x_x = eri_drv.apply_ket_hrr(R4CTerm(t_0_x_x_x), rints);
    
    const auto rr_0_y_0_yy_1 = eri_drv.apply_bra_vrr(R4CTerm(t_0_y_0_yy_1), rints);
    
    const auto rr_0_y_0_xy_1 = eri_drv.apply_bra_vrr(R4CTerm(t_0_y_0_xy_1), rints);
    
    const auto rr_0_x_0_xx_1 = eri_drv.apply_bra_vrr(R4CTerm(t_0_x_0_xx_1), rints);
    
    const auto rr_0_y_0_yy = eri_drv.apply_bra_vrr(R4CTerm(t_0_y_0_yy), rints);
    
    const auto rr_0_y_0_xy = eri_drv.apply_bra_vrr(R4CTerm(t_0_y_0_xy), rints);
    
    const auto rr_0_x_0_xx = eri_drv.apply_bra_vrr(R4CTerm(t_0_x_0_xx), rints);
    
    const auto rr_0_y_0_y_1 = eri_drv.apply_bra_vrr(R4CTerm(t_0_y_0_y_1), rints);
    
    const auto rr_0_x_0_x_1 = eri_drv.apply_bra_vrr(R4CTerm(t_0_x_0_x_1), rints);
    
    const auto rr_0_y_0_y = eri_drv.apply_bra_vrr(R4CTerm(t_0_y_0_y), rints);
    
    const auto rr_0_x_0_x = eri_drv.apply_bra_vrr(R4CTerm(t_0_x_0_x), rints);
    
    const auto rr_0_y_0_0_1 = eri_drv.apply_bra_vrr(R4CTerm(t_0_y_0_0_1), rints);
    
    const auto rr_0_x_0_0_1 = eri_drv.apply_bra_vrr(R4CTerm(t_0_x_0_0_1), rints);
    
    const auto rr_0_0_0_yy_2 = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_yy_2), rints);
    
    const auto rr_0_0_0_xy_2 = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_xy_2), rints);
    
    const auto rr_0_0_0_xx_2 = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_xx_2), rints);
    
    const auto rr_0_0_0_yy_1 = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_yy_1), rints);
    
    const auto rr_0_0_0_xy_1 = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_xy_1), rints);
    
    const auto rr_0_0_0_xx_1 = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_xx_1), rints);
    
    const auto rr_0_0_0_yy = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_yy), rints);
    
    const auto rr_0_0_0_xy = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_xy), rints);
    
    const auto rr_0_0_0_xx = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_xx), rints);
    
    const auto rr_0_0_0_y_3 = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_y_3), rints);
    
    const auto rr_0_0_0_x_3 = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_x_3), rints);
    
    const auto rr_0_0_0_y_2 = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_y_2), rints);
    
    const auto rr_0_0_0_x_2 = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_x_2), rints);
    
    const auto rr_0_0_0_y_1 = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_y_1), rints);
    
    const auto rr_0_0_0_x_1 = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_x_1), rints);
    
    const auto rr_0_0_0_y = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_y), rints);
    
    const auto rr_0_0_0_x = eri_drv.apply_ket_vrr(R4CTerm(t_0_0_0_x), rints);
    
    // compare vertices and edges of graph
    
    EXPECT_EQ(rgraph.vertices(), 22);
    
    EXPECT_EQ(rgraph[0], R4Group({rr_x_x_x_x, rr_x_y_x_y, rr_y_y_y_y,}));
    
    EXPECT_EQ(rgraph[1], R4Group({rr_0_xx_x_x, rr_0_xy_x_y, rr_0_yy_y_y,}));
    
    EXPECT_EQ(rgraph[2], R4Group({rr_0_xx_0_xx, rr_0_xy_0_xy, rr_0_yy_0_yy,}));
    
    EXPECT_EQ(rgraph[3], R4Group({rr_0_xx_0_x, rr_0_xy_0_y, rr_0_yy_0_y,}));
    
    EXPECT_EQ(rgraph[4], R4Group({rr_0_x_x_x, rr_0_y_x_y, rr_0_y_y_y,}));
    
    EXPECT_EQ(rgraph[5], R4Group({rr_0_x_0_xx_1, rr_0_y_0_xy_1, rr_0_y_0_yy_1,}));
    
    EXPECT_EQ(rgraph[6], R4Group({rr_0_x_0_xx, rr_0_y_0_xy, rr_0_y_0_yy,}));
    
    EXPECT_EQ(rgraph[7], R4Group({rr_0_x_0_x_1, rr_0_y_0_y_1}));
    
    EXPECT_EQ(rgraph[8], R4Group({rr_0_x_0_x, rr_0_y_0_y}));
    
    EXPECT_EQ(rgraph[9], R4Group({rr_0_x_0_0_1, rr_0_y_0_0_1}));
    
    EXPECT_EQ(rgraph[10], R4Group({rr_0_0_0_xx_2, rr_0_0_0_xy_2, rr_0_0_0_yy_2,}));
    
    EXPECT_EQ(rgraph[11], R4Group({rr_0_0_0_xx_1, rr_0_0_0_xy_1, rr_0_0_0_yy_1,}));
    
    EXPECT_EQ(rgraph[12], R4Group({rr_0_0_0_xx, rr_0_0_0_xy, rr_0_0_0_yy,}));
    
    EXPECT_EQ(rgraph[13], R4Group({rr_0_0_0_x_3, rr_0_0_0_y_3,}));
    
    EXPECT_EQ(rgraph[14], R4Group({rr_0_0_0_x_2, rr_0_0_0_y_2,}));
    
    EXPECT_EQ(rgraph[15], R4Group({rr_0_0_0_x_1, rr_0_0_0_y_1,}));
    
    EXPECT_EQ(rgraph[16], R4Group({rr_0_0_0_x, rr_0_0_0_y,}));
    
    EXPECT_EQ(rgraph[17], R4Group({R4CTerm(t_0_0_0_0_4),}));
    
    EXPECT_EQ(rgraph[18], R4Group({R4CTerm(t_0_0_0_0_3), }));
    
    EXPECT_EQ(rgraph[19], R4Group({R4CTerm(t_0_0_0_0_2), }));
    
    EXPECT_EQ(rgraph[20], R4Group({R4CTerm(t_0_0_0_0_1), }));
    
    EXPECT_EQ(rgraph[21], R4Group({R4CTerm(t_0_0_0_0_0), }));
    
    EXPECT_EQ(sints, rints);
    
    EXPECT_EQ(rgraph.edge(0), std::set<int>({1, 4}));
    
    EXPECT_EQ(rgraph.edge(1), std::set<int>({2, 3}));
    
    EXPECT_EQ(rgraph.edge(2), std::set<int>({5, 6, 7, 11, 12}));
    
    EXPECT_EQ(rgraph.edge(3), std::set<int>({7, 8, 9, 15, 16}));
    
    EXPECT_EQ(rgraph.edge(4), std::set<int>({6, 8}));
    
    EXPECT_EQ(rgraph.edge(5), std::set<int>({10, 14}));
    
    EXPECT_EQ(rgraph.edge(6), std::set<int>({11, 12, 15}));
    
    EXPECT_EQ(rgraph.edge(7), std::set<int>({14, 15, 19}));
    
    EXPECT_EQ(rgraph.edge(8), std::set<int>({15, 16, 20}));
    
    EXPECT_EQ(rgraph.edge(9), std::set<int>({19, 20}));
    
    EXPECT_EQ(rgraph.edge(10), std::set<int>({13, 14, 18, 19}));
    
    EXPECT_EQ(rgraph.edge(11), std::set<int>({14, 15, 19, 20}));
    
    EXPECT_EQ(rgraph.edge(12), std::set<int>({15, 16, 20, 21}));
    
    EXPECT_EQ(rgraph.edge(13), std::set<int>({17, 18}));
    
    EXPECT_EQ(rgraph.edge(14), std::set<int>({18, 19}));
    
    EXPECT_EQ(rgraph.edge(15), std::set<int>({19, 20}));
    
    EXPECT_EQ(rgraph.edge(16), std::set<int>({20, 21}));
    
    EXPECT_EQ(rgraph.edge(17), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(18), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(19), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(20), std::set<int>({}));
    
    EXPECT_EQ(rgraph.edge(21), std::set<int>({}));
}

TEST_F(EriDriverTest, CreateGraphWithDiagonal)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto p_z = TensorComponent(0, 0, 1);

    // bra and ket pairs
    
    const auto b_0_x = T2CPair({"GA", "GB"}, {s_0, p_x});
    
    const auto b_0_y = T2CPair({"GA", "GB"}, {s_0, p_y});
    
    const auto b_0_z = T2CPair({"GA", "GB"}, {s_0, p_z});
    
    const auto k_0_x = T2CPair({"GC", "GD"}, {s_0, p_x});
    
    const auto k_0_y = T2CPair({"GC", "GD"}, {s_0, p_y});
    
    const auto k_0_z = T2CPair({"GC", "GD"}, {s_0, p_z});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto t_0_x_0_x = T4CIntegral(b_0_x, k_0_x, operi);
    
    const auto t_0_y_0_y = T4CIntegral(b_0_y, k_0_y, operi);
    
    const auto t_0_z_0_z = T4CIntegral(b_0_z, k_0_z, operi);
    
    // reference generate graph
    
    const auto rd_0_x_0_x = R4CDist(R4CTerm(t_0_x_0_x));
    
    const auto rd_0_y_0_y = R4CDist(R4CTerm(t_0_y_0_y));
    
    const auto rd_0_z_0_z = R4CDist(R4CTerm(t_0_z_0_z));
    
    R4Graph rgraph(R4Group({rd_0_x_0_x, rd_0_y_0_y, rd_0_z_0_z}));
    
    std::set<T4CIntegral> sints;
    
    eri_drv.apply_recursion(rgraph, sints);
    
    // check create graph
    
    const auto tgraph = eri_drv.create_graph(0, 1, 0, 1, true);
    
    EXPECT_EQ(rgraph, tgraph);
}

TEST_F(EriDriverTest, CreateGraph)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto p_z = TensorComponent(0, 0, 1);

    // bra and ket pairs
    
    const auto b_0_x = T2CPair({"GA", "GB"}, {s_0, p_x});
    
    const auto b_0_y = T2CPair({"GA", "GB"}, {s_0, p_y});
    
    const auto b_0_z = T2CPair({"GA", "GB"}, {s_0, p_z});
    
    const auto k_0_x = T2CPair({"GC", "GD"}, {s_0, p_x});
    
    const auto k_0_y = T2CPair({"GC", "GD"}, {s_0, p_y});
    
    const auto k_0_z = T2CPair({"GC", "GD"}, {s_0, p_z});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto t_0_x_0_x = T4CIntegral(b_0_x, k_0_x, operi);
    
    const auto t_0_x_0_y = T4CIntegral(b_0_x, k_0_y, operi);
    
    const auto t_0_x_0_z = T4CIntegral(b_0_x, k_0_z, operi);
    
    const auto t_0_y_0_x = T4CIntegral(b_0_y, k_0_x, operi);
    
    const auto t_0_y_0_y = T4CIntegral(b_0_y, k_0_y, operi);
    
    const auto t_0_y_0_z = T4CIntegral(b_0_y, k_0_z, operi);
    
    const auto t_0_z_0_x = T4CIntegral(b_0_z, k_0_x, operi);
    
    const auto t_0_z_0_y = T4CIntegral(b_0_z, k_0_y, operi);
    
    const auto t_0_z_0_z = T4CIntegral(b_0_z, k_0_z, operi);
    
    // reference generate graph
    
    const auto rd_0_x_0_x = R4CDist(R4CTerm(t_0_x_0_x));
    
    const auto rd_0_x_0_y = R4CDist(R4CTerm(t_0_x_0_y));
    
    const auto rd_0_x_0_z = R4CDist(R4CTerm(t_0_x_0_z));
    
    const auto rd_0_y_0_x = R4CDist(R4CTerm(t_0_y_0_x));
    
    const auto rd_0_y_0_y = R4CDist(R4CTerm(t_0_y_0_y));
    
    const auto rd_0_y_0_z = R4CDist(R4CTerm(t_0_y_0_z));
    
    const auto rd_0_z_0_x = R4CDist(R4CTerm(t_0_z_0_x));
    
    const auto rd_0_z_0_y = R4CDist(R4CTerm(t_0_z_0_y));
    
    const auto rd_0_z_0_z = R4CDist(R4CTerm(t_0_z_0_z));
    
    R4Graph rgraph(R4Group({rd_0_x_0_x, rd_0_x_0_y, rd_0_x_0_z,
                            rd_0_y_0_x, rd_0_y_0_y, rd_0_y_0_z,
                            rd_0_z_0_x, rd_0_z_0_y, rd_0_z_0_z,}));
    
    std::set<T4CIntegral> sints;
    
    eri_drv.apply_recursion(rgraph, sints);
    
    // check create graph
    
    const auto tgraph = eri_drv.create_graph(0, 1, 0, 1, false);
    
    EXPECT_EQ(rgraph, tgraph);
}

TEST_F(EriDriverTest, CreateGraphs)
{
    EriDriver eri_drv;
    
    const auto vgraphs = eri_drv.create_graphs(1, 1, 1, 1, false);
    
    EXPECT_EQ(vgraphs[0], eri_drv.create_graph(0, 0, 0, 0, false));
    
    EXPECT_EQ(vgraphs[1], eri_drv.create_graph(0, 0, 0, 1, false));
    
    EXPECT_EQ(vgraphs[2], eri_drv.create_graph(0, 0, 1, 0, false));
    
    EXPECT_EQ(vgraphs[3], eri_drv.create_graph(0, 0, 1, 1, false));
    
    EXPECT_EQ(vgraphs[4], eri_drv.create_graph(0, 1, 0, 0, false));
    
    EXPECT_EQ(vgraphs[5], eri_drv.create_graph(0, 1, 0, 1, false));
    
    EXPECT_EQ(vgraphs[6], eri_drv.create_graph(0, 1, 1, 0, false));
    
    EXPECT_EQ(vgraphs[7], eri_drv.create_graph(0, 1, 1, 1, false));
    
    EXPECT_EQ(vgraphs[8], eri_drv.create_graph(1, 0, 0, 0, false));
    
    EXPECT_EQ(vgraphs[9], eri_drv.create_graph(1, 0, 0, 1, false));
    
    EXPECT_EQ(vgraphs[10], eri_drv.create_graph(1, 0, 1, 0, false));
    
    EXPECT_EQ(vgraphs[11], eri_drv.create_graph(1, 0, 1, 1, false));
    
    EXPECT_EQ(vgraphs[12], eri_drv.create_graph(1, 1, 0, 0, false));
    
    EXPECT_EQ(vgraphs[13], eri_drv.create_graph(1, 1, 0, 1, false));
    
    EXPECT_EQ(vgraphs[14], eri_drv.create_graph(1, 1, 1, 0, false));
    
    EXPECT_EQ(vgraphs[15], eri_drv.create_graph(1, 1, 1, 1, false));
}

TEST_F(EriDriverTest, GraphSignaturesMap)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);
    
    // bra and ket pairs
    
    const auto b_x_x = T2CPair({"GA", "GB"}, {p_x, p_x});
    
    const auto b_x_y = T2CPair({"GA", "GB"}, {p_x, p_y});
    
    const auto b_y_y = T2CPair({"GA", "GB"}, {p_y, p_y});
    
    const auto b_0_xx = T2CPair({"GA", "GB"}, {s_0, d_xx});
    
    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto b_0_x = T2CPair({"GA", "GB"}, {s_0, p_x});
    
    const auto b_0_y = T2CPair({"GA", "GB"}, {s_0, p_y});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto t_x_x = T4CIntegral(b_x_x, k_0_0, operi);
    
    const auto t_x_y = T4CIntegral(b_x_y, k_0_0, operi);
    
    const auto t_y_y = T4CIntegral(b_y_y, k_0_0, operi);
    
    const auto t_0_xx = T4CIntegral(b_0_xx, k_0_0, operi);
    
    const auto t_0_xy = T4CIntegral(b_0_xy, k_0_0, operi);
    
    const auto t_0_yy = T4CIntegral(b_0_yy, k_0_0, operi);
    
    const auto t_0_x = T4CIntegral(b_0_x, k_0_0, operi);
    
    const auto t_0_y = T4CIntegral(b_0_y, k_0_0, operi);
    
    // generate graph
    
    const auto rd_x_x = R4CDist(R4CTerm(t_x_x));
    
    const auto rd_x_y = R4CDist(R4CTerm(t_x_y));
    
    const auto rd_y_y = R4CDist(R4CTerm(t_y_y));
    
    R4Graph rgraph(R4Group({rd_x_x, rd_x_y, rd_y_y}));
    
    std::set<T4CIntegral> sints;
    
    eri_drv.apply_bra_hrr(rgraph, sints);
    
    // compare signature maps
    
    const auto smap = rgraph.signatures<T4CIntegral>();
    
    std::map<Signature<T4CIntegral>, R4Group> rmap;
    
    for (int i = 0; i < 3; i++)
    {
        rmap[rgraph[i].signature()] = rgraph[i];
    }
    
    EXPECT_EQ(rmap, smap);
}

TEST_F(EriDriverTest, RepositoryAdd)
{
    EriDriver eri_drv;
    
    // tensor components
    
    const auto s_0 = TensorComponent(0, 0, 0);
    
    const auto p_x = TensorComponent(1, 0, 0);
    
    const auto p_y = TensorComponent(0, 1, 0);
    
    const auto d_xx = TensorComponent(2, 0, 0);
    
    const auto d_xy = TensorComponent(1, 1, 0);
    
    const auto d_yy = TensorComponent(0, 2, 0);
    
    // bra and ket pairs
    
    const auto b_x_x = T2CPair({"GA", "GB"}, {p_x, p_x});
    
    const auto b_x_y = T2CPair({"GA", "GB"}, {p_x, p_y});
    
    const auto b_y_y = T2CPair({"GA", "GB"}, {p_y, p_y});
    
    const auto b_0_xx = T2CPair({"GA", "GB"}, {s_0, d_xx});
    
    const auto b_0_xy = T2CPair({"GA", "GB"}, {s_0, d_xy});
    
    const auto b_0_yy = T2CPair({"GA", "GB"}, {s_0, d_yy});
    
    const auto b_0_x = T2CPair({"GA", "GB"}, {s_0, p_x});
    
    const auto b_0_y = T2CPair({"GA", "GB"}, {s_0, p_y});
    
    const auto k_0_0 = T2CPair({"GC", "GD"}, {s_0, s_0});
    
    // operator
    
    const auto operi = OperatorComponent("1/|r-r'|");
    
    // integral components
    
    const auto t_x_x = T4CIntegral(b_x_x, k_0_0, operi);
    
    const auto t_x_y = T4CIntegral(b_x_y, k_0_0, operi);
    
    const auto t_y_y = T4CIntegral(b_y_y, k_0_0, operi);
    
    const auto t_0_xx = T4CIntegral(b_0_xx, k_0_0, operi);
    
    const auto t_0_xy = T4CIntegral(b_0_xy, k_0_0, operi);
    
    const auto t_0_yy = T4CIntegral(b_0_yy, k_0_0, operi);
    
    const auto t_0_x = T4CIntegral(b_0_x, k_0_0, operi);
    
    const auto t_0_y = T4CIntegral(b_0_y, k_0_0, operi);
    
    // generate graph
    
    const auto rd_x_x = R4CDist(R4CTerm(t_x_x));
    
    const auto rd_x_y = R4CDist(R4CTerm(t_x_y));
    
    const auto rd_y_y = R4CDist(R4CTerm(t_y_y));
    
    R4Graph rgraph(R4Group({rd_x_x, rd_x_y, rd_y_y}));
    
    std::set<T4CIntegral> sints;
    
    eri_drv.apply_bra_hrr(rgraph, sints);
    
    // create repository for integrals
    
    Repository<R4Group, T4CIntegral> repo;
    
    repo.add(V4Graphs({rgraph,}));
    
    const auto smap = rgraph.signatures<T4CIntegral>();
    
    const auto ref_repo = Repository<R4Group, T4CIntegral>(V4Graphs({rgraph,}), smap); 
    
    EXPECT_EQ(repo, ref_repo);
    
//    // compare signature maps
//
//    const auto smap = rgraph.signatures<T4CIntegral>();
//
//    std::map<Signature<T4CIntegral>, R4Group> rmap;
//
//    for (int i = 0; i < 3; i++)
//    {
//        rmap[rgraph[i].signature()] = rgraph[i];
//    }
//
//    EXPECT_EQ(rmap, smap);
}

//const auto nverts = rgraph.vertices();
//
//std::set<I4CIntegral> gints;
//
//for (int i = 0; i < nverts; i++)
//{
//    std::cout << "Vertice: " << i << std::endl;
//
//    for (const auto& tval : rgraph[i].roots())
//    {
//        auto gint = I4CIntegral(tval.integral());
//        
//        gints.insert(gint);
//        
//        std::cout << tval.integral().label(true) << "->" << gint.label(true) << " ";
//    }
//    
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
//
//for (const auto& tval : gints)
//{
//    std::cout << tval.label(true) << " ";
//}
