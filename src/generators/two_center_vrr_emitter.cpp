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

#include "two_center_vrr_emitter.hpp"

#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "operator.hpp"
#include "spherical_harmonics.hpp"
#include "t2c_defs.hpp"
#include "t2c_ovl_driver.hpp"
#include "tensor.hpp"

namespace {  // two-center VRR emitter helpers

/// The lowercase spectroscopic shell label of an angular momentum (s, p, d, ...).
std::string
shell_label(const int l)
{
    static const std::string labels = "spdfghiklmn";

    return std::string(1, labels[l]);
}

/// The number of Cartesian components of a shell of angular momentum l.
int
cartesian_count(const int l)
{
    return (l + 1) * (l + 2) / 2;
}

/// The canonical index of a Cartesian component within its shell.
int
component_index(const int l, const TensorComponent& component)
{
    int index = 0;

    for (const auto& c : Tensor(l).components())
    {
        if (c == component) return index;

        index++;
    }

    return -1;
}

/// The signature "void compute_<lb>(...)" of the single-step Cartesian VRR kernel
/// (no body, no terminator): the pair, the lower (s|lb-1)/(s|lb-2) integrals, the
/// Pc distances, and the Cartesian target out-parameter.
std::string
signature_cartesian_text(const int lb)
{
    const auto name = "compute_" + shell_label(lb);

    const auto indent = std::string(name.size() + 6, ' ');  // align under "void name("

    std::ostringstream os;

    os << "void " << name << "(const osfunc::CBasisFunctionPair& pair,\n";
    os << indent << "const osfunc::CArray<double>& s" << shell_label(lb - 1) << ",\n";

    if (lb >= 2) os << indent << "const osfunc::CArray<double>& s" << shell_label(lb - 2) << ",\n";

    os << indent << "const osfunc::CArray<double>& pc,\n";
    os << indent << "      osfunc::CArray<double>& s" << shell_label(lb) << ")";

    return os.str();
}

/// The signature "void compute_<lb>_sph(...)" of the full spherical VRR kernel (no
/// body, no terminator): the pair, the (s|s) seed, the Pc distances, and the
/// spherical target out-parameter.
std::string
signature_spherical_text(const int lb)
{
    const auto name = "compute_" + shell_label(lb) + "_sph";

    const auto indent = std::string(name.size() + 6, ' ');  // align under "void name("

    std::ostringstream os;

    os << "void " << name << "(const osfunc::CBasisFunctionPair& pair,\n";
    os << indent << "const osfunc::CArray<double>& ss,\n";
    os << indent << "const osfunc::CArray<double>& pc,\n";
    os << indent << "      osfunc::CArray<double>& s" << shell_label(lb) << ")";

    return os.str();
}

/// Whether a rational terminates as a decimal (denominator only 2s and 5s).
bool
is_terminating(const Fraction& value)
{
    long long den = value.denominator();

    while (den % 2 == 0) den /= 2;

    while (den % 5 == 0) den /= 5;

    return den == 1;
}

/// The exact terminating decimal of a non-negative rational, e.g. 3/8 -> "0.375".
std::string
terminating_decimal(const Fraction& value)
{
    const long long num = value.numerator();

    const long long den = value.denominator();

    if (den == 1) return std::to_string(num) + ".0";

    long long residue = den;

    int twos = 0, fives = 0;

    while (residue % 2 == 0) { residue /= 2; twos++; }

    while (residue % 5 == 0) { residue /= 5; fives++; }

    const int places = std::max(twos, fives);

    long long pow10 = 1;

    for (int i = 0; i < places; i++) pow10 *= 10;

    auto digits = std::to_string(num * (pow10 / den));

    if (static_cast<int>(digits.size()) <= places)
    {
        digits = std::string(places + 1 - digits.size(), '0') + digits;
    }

    auto fraction = digits.substr(digits.size() - places);

    while (fraction.size() > 1 && fraction.back() == '0') fraction.pop_back();

    return digits.substr(0, digits.size() - places) + "." + fraction;
}

/// The hoisted-constant name of a non-terminating rational, e.g. 1/3 -> "q1_3".
std::string
rational_name(const Fraction& value)
{
    return "q" + std::to_string(value.numerator()) + "_" + std::to_string(value.denominator());
}

/// The positive magnitude of a rational.
Fraction
magnitude(const Fraction& value)
{
    return Fraction((value.numerator() < 0) ? -value.numerator() : value.numerator(),
                    value.denominator());
}

/// Fully reduces a Cartesian (s|lb) component to the (s|s) seed by repeated ket
/// vertical recurrence, accumulating the Pc and fe factors and the prefactor.
std::vector<R2CTerm>
fully_reduce_ket(const T2COverlapDriver& drv, const R2CTerm& start)
{
    std::vector<R2CTerm> base;

    std::vector<R2CTerm> work{start};

    while (!work.empty())
    {
        std::vector<R2CTerm> next;

        for (const auto& term : work)
        {
            if (term[1].order() == 0)
            {
                base.push_back(term);

                continue;
            }

            const auto dist = drv.apply_ket_vrr(term);

            for (std::size_t j = 0; j < dist.terms(); j++) next.push_back(dist[j]);
        }

        work = next;
    }

    return base;
}

/// One contribution to a spherical (s|lb) row: a signed coefficient, a product of
/// Pc-distance pointers, and a power of the (scalar) fe factor. The (s|s) seed and
/// the contraction coefficient are common to every contribution and factored out.
struct Contribution
{
    Fraction                 coeff;     // transform * recurrence coefficient
    int                      radicand;  // square-free radical of the transform coefficient
    std::vector<std::string> pc;        // Pc-distance pointers, e.g. {"pc_x", "pc_x"}
    int                      fe_power;  // power of the fe scalar
};

/// The Pc-distance pointer of a PB factor, mapping the driver's "pb_x" label onto
/// the kernel's generic "pc_x".
std::string
pc_pointer(const Factor& factor)
{
    auto label = factor.label();  // "pb_x" / "pb_y" / "pb_z"

    label.replace(0, 2, "pc");

    return label;
}

/// The factor part of a contribution: the fe powers and the Pc product, with no
/// coefficient, e.g. "fe * pc_x[i] * pc_y[i]" (the coefficient is emitted by the
/// caller, the seed t_ss[i] appended afterwards).
std::string
factor_body(const Contribution& c)
{
    std::string body;

    for (int n = 0; n < c.fe_power; n++) body += (body.empty() ? "" : " * ") + std::string("fe");

    for (const auto& a : c.pc) body += (body.empty() ? "" : " * ") + a + "[i]";

    if (body.empty()) body = "1.0";  // pure constant (degenerate)

    return body;
}

/// The positive greatest common divisor of a component's rational coefficients.
/// Every coefficient is an integer multiple of it, so dividing it out leaves
/// integer inner coefficients.
Fraction
group_gcd(const std::vector<Contribution>& terms)
{
    long long gnum = 0, gden = 1;

    for (const auto& t : terms)
    {
        const auto mag = magnitude(t.coeff);

        gnum = std::gcd(gnum, static_cast<long long>(mag.numerator()));

        const long long d = mag.denominator();

        gden = gden / std::gcd(gden, d) * d;
    }

    if (gnum == 0) gnum = 1;

    return Fraction(static_cast<int>(gnum), static_cast<int>(gden));
}

}  // namespace

std::string
format_vrr_spherical_kernel(const int lb)
{
    const auto target = "s" + shell_label(lb);

    const auto nspher = 2 * lb + 1;

    // the full ket VRR reduction of every Cartesian (s|lb) component, keyed by the
    // ket Cartesian component.

    const T2COverlapDriver drv;

    const auto integral = I2CIntegral(I1CPair("GA", 0), I1CPair("GB", lb), Operator("1"), 0, {});

    std::map<TensorComponent, std::vector<R2CTerm>> reductions;

    for (const auto& comp : integral.components<T1CPair, T1CPair>())
    {
        reductions[comp[1]] = fully_reduce_ket(drv, R2CTerm(comp));
    }

    // fold the Cartesian->spherical transform into the recurrence: for each
    // spherical component accumulate the (Pc product, fe power) contributions.

    std::vector<std::vector<Contribution>> rows(nspher);

    for (int c = 0; c < nspher; c++)
    {
        std::map<std::pair<std::string, int>, Contribution> acc;

        for (const auto& term : sphar::two_center_spherical_factors(0, lb, 0, c))
        {
            for (const auto& rterm : reductions.at(term.ket))
            {
                std::vector<std::string> pc;

                int fe_power = 0;

                for (const auto& fact : rterm.factors())
                {
                    if (fact.name() == "PB")
                    {
                        for (int n = 0; n < rterm.factor_order(fact); n++) pc.push_back(pc_pointer(fact));
                    }
                    else if (fact.name() == "1/eta")
                    {
                        fe_power += rterm.factor_order(fact);
                    }
                }

                std::sort(pc.begin(), pc.end());

                std::string pc_key;

                for (const auto& a : pc) pc_key += a + "*";

                const auto coeff = term.factor.factor * rterm.prefactor();

                const auto key = std::make_pair(pc_key, fe_power);

                if (const auto it = acc.find(key); it != acc.end())
                    it->second.coeff = it->second.coeff + coeff;
                else
                    acc.emplace(key, Contribution{coeff, term.factor.radicand, pc, fe_power});
            }
        }

        for (const auto& [key, contrib] : acc)
        {
            if (contrib.coeff == Fraction(0)) continue;

            rows[c].push_back(contrib);
        }
    }

    // emit the kernel.

    std::ostringstream os;

    os << signature_spherical_text(lb) << "\n";
    os << "{\n";

    os << "    // bra and ket contracted basis functions\n";
    os << "    const auto bra = pair.bra();\n";
    os << "    const auto ket = pair.ket();\n\n";

    os << "    // number of atom pairs (columns)\n";
    os << "    const auto npairs = " << target << ".ncols();\n\n";

    // hoist the square-root transformation factors.
    std::set<int> radicals;

    bool uses_fe = false;

    for (const auto& row : rows)
    {
        for (const auto& contrib : row)
        {
            if (contrib.radicand != 1) radicals.insert(contrib.radicand);

            if (contrib.fe_power > 0) uses_fe = true;
        }
    }

    if (!radicals.empty())
    {
        os << "    // transformation factors\n";

        for (const auto r : radicals) os << "    const double f" << r << " = std::sqrt(" << r << ".0);\n";

        os << "\n";
    }

    // the fe = 1/(2 eta) factor only survives if the spherical transform did not
    // cancel the trace (it always does for the traceless harmonics), so the
    // exponents are pulled out only when actually needed.
    if (uses_fe)
    {
        os << "    // primitive exponents (for fe = 1 / (2 eta))\n";
        os << "    const auto bra_exps = bra.get_exponents();\n";
        os << "    const auto ket_exps = ket.get_exponents();\n\n";
    }

    // the spherical result rows are contracted, so they are set up once and
    // accumulated into over the primitive pairs (the buffer is zero-initialized).
    os << "    // spherical (s|" << shell_label(lb) << ") result rows (accumulated over primitives)\n";

    for (int c = 0; c < nspher; c++)
    {
        os << "    auto " << target << "_" << c << " = " << target << ".row(" << c << ");\n";
    }

    os << "\n";

    os << "    // outer loop over primitive basis-function pairs\n";
    os << "    for (std::size_t p = 0; p < bra.number_of_primitive_functions(); p++)\n";
    os << "    {\n";
    os << "        for (std::size_t q = 0; q < ket.number_of_primitive_functions(); q++)\n";
    os << "        {\n";

    os << "            const auto ip = p * ket.number_of_primitive_functions() + q;\n\n";

    if (uses_fe)
    {
        os << "            const auto fe = 1.0 / (2.0 * (bra_exps[p] + ket_exps[q]));   // 1 / (2 eta)\n\n";
    }

    os << "            auto t_ss = ss.row(ip);          // (s|s) primitive overlap (contraction folded in)\n";
    os << "            auto pc_x = pc.row(ip * 3 + 0);\n";
    os << "            auto pc_y = pc.row(ip * 3 + 1);\n";
    os << "            auto pc_z = pc.row(ip * 3 + 2);\n";

    // one SIMD accumulation loop per spherical component.
    for (int c = 0; c < nspher; c++)
    {
        std::set<std::string> used;

        for (const auto& contrib : rows[c])
        {
            for (const auto& a : contrib.pc) used.insert(a);
        }

        std::vector<std::string> aligned(used.begin(), used.end());

        aligned.push_back("t_ss");

        aligned.push_back(target + "_" + std::to_string(c));

        os << "\n";
        os << "            // ket spherical component " << c << "\n";
        os << "            #pragma omp simd aligned(";

        for (std::size_t n = 0; n < aligned.size(); n++) os << (n ? ", " : "") << aligned[n];

        os << " : 64)\n";

        os << "            for (std::size_t i = 0; i < npairs; i++)\n";
        os << "            {\n";

        // the (s|s) seed (with contraction folded in) is common: sd_c[i] += (...) * t_ss[i].
        // the common rational GCD and radical are factored out of the component, so
        // the inner coefficients are integers.
        const auto& terms = rows[c];

        // build the assignment prefix up to the inner sum, so continuation lines can
        // align under the first inner term.
        std::string prefix = std::string("                ") + target + "_" + std::to_string(c) + "[i] += ";

        if (terms.empty())
        {
            os << prefix << "0.0;\n";

            os << "            }\n";

            continue;
        }

        const auto g = group_gcd(terms);

        const int radicand = terms.front().radicand;

        const bool first_neg = terms.front().coeff.numerator() < 0;

        const long long gnum = first_neg ? -g.numerator() : g.numerator();

        const long long gden = g.denominator();

        // the factored-out coefficient: sign, rational GCD, radical.
        if (first_neg) prefix += "-";

        if (!(g == Fraction(1))) prefix += terminating_decimal(g) + " * ";

        if (radicand != 1) prefix += "f" + std::to_string(radicand) + " * ";

        if (terms.size() == 1)
        {
            os << prefix << factor_body(terms[0]) << " * t_ss[i];\n";
        }
        else
        {
            // one inner term per line, the operator leading each continuation line.
            const auto lead = prefix + "(";

            const auto hang = std::string(lead.size() - 2, ' ');

            os << lead;

            for (std::size_t t = 0; t < terms.size(); t++)
            {
                const auto& ct = terms[t].coeff;

                const Fraction ic(static_cast<int>(ct.numerator() * gden),
                                  static_cast<int>(ct.denominator() * gnum));

                const bool ineg = ic.numerator() < 0;

                if (t == 0)
                {
                    if (ineg) os << "-";
                }
                else
                {
                    os << "\n" << hang << (ineg ? "- " : "+ ");
                }

                const auto imag = magnitude(ic);

                if (!(imag == Fraction(1))) os << terminating_decimal(imag) << " * ";

                os << factor_body(terms[t]);
            }

            os << ") * t_ss[i];\n";
        }

        os << "            }\n";
    }

    os << "        }\n";
    os << "    }\n";
    os << "}\n";

    return os.str();
}

namespace {  // single-step Cartesian VRR helpers

/// The row-pointer name of a lower (s|.) integral component, e.g. (s|p_y) -> "sp_1".
std::string
lower_row_name(const T2CIntegral& integral)
{
    const auto l = integral[1].order();

    return "s" + shell_label(l) + "_" + std::to_string(component_index(l, integral[1]));
}

}  // namespace

std::string
format_vrr_cartesian_kernel(const int lb)
{
    const auto target = "s" + shell_label(lb);

    const auto ntarget = cartesian_count(lb);

    // the lower (s|lb-1) and (s|lb-2) Cartesian integrals the single step consumes.

    std::vector<std::pair<int, std::string>> inputs;

    inputs.push_back({lb - 1, "s" + shell_label(lb - 1)});

    if (lb >= 2) inputs.push_back({lb - 2, "s" + shell_label(lb - 2)});

    // the single ket VRR step of every Cartesian (s|lb) component.

    const T2COverlapDriver drv;

    const auto integral = I2CIntegral(I1CPair("GA", 0), I1CPair("GB", lb), Operator("1"), 0, {});

    const auto comps = integral.components<T1CPair, T1CPair>();

    std::ostringstream os;

    os << signature_cartesian_text(lb) << "\n";
    os << "{\n";

    os << "    // bra and ket contracted basis functions\n";
    os << "    const auto bra = pair.bra();\n";
    os << "    const auto ket = pair.ket();\n\n";

    os << "    // primitive exponents (for fe = 1 / (2 eta))\n";
    os << "    const auto bra_exps = bra.get_exponents();\n";
    os << "    const auto ket_exps = ket.get_exponents();\n\n";

    os << "    // number of atom pairs (columns)\n";
    os << "    const auto npairs = " << target << ".ncols();\n\n";

    os << "    // loop over primitive basis-function pairs (the result is primitive)\n";
    os << "    for (std::size_t p = 0; p < bra.number_of_primitive_functions(); p++)\n";
    os << "    {\n";
    os << "        for (std::size_t q = 0; q < ket.number_of_primitive_functions(); q++)\n";
    os << "        {\n";

    os << "            const auto ip = p * ket.number_of_primitive_functions() + q;\n\n";

    os << "            const auto fe = 1.0 / (2.0 * (bra_exps[p] + ket_exps[q]));   // 1 / (2 eta)\n\n";

    // input and output row pointers, offset by ip times the integral's component count.
    for (const auto& [order, label] : inputs)
    {
        const auto ncomp = cartesian_count(order);

        os << "            // lower integral (s|" << shell_label(order) << ")\n";

        for (int r = 0; r < ncomp; r++)
        {
            os << "            auto " << label << "_" << r << " = " << label << ".row(ip * " << ncomp
               << " + " << r << ");\n";
        }
    }

    os << "\n";
    os << "            auto pc_x = pc.row(ip * 3 + 0);\n";
    os << "            auto pc_y = pc.row(ip * 3 + 1);\n";
    os << "            auto pc_z = pc.row(ip * 3 + 2);\n\n";

    os << "            // Cartesian (s|" << shell_label(lb) << ") result\n";

    for (int k = 0; k < ntarget; k++)
    {
        os << "            auto " << target << "_" << k << " = " << target << ".row(ip * " << ntarget
           << " + " << k << ");\n";
    }

    // one SIMD step per Cartesian component.
    for (int k = 0; k < ntarget; k++)
    {
        const auto dist = drv.apply_ket_vrr(R2CTerm(comps[k]));

        // gather the pointers this component touches.
        std::set<std::string> used;

        for (std::size_t t = 0; t < dist.terms(); t++)
        {
            used.insert(lower_row_name(dist[t].integral()));

            for (const auto& fact : dist[t].factors())
            {
                if (fact.name() == "PB") used.insert(pc_pointer(fact));
            }
        }

        std::vector<std::string> aligned(used.begin(), used.end());

        aligned.push_back(target + "_" + std::to_string(k));

        os << "\n";
        os << "            #pragma omp simd aligned(";

        for (std::size_t n = 0; n < aligned.size(); n++) os << (n ? ", " : "") << aligned[n];

        os << " : 64)\n";

        os << "            for (std::size_t i = 0; i < npairs; i++)\n";
        os << "            {\n";

        os << "                " << target << "_" << k << "[i] = ";

        for (std::size_t t = 0; t < dist.terms(); t++)
        {
            const auto& rterm = dist[t];

            std::string pc_factor;

            int fe_power = 0;

            for (const auto& fact : rterm.factors())
            {
                if (fact.name() == "PB") pc_factor = pc_pointer(fact);

                else if (fact.name() == "1/eta") fe_power += rterm.factor_order(fact);
            }

            const auto mag = magnitude(rterm.prefactor());

            const bool neg = rterm.prefactor().numerator() < 0;

            if (t == 0) { if (neg) os << "-"; }
            else os << (neg ? " - " : " + ");

            if (!(mag == Fraction(1))) os << terminating_decimal(mag) << " * ";

            if (!pc_factor.empty()) os << pc_factor << "[i] * ";

            for (int n = 0; n < fe_power; n++) os << "fe * ";

            os << lower_row_name(rterm.integral()) << "[i]";
        }

        os << ";\n";
        os << "            }\n";
    }

    os << "        }\n";
    os << "    }\n";
    os << "}\n";

    return os.str();
}

std::string
format_vrr_cartesian_signature(const int lb)
{
    return signature_cartesian_text(lb);
}

std::string
format_vrr_spherical_signature(const int lb)
{
    return signature_spherical_text(lb);
}
