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

#include "two_center_hrr_emitter.hpp"

#include <algorithm>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <utility>
#include <vector>

#include "operator.hpp"
#include "spherical_harmonics.hpp"
#include "t2c_defs.hpp"
#include "t2c_hrr_driver.hpp"
#include "tensor.hpp"

namespace {  // HRR emitter helpers

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

/// The row-pointer name of a base integral component, e.g. (s|d_xy) -> "sd_1".
/// The flat row index is bra-major: bra_index * cartesian_count(ket) + ket_index.
std::string
base_row_name(const T2CIntegral& integral)
{
    const auto bl = integral[0].order();

    const auto kl = integral[1].order();

    const auto bra = component_index(bl, integral[0]);

    const auto ket = component_index(kl, integral[1]);

    const auto flat = bra * cartesian_count(kl) + ket;

    return shell_label(bl) + shell_label(kl) + "_" + std::to_string(flat);
}

/// Fully reduces a Cartesian target component to the base integrals by repeated
/// single-step horizontal recurrence on the incremented side, until that side
/// reaches angular momentum 0 (bra for la <= lb, ket otherwise). The returned
/// terms each carry the accumulated AB factors and rational prefactor.
std::vector<R2CTerm>
fully_reduce(const T2CHRRDriver& drv, const R2CTerm& start, const bool bra_incremented)
{
    std::vector<R2CTerm> base;

    std::vector<R2CTerm> work{start};

    while (!work.empty())
    {
        std::vector<R2CTerm> next;

        for (const auto& term : work)
        {
            const auto order = bra_incremented ? term[0].order() : term[1].order();

            if (order == 0)
            {
                base.push_back(term);

                continue;
            }

            const auto dist = bra_incremented ? drv.apply_bra_vrr(term) : drv.apply_ket_vrr(term);

            for (std::size_t j = 0; j < dist.terms(); j++) next.push_back(dist[j]);
        }

        work = next;
    }

    return base;
}

/// The sorted AB-distance pointer factors of a recurrence term, with repetition
/// for powers, e.g. AB_x^2 -> {"ab_x", "ab_x"}.
std::vector<std::string>
ab_factors_of(const R2CTerm& term)
{
    std::vector<std::string> factors;

    for (const auto& fact : term.factors())
    {
        if (fact.name() == "AB")
        {
            for (int n = 0; n < term.factor_order(fact); n++) factors.push_back(fact.label());
        }
    }

    std::sort(factors.begin(), factors.end());

    return factors;
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
/// Assumes the rational terminates (see is_terminating).
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

/// One contribution to a spherical target row: a signed coefficient times a base
/// row pointer, times a product of AB-distance pointers.
struct Contribution
{
    Fraction                 coeff;       // combined transform * recurrence coefficient
    int                      radicand;    // square-free radical of the transform coefficient
    std::string              row;         // base integral row pointer (e.g. "sg_4")
    std::vector<std::string> ab_factors;  // AB-distance pointers (e.g. {"ab_x", "ab_x"})
};

/// The coefficient-and-row part of a contribution (no sign, no AB factors), e.g.
/// "0.5 * f3 * sd_0[i]". Terminating rationals are inlined as decimals (1/2 ->
/// 0.5); non-terminating rationals and square roots use the hoisted constants
/// (q1_3, f3) declared once outside the loops.
std::string
magnitude_body(const Contribution& c)
{
    std::string body;

    const auto mag = magnitude(c.coeff);

    if (!(mag == Fraction(1)))
    {
        body += (is_terminating(mag) ? terminating_decimal(mag) : rational_name(mag)) + " * ";
    }

    if (c.radicand != 1) body += "f" + std::to_string(c.radicand) + " * ";

    body += c.row + "[i]";

    return body;
}

/// The AB-distance product as code, e.g. {"ab_x", "ab_x"} -> "ab_x[i] * ab_x[i]".
std::string
ab_product_body(const std::vector<std::string>& ab_factors)
{
    std::string body;

    for (std::size_t i = 0; i < ab_factors.size(); i++)
    {
        body += (i ? " * " : "") + ab_factors[i] + "[i]";
    }

    return body;
}

/// The full magnitude (no sign) of a contribution, AB factors appended.
std::string
contribution_body(const Contribution& c)
{
    std::string body = magnitude_body(c);

    for (const auto& ab : c.ab_factors) body += " * " + ab + "[i]";

    return body;
}

/// The positive greatest common divisor of a group's rational coefficients, as a
/// fraction (gcd of numerators over lcm of denominators). Every coefficient is an
/// integer multiple of it, so dividing it out leaves integer inner coefficients.
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

/// The base integrals the recurrence consumes for a (la|lb) target, one CArray
/// parameter each: (s|lb..lb+la) when the bra is incremented (la <= lb), or
/// (la..la+lb|s) when the ket is incremented.
std::vector<std::pair<int, int>>
hrr_bases(const int la, const int lb)
{
    std::vector<std::pair<int, int>> bases;

    if (la <= lb)
    {
        for (int k = lb; k <= lb + la; k++) bases.push_back({0, k});
    }
    else
    {
        for (int k = la; k <= la + lb; k++) bases.push_back({k, 0});
    }

    return bases;
}

/// The kernel signature "void compute_<la>_<lb>(<inputs>)" (no body, no ';'), with
/// the base-integral CArrays, the AB distances, and the target out-parameter.
std::string
signature_text(const int la, const int lb)
{
    const auto name = "compute_" + shell_label(la) + "_" + shell_label(lb);

    const auto indent = std::string(name.size() + 6, ' ');  // align under "void name("

    std::ostringstream os;

    os << "void " << name << "(";

    for (const auto& [bl, kl] : hrr_bases(la, lb))
    {
        os << "const osfunc::CArray<double>& " << shell_label(bl) << shell_label(kl) << ",\n"
           << indent;
    }

    os << "const osfunc::CArray<double>& ab,\n";
    os << indent << "      osfunc::CArray<double>& " << shell_label(la) << shell_label(lb) << ")";

    return os.str();
}

}  // namespace

std::string
format_hrr_signature(const int la, const int lb)
{
    return signature_text(la, lb);
}

std::string
format_hrr_kernel(const int la, const int lb)
{
    const bool bra_incremented = (la <= lb);

    // the target carries (2*la + 1) * (2*lb + 1) spherical components per block.

    const int target_size = (2 * la + 1) * (2 * lb + 1);

    // the base integrals the recurrence consumes, one CArray parameter each.

    const auto bases = hrr_bases(la, lb);

    const auto target = shell_label(la) + shell_label(lb);

    // the horizontal recurrence of every Cartesian target component, keyed by its
    // (bra, ket) Cartesian components.

    const T2CHRRDriver drv;

    const auto integral = I2CIntegral(I1CPair("GA", la), I1CPair("GB", lb), Operator("1"), 0, {});

    std::map<std::pair<TensorComponent, TensorComponent>, std::vector<R2CTerm>> hrr;

    for (const auto& comp : integral.components<T1CPair, T1CPair>())
    {
        hrr.emplace(std::make_pair(comp[0], comp[1]),
                    fully_reduce(drv, R2CTerm(comp), bra_incremented));
    }

    // for each spherical target component, fold the Cartesian->spherical transform
    // into the recurrence: accumulate the base-row contributions, keyed by (row,
    // AB) so equal terms combine.

    std::vector<std::vector<Contribution>> rows(target_size);

    for (int c = 0; c < target_size; c++)
    {
        const int bra_sph = c / (2 * lb + 1);

        const int ket_sph = c % (2 * lb + 1);

        std::map<std::pair<std::string, std::string>, Contribution> acc;

        for (const auto& term : sphar::two_center_spherical_factors(la, lb, bra_sph, ket_sph))
        {
            for (const auto& rterm : hrr.at(std::make_pair(term.bra, term.ket)))
            {
                const auto row = base_row_name(rterm.integral());

                const auto ab_factors = ab_factors_of(rterm);

                std::string ab_key;

                for (const auto& a : ab_factors) ab_key += a + "*";

                const auto coeff = term.factor.factor * rterm.prefactor();

                const auto key = std::make_pair(row, ab_key);

                if (const auto it = acc.find(key); it != acc.end())
                    it->second.coeff = it->second.coeff + coeff;
                else
                    acc.emplace(key, Contribution{coeff, term.factor.radicand, row, ab_factors});
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

    os << signature_text(la, lb) << "\n";
    os << "{\n";

    os << "    // number of spherical components in the target integral\n";
    os << "    const std::size_t ncomps = " << target_size << ";\n\n";

    os << "    // integral blocks are packed one after another in the rows\n";
    os << "    const auto nblocks = " << target << ".nrows() / ncomps;\n\n";

    os << "    // number of atom pairs (columns)\n";
    os << "    const auto npairs = " << target << ".ncols();\n\n";

    // hoist the transformation factors that do not reduce to a clean decimal: the
    // square roots and the non-terminating rationals.

    std::set<int> radicals;

    std::map<std::pair<long long, long long>, Fraction> rationals;

    for (const auto& row : rows)
    {
        for (const auto& contrib : row)
        {
            if (contrib.radicand != 1) radicals.insert(contrib.radicand);

            const auto mag = magnitude(contrib.coeff);

            if (!(mag == Fraction(1)) && !is_terminating(mag))
            {
                rationals[{mag.numerator(), mag.denominator()}] = mag;
            }
        }
    }

    if (!radicals.empty() || !rationals.empty())
    {
        os << "    // transformation factors\n";

        for (const auto r : radicals)
        {
            os << "    const double f" << r << " = std::sqrt(" << r << ".0);\n";
        }

        for (const auto& [key, value] : rationals)
        {
            os << "    const double " << rational_name(value) << " = (" << value.numerator()
               << ".0 / " << value.denominator() << ".0);\n";
        }

        os << "\n";
    }

    os << "    // AB distances (per atom-pair column)\n";
    os << "    auto ab_x = ab.row(0);\n";
    os << "    auto ab_y = ab.row(1);\n";
    os << "    auto ab_z = ab.row(2);\n\n";

    os << "    // outermost loop runs over the integral blocks\n";
    os << "    for (std::size_t iblock = 0; iblock < nblocks; iblock++)\n";
    os << "    {\n";

    for (const auto& [bl, kl] : bases)
    {
        const auto label = shell_label(bl) + shell_label(kl);

        const auto nrows = cartesian_count(bl) * cartesian_count(kl);

        os << "        // base integral (" << shell_label(bl) << "|" << shell_label(kl) << "): "
           << nrows << " Cartesian components\n";
        os << "        const auto " << label << "_off = iblock * " << nrows << ";\n";

        for (int r = 0; r < nrows; r++)
        {
            os << "        auto " << label << "_" << r << " = " << label << ".row(" << label
               << "_off + " << r << ");\n";
        }

        os << "\n";
    }

    os << "        // target (" << shell_label(la) << "|" << shell_label(lb) << "): " << target_size
       << " spherical components\n";
    os << "        const auto " << target << "_off = iblock * " << target_size << ";\n";

    for (int r = 0; r < target_size; r++)
    {
        os << "        auto " << target << "_" << r << " = " << target << ".row(" << target
           << "_off + " << r << ");\n";
    }

    // the SIMD contraction over atom-pair columns, one loop per spherical component
    // of the incremented side (bra when la <= lb, ket otherwise), each sweeping the
    // other side's components. Grouping on the incremented side pins one AB axis.

    const int nbra = 2 * la + 1;

    const int nket = 2 * lb + 1;

    const int ngroups = bra_incremented ? nbra : nket;

    for (int g = 0; g < ngroups; g++)
    {
        // the target component indices in this group (bra-major flat index); the
        // incremented side indexes the group, the other side runs inside it.

        std::vector<int> comps;

        if (bra_incremented)
        {
            for (int k = 0; k < nket; k++) comps.push_back(g * nket + k);
        }
        else
        {
            for (int b = 0; b < nbra; b++) comps.push_back(b * nket + g);
        }

        // one SIMD loop per target component.

        for (const int c : comps)
        {
            std::set<std::string> used;

            for (const auto& contrib : rows[c])
            {
                used.insert(contrib.row);

                for (const auto& ab : contrib.ab_factors) used.insert(ab);
            }

            std::vector<std::string> aligned(used.begin(), used.end());

            aligned.push_back(target + "_" + std::to_string(c));

            os << "\n";
            os << "        // " << (bra_incremented ? "bra" : "ket") << " spherical component " << g
               << ", target row " << c << "\n";
            os << "        #pragma omp simd aligned(";

            for (std::size_t n = 0; n < aligned.size(); n++) os << (n ? ", " : "") << aligned[n];

            os << " : 64)\n";

            os << "        for (std::size_t i = 0; i < npairs; i++)\n";
            os << "        {\n";

            const auto lead = std::string("            ") + target + "_" + std::to_string(c) + "[i] = ";

            const auto hang = std::string(lead.size() - 2, ' ');

            // group the contributions by their AB-distance product, ordered by
            // descending degree then lexicographically (so the bare, AB-free terms
            // come last).
            std::map<std::vector<std::string>, std::vector<Contribution>> by_ab;

            for (const auto& contrib : rows[c]) by_ab[contrib.ab_factors].push_back(contrib);

            std::vector<std::pair<std::vector<std::string>, std::vector<Contribution>>> groups(
                by_ab.begin(), by_ab.end());

            std::sort(groups.begin(), groups.end(), [](const auto& a, const auto& b) {
                if (a.first.size() != b.first.size()) return a.first.size() > b.first.size();

                return a.first < b.first;
            });

            // a summand is one line of the sum: a factored AB group "(...) * ab"
            // (>= 2 terms sharing the same AB product) or a single bare term.
            struct Summand
            {
                bool        neg;
                std::string body;
            };

            std::vector<Summand> summands;

            for (const auto& [ab, terms] : groups)
            {
                const auto g = group_gcd(terms);

                const int  radicand = terms.front().radicand;

                const bool first_neg = terms.front().coeff.numerator() < 0;

                // a group is worth factoring when it has >= 2 terms and shares a
                // common coefficient, radical, or AB product to pull out.
                const bool factorable =
                    terms.size() > 1 && (!(g == Fraction(1)) || radicand != 1 || !ab.empty());

                if (factorable)
                {
                    // divide out the signed common factor (sign from the first term
                    // so the leading inner term is positive); inner coeffs are integers.
                    const long long gnum = first_neg ? -g.numerator() : g.numerator();

                    const long long gden = g.denominator();

                    std::string inner;

                    for (std::size_t t = 0; t < terms.size(); t++)
                    {
                        const auto& ct = terms[t].coeff;

                        const Fraction ic(static_cast<int>(ct.numerator() * gden),
                                          static_cast<int>(ct.denominator() * gnum));

                        const bool ineg = ic.numerator() < 0;

                        inner += (t == 0) ? (ineg ? "-" : "") : (ineg ? " - " : " + ");

                        const auto imag = magnitude(ic);

                        if (!(imag == Fraction(1))) inner += terminating_decimal(imag) + " * ";

                        inner += terms[t].row + "[i]";
                    }

                    std::string outer;

                    if (!(g == Fraction(1)))
                    {
                        outer += (is_terminating(g) ? terminating_decimal(g) : rational_name(g)) + " * ";
                    }

                    if (radicand != 1) outer += "f" + std::to_string(radicand) + " * ";

                    std::string body = outer + "(" + inner + ")";

                    if (!ab.empty()) body += " * " + ab_product_body(ab);

                    summands.push_back({first_neg, body});
                }
                else
                {
                    for (const auto& contrib : terms)
                    {
                        summands.push_back({contrib.coeff.numerator() < 0, contribution_body(contrib)});
                    }
                }
            }

            os << lead;

            if (summands.empty()) os << "0.0";

            // one summand per line, the operator leading each continuation line.
            for (std::size_t s = 0; s < summands.size(); s++)
            {
                if (s == 0)
                {
                    if (summands[s].neg) os << "-";
                }
                else
                {
                    os << "\n" << hang << (summands[s].neg ? "- " : "+ ");
                }

                os << summands[s].body;
            }

            os << ";\n";

            os << "        }\n";
        }
    }

    os << "    }\n";
    os << "}\n";

    return os.str();
}
