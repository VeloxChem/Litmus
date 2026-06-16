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

#include "spherical_harmonics.hpp"

#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

#include "tensor.hpp"

namespace sphar {  // real solid harmonics namespace

namespace {  // exact-arithmetic helpers

/// Splits a non-negative integer n into n = g * g * s, where s is square-free.
/// @param n The non-negative integer to factor.
/// @return The pair (g, s) with g the extracted root and s the square-free part.
std::pair<long long, long long>
reduce_radical(long long n)
{
    if (n == 0) return {0, 1};

    long long g = 1, s = 1;

    for (long long p = 2; p * p <= n; p++)
    {
        int e = 0;

        while (n % p == 0)
        {
            n /= p;
            e++;
        }

        for (int i = 0; i < e / 2; i++) g *= p;

        if (e % 2 == 1) s *= p;
    }

    s *= n;  // the remaining factor is 1 or a single prime (odd exponent)

    return {g, s};
}

/// An exact sum of rational multiples of square roots of square-free integers,
/// i.e. a value of the form sum_k c_k * sqrt(r_k). This is closed under the
/// addition and multiplication the solid-harmonic recurrences require; each map
/// entry maps a square-free radicand (>= 1) to its non-zero rational coefficient.
class RadicalSum
{
    /// The non-zero terms, keyed by square-free radicand.
    std::map<int, Fraction> _terms;

    /// Adds c * sqrt(r) into the sum, dropping the entry if it cancels to zero.
    void _add_term(const int radicand, const Fraction& coeff)
    {
        if (coeff == Fraction(0)) return;

        if (const auto it = _terms.find(radicand); it != _terms.end())
        {
            const auto sum = it->second + coeff;

            if (sum == Fraction(0))
            {
                _terms.erase(it);
            }
            else
            {
                it->second = sum;
            }
        }
        else
        {
            _terms[radicand] = coeff;
        }
    }

public:
    RadicalSum() = default;

    /// @return The rational value f as the RadicalSum f * sqrt(1).
    static RadicalSum from_fraction(const Fraction& f)
    {
        RadicalSum sum;

        sum._add_term(1, f);

        return sum;
    }

    /// @return The integer n as the RadicalSum n * sqrt(1).
    static RadicalSum from_int(const int n) { return from_fraction(Fraction(n)); }

    /// Builds sqrt(r) for a non-negative rational r = p / q, normalized to
    /// (g / q) * sqrt(s) where p * q = g * g * s and s is square-free.
    /// @param r The non-negative rational under the square root.
    /// @return The RadicalSum representing sqrt(r).
    static RadicalSum from_sqrt(const Fraction& r)
    {
        RadicalSum sum;

        if (r.numerator() == 0) return sum;

        const auto [g, s] = reduce_radical(static_cast<long long>(r.numerator()) * r.denominator());

        sum._add_term(static_cast<int>(s), Fraction(static_cast<int>(g), r.denominator()));

        return sum;
    }

    /// @return This sum with every term scaled by the rational f.
    RadicalSum scaled(const Fraction& f) const
    {
        RadicalSum out;

        for (const auto& [radicand, coeff] : _terms) out._add_term(radicand, coeff * f);

        return out;
    }

    /// @return The sum of this and another RadicalSum.
    RadicalSum operator+(const RadicalSum& other) const
    {
        RadicalSum out = *this;

        for (const auto& [radicand, coeff] : other._terms) out._add_term(radicand, coeff);

        return out;
    }

    /// @return The difference of this and another RadicalSum.
    RadicalSum operator-(const RadicalSum& other) const
    {
        return *this + other.scaled(Fraction(-1));
    }

    /// @return The product of this and another RadicalSum (sqrt(a) sqrt(b) =
    ///         sqrt(a b), then re-normalized to square-free radicands).
    RadicalSum operator*(const RadicalSum& other) const
    {
        RadicalSum out;

        for (const auto& [r1, c1] : _terms)
        {
            for (const auto& [r2, c2] : other._terms)
            {
                const auto [g, s] = reduce_radical(static_cast<long long>(r1) * r2);

                out._add_term(static_cast<int>(s), c1 * c2 * Fraction(static_cast<int>(g)));
            }
        }

        return out;
    }

    /// @return True if the sum has no terms (is exactly zero).
    bool is_zero() const { return _terms.empty(); }

    /// Collapses the sum to a single coefficient factor * sqrt(radicand). Every
    /// real solid-harmonic coefficient is a single such term; a leftover sum of
    /// distinct radicals would mean a derivation bug, so it is rejected.
    /// @return The single-radical coefficient.
    SphericalFactor to_factor() const
    {
        if (_terms.empty()) return SphericalFactor();

        if (_terms.size() > 1)
        {
            throw std::logic_error("sphar: coefficient is not a single square root");
        }

        const auto& [radicand, coeff] = *_terms.begin();

        return SphericalFactor(coeff, radicand);
    }
};

/// A symbolic polynomial in x, y, z with RadicalSum coefficients, keyed by the
/// Cartesian monomial x^i y^j z^k it multiplies.
using Polynomial = std::map<TensorComponent, RadicalSum>;

/// Adds c * (monomial) into a polynomial, pruning the entry if it cancels.
void
add_monomial(Polynomial& poly, const TensorComponent& monomial, const RadicalSum& coeff)
{
    if (coeff.is_zero()) return;

    if (const auto it = poly.find(monomial); it != poly.end())
    {
        auto sum = it->second + coeff;

        if (sum.is_zero())
        {
            poly.erase(it);
        }
        else
        {
            it->second = std::move(sum);
        }
    }
    else
    {
        poly[monomial] = coeff;
    }
}

/// @return The sum of two polynomials.
Polynomial
add(const Polynomial& lhs, const Polynomial& rhs)
{
    Polynomial out = lhs;

    for (const auto& [monomial, coeff] : rhs) add_monomial(out, monomial, coeff);

    return out;
}

/// @return The difference of two polynomials.
Polynomial
subtract(const Polynomial& lhs, const Polynomial& rhs)
{
    Polynomial out = lhs;

    for (const auto& [monomial, coeff] : rhs) add_monomial(out, monomial, coeff.scaled(Fraction(-1)));

    return out;
}

/// @return The polynomial with every coefficient multiplied by a RadicalSum.
Polynomial
scale(const Polynomial& poly, const RadicalSum& coeff)
{
    Polynomial out;

    for (const auto& [monomial, mcoeff] : poly) add_monomial(out, monomial, mcoeff * coeff);

    return out;
}

/// @return The polynomial multiplied by the Cartesian variable on the given axis
///         (raises that axis' exponent by one).
Polynomial
multiply_by_axis(const Polynomial& poly, const char axis)
{
    Polynomial out;

    for (const auto& [monomial, coeff] : poly)
    {
        add_monomial(out, monomial.shift(axis, 1).value(), coeff);
    }

    return out;
}

/// @return The polynomial multiplied by r^2 = x^2 + y^2 + z^2.
Polynomial
multiply_by_rsq(const Polynomial& poly)
{
    Polynomial out;

    for (const auto& [monomial, coeff] : poly)
    {
        add_monomial(out, TensorComponent(monomial['x'] + 2, monomial['y'], monomial['z']), coeff);
        add_monomial(out, TensorComponent(monomial['x'], monomial['y'] + 2, monomial['z']), coeff);
        add_monomial(out, TensorComponent(monomial['x'], monomial['y'], monomial['z'] + 2), coeff);
    }

    return out;
}

/// Computes the Cartesian polynomial of the real solid harmonic S_{l,m},
/// following the recurrences of SolidHarmonics.py (Helgaker-Jorgensen-Olsen
/// convention). Results are memoized across the recursion.
/// @param l The angular momentum.
/// @param m The order.
/// @param memo The memoization table keyed by (l, m).
/// @return The Cartesian polynomial of S_{l,m} (empty when (l, m) is invalid).
const Polynomial&
solid_harmonic(const int l, const int m, std::map<std::pair<int, int>, Polynomial>& memo)
{
    static const Polynomial zero;

    if (l < 0 || m < -l || m > l) return zero;

    const auto key = std::make_pair(l, m);

    if (const auto it = memo.find(key); it != memo.end()) return it->second;

    Polynomial out;

    if (l == 0)
    {
        out[TensorComponent(0, 0, 0)] = RadicalSum::from_int(1);
    }
    else if (l == m || l == -m)
    {
        // diagonal recurrence: grow the top/bottom harmonic from level l - 1.

        const int  delta = (l == 1) ? 1 : 0;
        const auto inside = (delta == 1) ? Fraction(2 * l - 1, 2 * l) * Fraction(2)
                                         : Fraction(2 * l - 1, 2 * l);
        const auto pf = RadicalSum::from_sqrt(inside);

        if (l == m)
        {
            out = scale(multiply_by_axis(solid_harmonic(l - 1, l - 1, memo), 'x'), pf);

            if (delta == 0)
            {
                const auto term = scale(multiply_by_axis(solid_harmonic(l - 1, -(l - 1), memo), 'y'), pf);

                out = subtract(out, term);
            }
        }
        else
        {
            out = scale(multiply_by_axis(solid_harmonic(l - 1, l - 1, memo), 'y'), pf);

            if (delta == 0)
            {
                const auto term = scale(multiply_by_axis(solid_harmonic(l - 1, -(l - 1), memo), 'x'), pf);

                out = add(out, term);
            }
        }
    }
    else
    {
        // vertical recurrence: step down z from level l - 1 and l - 2.

        const auto pf = RadicalSum::from_sqrt(Fraction(1, (l + m) * (l - m)));

        const auto term1 = scale(multiply_by_axis(solid_harmonic(l - 1, m, memo), 'z'),
                                 pf.scaled(Fraction(2 * l - 1)));

        const auto coeff2 = pf * RadicalSum::from_sqrt(Fraction((l - 1 + m) * (l - 1 - m)));

        const auto term2 = scale(multiply_by_rsq(solid_harmonic(l - 2, m, memo)), coeff2);

        out = subtract(term1, term2);
    }

    return memo.emplace(key, std::move(out)).first->second;
}

}  // namespace

SphericalFactor::SphericalFactor()

    : factor(0)

    , radicand(1)
{
}

SphericalFactor::SphericalFactor(const Fraction& factor, const int radicand)

    : factor(factor)

    , radicand(radicand)
{
}

bool
SphericalFactor::operator==(const SphericalFactor& other) const
{
    return (factor == other.factor) && (radicand == other.radicand);
}

bool
SphericalFactor::operator!=(const SphericalFactor& other) const
{
    return !((*this) == other);
}

bool
SphericalFactor::is_zero() const
{
    return factor == Fraction(0);
}

double
SphericalFactor::value() const
{
    return (static_cast<double>(factor.numerator()) / factor.denominator()) *
           std::sqrt(static_cast<double>(radicand));
}

std::string
SphericalFactor::to_string() const
{
    if (radicand == 1) return factor.to_string();

    return factor.to_string() + " * sqrt(" + std::to_string(radicand) + ")";
}

VSphericalTerms
spherical_factors(const int l, const int m)
{
    if (l < 0 || m < -l || m > l) return {};

    std::map<std::pair<int, int>, Polynomial> memo;

    const auto& poly = solid_harmonic(l, m, memo);

    // emit the non-zero terms in canonical tensor order (the order Litmus and
    // VeloxChem index Cartesian components by).

    VSphericalTerms terms;

    for (const auto& component : Tensor(l).components())
    {
        if (const auto it = poly.find(component); it != poly.end())
        {
            terms.emplace_back(component, it->second.to_factor());
        }
    }

    return terms;
}

VSphericalTerms
spherical_component_factors(const int l, const int component)
{
    if (l < 0 || component < 0 || component > 2 * l) return {};

    return spherical_factors(l, component - l);
}

}  // namespace sphar
