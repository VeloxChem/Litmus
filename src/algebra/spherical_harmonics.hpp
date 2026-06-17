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

#ifndef spherical_harmonics_hpp
#define spherical_harmonics_hpp

#include <string>
#include <utility>
#include <vector>

#include "fraction.hpp"
#include "tensor_component.hpp"

namespace sphar {  // real solid harmonics namespace

/// An exact Cartesian-to-spherical transformation coefficient, stored as a
/// rational factor times the square root of a square-free integer radicand:
///
///     value = factor * sqrt(radicand)
///
/// This is the closed form every real solid-harmonic coefficient takes in the
/// Helgaker-Jorgensen-Olsen convention (a radicand of 1 is a plain rational).
struct SphericalFactor
{
    /// The rational factor in front of the square root.
    Fraction factor;

    /// The square-free radicand under the square root (>= 1).
    int radicand;

    /// Creates a zero coefficient (0 * sqrt(1)).
    SphericalFactor();

    /// Creates a coefficient factor * sqrt(radicand).
    /// @param factor The rational factor in front of the square root.
    /// @param radicand The square-free radicand under the square root (>= 1).
    SphericalFactor(const Fraction& factor, const int radicand);

    /// Compares this coefficient with another for equality.
    /// @param other The other coefficient to compare.
    /// @return True if the coefficients are equal, false otherwise.
    bool operator==(const SphericalFactor& other) const;

    /// Compares this coefficient with another for inequality.
    /// @param other The other coefficient to compare.
    /// @return True if the coefficients differ, false otherwise.
    bool operator!=(const SphericalFactor& other) const;

    /// @return True if the coefficient is exactly zero.
    bool is_zero() const;

    /// @return The numerical value factor * sqrt(radicand).
    double value() const;

    /// @return A textual representation, e.g. "1/2 * sqrt(3)" or "-1/2".
    std::string to_string() const;
};

/// Multiplies two transformation coefficients, reducing the combined radical back
/// to square-free form, e.g. sqrt(3) * sqrt(3) -> 3 and sqrt(2) * sqrt(6) ->
/// 2 * sqrt(3). This is the coefficient algebra a two-center transform needs, as
/// the bra and ket factors multiply per Cartesian term.
/// @param lhs The left coefficient.
/// @param rhs The right coefficient.
/// @return The product factor_l * factor_r * sqrt(radicand_l * radicand_r).
SphericalFactor operator*(const SphericalFactor& lhs, const SphericalFactor& rhs);

/// The non-zero Cartesian terms of a real solid harmonic: each entry pairs a
/// Cartesian component x^i y^j z^k with its transformation coefficient, in
/// canonical tensor order.
using VSphericalTerms = std::vector<std::pair<TensorComponent, SphericalFactor>>;

/// Computes the Cartesian expansion of the real solid harmonic S_{l,m} in the
/// Helgaker-Jorgensen-Olsen convention.
/// @param l The angular momentum (l >= 0).
/// @param m The order, with -l <= m <= l.
/// @return The non-zero (Cartesian component, coefficient) pairs in canonical
///         tensor order, or an empty vector when (l, m) is out of range.
VSphericalTerms spherical_factors(const int l, const int m);

/// Computes the Cartesian expansion of the real solid harmonic for a spherical
/// component index. The components run m = -l, ..., +l, so the index is m + l
/// (the ordering used by VeloxChem's spherical-momentum transformation).
/// @param l The angular momentum (l >= 0).
/// @param component The spherical component index, 0 <= component <= 2 * l.
/// @return The non-zero (Cartesian component, coefficient) pairs in canonical
///         tensor order, or an empty vector when the index is out of range.
VSphericalTerms spherical_component_factors(const int l, const int component);

/// One term of a two-center Cartesian-to-spherical transform: the bra (center A)
/// and ket (center B) Cartesian components and the combined exact coefficient.
struct SphericalPairTerm
{
    /// The bra (center A) Cartesian component.
    TensorComponent bra;

    /// The ket (center B) Cartesian component.
    TensorComponent ket;

    /// The combined transformation coefficient (bra coefficient * ket coefficient).
    SphericalFactor factor;

    /// Compares this pair term with another for equality.
    /// @param other The other pair term to compare.
    /// @return True if the components and coefficient are equal, false otherwise.
    bool operator==(const SphericalPairTerm& other) const;
};

/// The non-zero Cartesian terms of a two-center spherical component pair, in
/// canonical (bra, ket) tensor order.
using VSphericalPairTerms = std::vector<SphericalPairTerm>;

/// Computes the exact Cartesian-to-spherical transform of one (la|lb) two-center
/// block component: a spherical bra/ket component pair expanded into its non-zero
/// Cartesian (bra, ket) terms, each carrying the combined coefficient
/// bra_coefficient * ket_coefficient. The transform is the tensor product of the
/// per-shell real solid harmonics, so it is identical for every two-center
/// integral operator (overlap, kinetic energy, nuclear potential, ...). The
/// terms run in canonical order: ket fastest within each bra component.
/// @param la The bra (center A) angular momentum.
/// @param lb The ket (center B) angular momentum.
/// @param bra_component The bra spherical component index, 0 <= . <= 2 * la.
/// @param ket_component The ket spherical component index, 0 <= . <= 2 * lb.
/// @return The non-zero (bra, ket, coefficient) terms, or an empty vector when
///         either component index is out of range.
VSphericalPairTerms two_center_spherical_factors(const int la,
                                                 const int lb,
                                                 const int bra_component,
                                                 const int ket_component);

}  // namespace sphar

#endif /* spherical_harmonics_hpp */
