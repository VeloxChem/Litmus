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

}  // namespace sphar

#endif /* spherical_harmonics_hpp */
