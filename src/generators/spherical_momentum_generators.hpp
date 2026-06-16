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

#ifndef spherical_momentum_generators_hpp
#define spherical_momentum_generators_hpp

#include <string>

/// Generates VeloxChem's SphericalMomentum.hpp: the spher_mom::transformation_factors<N>
/// template that returns, for each spherical component, the Cartesian-to-spherical
/// transformation factors of the real solid harmonics. The coefficients come from
/// the exact sphar:: solid-harmonic expansion, so the emitted header no longer has
/// to be hand-maintained.
class SphericalMomentumGenerator
{
public:
    /// Creates a spherical-momentum header generator.
    SphericalMomentumGenerator() = default;

    /// Writes SphericalMomentum.hpp for angular momenta 0..max_ang_mom into the
    /// current working directory.
    /// @param max_ang_mom The maximum angular momentum to tabulate (>= 0).
    void generate(const int max_ang_mom) const;
};

/// Builds the text of the SphericalMomentum.hpp header for angular momenta
/// 0..max_ang_mom. Exposed (apart from generate) so the emitted source can be
/// unit-tested without touching the file system.
/// @param max_ang_mom The maximum angular momentum to tabulate (>= 0).
/// @return The full contents of the generated header.
std::string format_spherical_momentum(const int max_ang_mom);

#endif /* spherical_momentum_generators_hpp */
