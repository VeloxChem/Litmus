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

#ifndef two_center_vrr_emitter_hpp
#define two_center_vrr_emitter_hpp

#include <string>

/// Builds the source of an os2c::ovl spherical two-center VRR kernel: it builds the
/// (s|lb) overlap from the (s|s) seed via the Obara-Saika ket vertical recurrence,
/// fully reduced, and folds in the Cartesian-to-spherical transform so the result
/// is the spherical (s|lb) block. The kernel contracts over the primitive basis
/// functions of the pair.
/// @param lb The ket angular momentum.
/// @return The generated kernel source.
std::string format_vrr_spherical_kernel(const int lb);

/// Builds the source of an os2c::vrr::ovl single-step Cartesian two-center VRR
/// kernel: it builds the primitive Cartesian (s|lb) overlap from the lower
/// (s|lb-1) and (s|lb-2) Cartesian integrals via one Obara-Saika ket step. The
/// kernel is per-primitive (the fe = 1/(2 eta) factor is retained), so the result
/// is contracted downstream.
/// @param lb The ket angular momentum.
/// @return The generated kernel source.
std::string format_vrr_cartesian_kernel(const int lb);

/// Builds the spherical VRR kernel signature "void compute_<lb>_sph(...)" (no
/// body, no terminator), for the declaration in the matching header.
/// @param lb The ket angular momentum.
/// @return The signature text.
std::string format_vrr_spherical_signature(const int lb);

/// Builds the Cartesian VRR kernel signature "void compute_<lb>(...)" (no body, no
/// terminator), for the declaration in the matching header.
/// @param lb The ket angular momentum.
/// @return The signature text.
std::string format_vrr_cartesian_signature(const int lb);

#endif /* two_center_vrr_emitter_hpp */
