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

#ifndef two_center_hrr_emitter_hpp
#define two_center_hrr_emitter_hpp

#include <string>

/// Builds the source of an os2c::hrr Cartesian horizontal-recurrence kernel for a
/// two-center target (la|lb). The kernel takes one contracted Cartesian CArray per
/// base integral the recurrence consumes, the AB distances CArray, and writes the
/// Cartesian target into an out-parameter. The angular momentum is incremented on
/// the bra side when la <= lb, on the ket side otherwise.
/// @param la The bra angular momentum.
/// @param lb The ket angular momentum.
/// @return The generated kernel source.
std::string format_hrr_kernel(const int la, const int lb);

/// Builds the kernel signature "void compute_<la>_<lb>(<inputs>)" (no body, no
/// terminator), for the declaration in the matching header.
/// @param la The bra angular momentum.
/// @param lb The ket angular momentum.
/// @return The signature text.
std::string format_hrr_signature(const int la, const int lb);

#endif /* two_center_hrr_emitter_hpp */
