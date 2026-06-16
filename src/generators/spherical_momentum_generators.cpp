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

#include "spherical_momentum_generators.hpp"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <vector>

#include "spherical_harmonics.hpp"
#include "tensor.hpp"

namespace {  // spherical-momentum emitter helpers

/// The VeloxChem source-license header reproduced verbatim, so the regenerated
/// header is a drop-in replacement.
const char* const license_header =
    R"(//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
)";

/// The spectroscopic shell label of an angular momentum (S, P, D, F, ...).
/// @param l The angular momentum.
/// @return The upper-case shell label.
std::string
shell_label(const int l)
{
    static const std::string labels = "SPDFGHIKLMNOQRTUV";

    return (l >= 0 && l < static_cast<int>(labels.size())) ? std::string(1, labels[l])
                                                           : ("L" + std::to_string(l));
}

/// Formats a rational as an exact terminating decimal literal, e.g. 3/8 -> "0.375".
/// Falls back to a parenthesized ratio when the rational does not terminate (which
/// does not happen for solid-harmonic coefficients, but keeps the output exact).
/// @param value The rational to format.
/// @return The decimal literal.
std::string
format_decimal(const Fraction& value)
{
    long long num = value.numerator();

    const long long den = value.denominator();

    if (den == 1) return std::to_string(num) + ".0";

    const bool negative = num < 0;

    if (negative) num = -num;

    // a rational terminates as a decimal iff its denominator is 2^a * 5^b.

    long long residue = den;

    int twos = 0, fives = 0;

    while (residue % 2 == 0)
    {
        residue /= 2;
        twos++;
    }

    while (residue % 5 == 0)
    {
        residue /= 5;
        fives++;
    }

    if (residue != 1)
    {
        return std::string(negative ? "-" : "") + "(" + std::to_string(num) + ".0 / " +
               std::to_string(den) + ".0)";
    }

    const int places = std::max(twos, fives);

    long long pow10 = 1;

    for (int i = 0; i < places; i++) pow10 *= 10;

    auto digits = std::to_string(num * (pow10 / den));  // exact: den divides 10^places

    if (static_cast<int>(digits.size()) <= places)
    {
        digits = std::string(places + 1 - digits.size(), '0') + digits;
    }

    auto fraction = digits.substr(digits.size() - places);

    while (fraction.size() > 1 && fraction.back() == '0') fraction.pop_back();

    return std::string(negative ? "-" : "") + digits.substr(0, digits.size() - places) + "." + fraction;
}

/// Formats one transformation factor as a "{index, coefficient}" entry, using a
/// shared "f<radicand>" constant for the square-root part.
/// @param index The Cartesian component index.
/// @param factor The exact transformation coefficient.
/// @return The brace-enclosed entry.
std::string
format_entry(const int index, const sphar::SphericalFactor& factor)
{
    std::string coefficient;

    if (factor.radicand == 1)
    {
        coefficient = format_decimal(factor.factor);
    }
    else
    {
        const auto root = "f" + std::to_string(factor.radicand);

        if (factor.factor == Fraction(1))
        {
            coefficient = root;
        }
        else if (factor.factor == Fraction(-1))
        {
            coefficient = "-" + root;
        }
        else
        {
            coefficient = format_decimal(factor.factor) + " * " + root;
        }
    }

    return "{" + std::to_string(index) + ", " + coefficient + "},";
}

/// Maps each Cartesian component of order l to its canonical index.
/// @param l The angular momentum.
/// @return The component-to-index map.
std::map<TensorComponent, int>
cartesian_indices(const int l)
{
    std::map<TensorComponent, int> indices;

    int index = 0;

    for (const auto& component : Tensor(l).components()) indices[component] = index++;

    return indices;
}

/// Emits the "if constexpr (N == l)" block tabulating the transformation factors
/// of all 2l+1 spherical components of angular momentum l.
/// @param os The output stream.
/// @param l The angular momentum.
void
write_block(std::ostringstream& os, const int l)
{
    const auto indices = cartesian_indices(l);

    // gather the spherical expansions and the distinct square-root radicands they
    // use, in order of first appearance (each component shares a single radicand).

    std::vector<sphar::VSphericalTerms> components;

    std::vector<int> radicands;

    for (int component = 0; component <= 2 * l; component++)
    {
        auto terms = sphar::spherical_component_factors(l, component);

        for (const auto& [tensor, factor] : terms)
        {
            if (factor.radicand != 1 &&
                std::find(radicands.begin(), radicands.end(), factor.radicand) == radicands.end())
            {
                radicands.push_back(factor.radicand);
            }
        }

        components.push_back(std::move(terms));
    }

    os << "    // " << shell_label(l) << " type real solid harmonics\n";
    os << "    if constexpr (N == " << l << ")\n";
    os << "    {\n";

    for (const auto radicand : radicands)
    {
        os << "        const double f" << radicand << " = std::sqrt(" << radicand << ".0);\n\n";
    }

    for (int component = 0; component <= 2 * l; component++)
    {
        os << "        if (component == " << component << ")\n";
        os << "            return {\n";

        for (const auto& [tensor, factor] : components[component])
        {
            os << "                " << format_entry(indices.at(tensor), factor) << "\n";
        }

        os << "            };\n";
    }

    os << "\n        return std::vector<std::pair<int, double>>();\n";
    os << "    }\n\n";
}

}  // namespace

std::string
format_spherical_momentum(const int max_ang_mom)
{
    std::ostringstream os;

    os << license_header << "\n";

    os << "#ifndef SphericalMomentum_hpp\n";
    os << "#define SphericalMomentum_hpp\n\n";

    os << "#include <cmath>\n";
    os << "#include <utility>\n";
    os << "#include <vector>\n\n";

    os << "// NOTE: This file is generated by Litmus. Do not edit by hand.\n\n";

    os << "namespace spher_mom {\n\n";

    os << "/// @brief Creates vector of Cartesian to spherical transformation factors.\n";
    os << "/// @tparam N The order of angular momentum tensor.\n";
    os << "/// @param component The component of spherical momentum to generate transformation factors for.\n";
    os << "template <int N>\n";
    os << "auto\n";
    os << "transformation_factors(const int component) -> std::vector<std::pair<int, double>>\n";
    os << "{\n";
    os << "    // Cartesian solid harmonics expansion factors are generated recursively\n";
    os << "    // using Eqs. A2-A6 (see Supporting Info).\n";
    os << "    // J. Chem. Theory Comput. 2020, https://doi.org/10.1021/acs.jctc.9b01296\n\n";

    for (int l = 0; l <= max_ang_mom; l++) write_block(os, l);

    os << "    // TODO: Add higher order transformation factors l > " << max_ang_mom << "\n\n";

    os << "    return std::vector<std::pair<int, double>>();\n";
    os << "}\n\n";

    os << "}  // namespace spher_mom\n\n";

    os << "#endif /* SphericalMomentum_hpp */\n";

    return os.str();
}

void
SphericalMomentumGenerator::generate(const int max_ang_mom) const
{
    std::ofstream fstream;

    fstream.open("SphericalMomentum.hpp", std::ios_base::trunc);

    fstream << format_spherical_momentum(max_ang_mom);

    fstream.close();
}
