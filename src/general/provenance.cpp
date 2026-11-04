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

#include "provenance.hpp"

#include <sstream>
#include <vector>

#include "litmus_version.hpp"

namespace ost {  // ost namespace

namespace {

/// The VeloxChem license header, carried verbatim by every generated file
/// because the generated code is distributed as part of VeloxChem.
const char* const veloxchem_license_header =
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

/// The meaningful lines of the recorded configuration: blank lines and
/// full-line comments are dropped, everything else is kept verbatim.
std::vector<std::string>&
config_lines()
{
    static std::vector<std::string> lines;

    return lines;
}

/// Removes leading and trailing whitespace.
std::string
trim(const std::string& str)
{
    const auto first = str.find_first_not_of(" \t\r");

    if (first == std::string::npos) return "";

    const auto last = str.find_last_not_of(" \t\r");

    return str.substr(first, last - first + 1);
}

}  // namespace

void
set_generation_config(const std::string& config_text)
{
    auto& lines = config_lines();

    lines.clear();

    std::istringstream stream(config_text);

    std::string line;

    while (std::getline(stream, line))
    {
        line = trim(line);

        if (line.empty() || (line[0] == '#')) continue;

        lines.push_back(line);
    }
}

std::string
litmus_revision()
{
    return LITMUS_GIT_REVISION;
}

std::string
generated_file_banner()
{
    std::ostringstream os;

    os << veloxchem_license_header
       << "\n"
       << "// This file was generated automatically by Litmus, the VeloxChem\n"
       << "// integrals generator (https://github.com/VeloxChem/Litmus).\n"
       << "// DO NOT EDIT: manual changes are lost when the file is regenerated.\n"
       << "//\n"
       << "// Litmus revision: " << litmus_revision() << "\n";

    const auto& lines = config_lines();

    if (!lines.empty())
    {
        os << "//\n"
           << "// To regenerate, build Litmus at the revision above and run\n"
           << "//\n"
           << "//     litmus.x run <config-file>\n"
           << "//\n"
           << "// in the output directory, with <config-file> containing:\n"
           << "//\n";

        for (const auto& line : lines)
        {
            os << "//     " << line << "\n";
        }
    }

    os << "\n";

    return os.str();
}

void
open_generated_file(std::ofstream& fstream, const std::string& fname)
{
    fstream.open(fname.c_str(), std::ios_base::trunc);

    fstream << generated_file_banner();
}

}  // namespace ost
