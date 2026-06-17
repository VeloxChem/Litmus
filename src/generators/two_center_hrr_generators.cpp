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

#include "two_center_hrr_generators.hpp"

#include <fstream>
#include <iostream>
#include <string>

#include "tensor.hpp"
#include "two_center_hrr_emitter.hpp"

namespace {  // two-center HRR generator helpers

/// The base file name (no extension) of a (la|lb) kernel, e.g. (1|1) -> "...PP".
std::string
kernel_file_name(const int la, const int lb)
{
    return "ObaraSaikaTwoCenterHrr" + Tensor(la).label() + Tensor(lb).label();
}

/// Whether a (la|lb) target is emitted for the requested recursion type. HRR needs
/// both sides >= 1; hrr_bra keeps the bra-incremented targets (la <= lb), hrr_ket
/// the ket-incremented ones (la > lb), hrr_bra_ket all of them.
bool
selected(const cfg::RecursionType type, const int la, const int lb)
{
    if (la < 1 || lb < 1) return false;

    switch (type)
    {
        case cfg::RecursionType::hrr_bra_ket:
            return true;

        case cfg::RecursionType::hrr_bra:
            return la <= lb;

        case cfg::RecursionType::hrr_ket:
            return la > lb;
    }

    return false;  // unreachable: every RecursionType is handled above
}

/// Writes the kernel declaration header (.hpp).
void
write_hpp(const int la, const int lb)
{
    const auto base = kernel_file_name(la, lb);

    const auto guard = base + "_hpp";

    std::ofstream fstream;

    fstream.open((base + ".hpp").c_str(), std::ios_base::trunc);

    fstream << "#ifndef " << guard << "\n";
    fstream << "#define " << guard << "\n\n";
    fstream << "#include \"Array.hpp\"\n\n";
    fstream << "namespace os2c::hrr {  // horizontal recurrence\n\n";
    fstream << format_hrr_signature(la, lb) << ";\n\n";
    fstream << "}  // namespace os2c::hrr\n\n";
    fstream << "#endif /* " << guard << " */\n";

    fstream.close();
}

/// Writes the kernel definition (.cpp).
void
write_cpp(const int la, const int lb)
{
    const auto base = kernel_file_name(la, lb);

    std::ofstream fstream;

    fstream.open((base + ".cpp").c_str(), std::ios_base::trunc);

    fstream << "#include \"" << base << ".hpp\"\n\n";
    fstream << "#include <cmath>\n\n";
    fstream << "namespace os2c::hrr {  // horizontal recurrence\n\n";
    fstream << format_hrr_kernel(la, lb) << "\n";
    fstream << "}  // namespace os2c::hrr\n";

    fstream.close();
}

}  // namespace

void
TwoCenterHrrGenerator::generate(const cfg::RunConfiguration& run_config) const
{
    const auto type = *run_config.recursion_type;

    int count = 0;

    for (int la = run_config.min_ang_mom; la <= run_config.max_ang_mom; la++)
    {
        for (int lb = run_config.min_ang_mom; lb <= run_config.max_ang_mom; lb++)
        {
            if (!selected(type, la, lb)) continue;

            write_hpp(la, lb);

            write_cpp(la, lb);

            std::cout << "Generated " << kernel_file_name(la, lb) << " kernel" << std::endl;

            count++;
        }
    }

    std::cout << "Generated " << count << " " << cfg::to_string(type)
              << " two-center HRR kernels." << std::endl;
}
