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

#include "two_center_vrr_generators.hpp"

#include <fstream>
#include <iostream>
#include <string>

#include "tensor.hpp"
#include "two_center_vrr_emitter.hpp"

namespace {  // two-center VRR generator helpers

/// The naming/emission traits of a VRR flavor: the file-name tag, the kernel
/// namespace and its documentation, whether the source needs <cmath>, and the
/// signature/definition emitters.
struct VrrFlavor
{
    std::string file_tag;   // "Cart" or "Sph"
    std::string ns;         // "os2c::vrr::ovl" or "os2c::ovl"
    std::string caption;    // namespace comment
    bool        needs_math;  // the spherical kernels use std::sqrt

    std::string (*signature)(const int);
    std::string (*kernel)(const int);
};

/// The flavor traits for the configured recursion type. The switch carries no
/// default, so a new RecursionType trips -Wswitch here.
/// @param type The recursion type (a vrr_* value).
/// @return The flavor traits.
VrrFlavor
flavor(const cfg::RecursionType type)
{
    switch (type)
    {
        case cfg::RecursionType::vrr_cartesian:
            return {"Cart", "os2c::vrr::ovl", "overlap Cartesian vertical recurrence",
                    false, format_vrr_cartesian_signature, format_vrr_cartesian_kernel};

        case cfg::RecursionType::vrr_spherical:
            return {"Sph", "os2c::ovl", "overlap spherical vertical recurrence",
                    true, format_vrr_spherical_signature, format_vrr_spherical_kernel};

        // the HRR recursion types are handled by the HRR generator, not here.
        case cfg::RecursionType::hrr_bra_ket:
        case cfg::RecursionType::hrr_bra:
        case cfg::RecursionType::hrr_ket:
            break;
    }

    throw cfg::ConfigError("two-center VRR generator: recursion_type '" + cfg::to_string(type) +
                           "' is not a vertical recurrence");
}

/// The base file name (no extension) of a (s|lb) kernel, e.g. lb = 2, "Cart" ->
/// "ObaraSaikaTwoCenterOverlapVrrCartD".
std::string
kernel_file_name(const VrrFlavor& flv, const int lb)
{
    return "ObaraSaikaTwoCenterOverlapVrr" + flv.file_tag + Tensor(lb).label();
}

/// Writes the kernel declaration header (.hpp).
void
write_hpp(const VrrFlavor& flv, const int lb)
{
    const auto base = kernel_file_name(flv, lb);

    const auto guard = base + "_hpp";

    std::ofstream fstream;

    fstream.open((base + ".hpp").c_str(), std::ios_base::trunc);

    fstream << "#ifndef " << guard << "\n";
    fstream << "#define " << guard << "\n\n";
    fstream << "#include \"Array.hpp\"\n";
    fstream << "#include \"BasisFunctionPair.hpp\"\n\n";
    fstream << "namespace " << flv.ns << " {  // " << flv.caption << "\n\n";
    fstream << flv.signature(lb) << ";\n\n";
    fstream << "}  // namespace " << flv.ns << "\n\n";
    fstream << "#endif /* " << guard << " */\n";

    fstream.close();
}

/// Writes the kernel definition (.cpp).
void
write_cpp(const VrrFlavor& flv, const int lb)
{
    const auto base = kernel_file_name(flv, lb);

    std::ofstream fstream;

    fstream.open((base + ".cpp").c_str(), std::ios_base::trunc);

    fstream << "#include \"" << base << ".hpp\"\n\n";

    if (flv.needs_math) fstream << "#include <cmath>\n\n";

    fstream << "namespace " << flv.ns << " {  // " << flv.caption << "\n\n";
    fstream << flv.kernel(lb) << "\n";
    fstream << "}  // namespace " << flv.ns << "\n";

    fstream.close();
}

}  // namespace

void
TwoCenterVrrGenerator::generate(const cfg::RunConfiguration& run_config) const
{
    const auto flv = flavor(*run_config.recursion_type);

    // the vertical recurrence builds the ket, so lb must be at least 1.
    const auto min_lb = (run_config.min_ang_mom < 1) ? 1 : run_config.min_ang_mom;

    int count = 0;

    for (int lb = min_lb; lb <= run_config.max_ang_mom; lb++)
    {
        write_hpp(flv, lb);

        write_cpp(flv, lb);

        std::cout << "Generated " << kernel_file_name(flv, lb) << " kernel" << std::endl;

        count++;
    }

    std::cout << "Generated " << count << " " << cfg::to_string(*run_config.recursion_type)
              << " two-center VRR kernels." << std::endl;
}
