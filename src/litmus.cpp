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

#include <array>
#include <chrono>
#include <cstddef>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "config.hpp"
#include "run_configuration.hpp"

#include "t2c_cpu_generators.hpp"
#include "t2c_geom_cpu_generators.hpp"
#include "t2c_geom_deriv_cpu_generators.hpp"
#include "t4c_cpu_generators.hpp"
#include "t4c_geom_cpu_generators.hpp"
#include "t4c_diag_cpu_generators.hpp"
#include "t4c_geom_deriv_cpu_generators.hpp"
#include "t4c_geom_hrr_cpu_generators.hpp"
#include "t4c_eri_tree_generators.hpp"
#include "t3c_cpu_generators.hpp"
#include "t3c_geom_cpu_generators.hpp"
#include "t3c_geom_hrr_cpu_generators.hpp"
#include "g2c_cpu_generators.hpp"
#include "t2c_ecp_cpu_generators.hpp"
#include "t2c_geom_ecp_generators.hpp"
#include "t2c_hrr_cpu_generators.hpp"
#include "t2c_proj_ecp_cpu_generators.hpp"
#include "t2c_geom_proj_ecp_cpu_generators.hpp"

#include "two_center_generators.hpp"
#include "spherical_momentum_generators.hpp"

namespace {  // run-driver helpers

/// The run-type families understood by the dispatcher (for help and errors).
const char* const valid_types =
    "t2c_cpu, t2c_hrr_cpu, t2c_geom_cpu, t4c_cpu, t4c_geom_cpu, t4c_geom_hrr_cpu, "
    "t4c_diag_cpu, t4c_call_tree, t3c_cpu, t3c_geom_hrr_cpu, g2c_cpu, t2c_ecp_cpu, "
    "t2c_proj_ecp_cpu, spherical_momentum";

/// Prints command-line usage and the configuration schema.
void
print_usage(std::ostream& os)
{
    os << "Litmus - an automated molecular integrals generator.\n\n"
       << "Usage:\n"
       << "  litmus run <config-file>   Generate integrals described by the config file.\n"
       << "  litmus --help              Show this help.\n\n"
       << "Generated source files are written to the current working directory.\n\n"
       << "Config file (minimal TOML subset: 'key = value', '#' comments). The\n"
       << "schema is chosen per config: an 'integral_type' or 'recursion_type' key\n"
       << "selects the new-style schema, otherwise the legacy 'type' schema is used.\n\n"
       << "Legacy schema (key 'type'):\n"
       << "  type       run-type family (required). One of:\n"
       << "             " << valid_types << "\n"
       << "  lmax       maximum angular momentum (int, default 0). For the\n"
       << "             spherical_momentum type this is the highest shell tabulated.\n"
       << "  integral   integral/operator label (string, default \"none\").\n"
       << "  geom       geometric-derivative orders (int array; arity depends on type:\n"
       << "             3 for t2c/t3c/ecp, 4 for t4c_geom, 5 for t4c).\n"
       << "  aux_lmax   auxiliary angular momentum for t3c types (int, default lmax+2).\n"
       << "  proj_lmax  projector angular momentum for t2c_proj_ecp (int, default 0).\n"
       << "  rec_form   recursion form for t2c types (2-entry int array, default [1, 0]).\n"
       << "  use_rs     range-separation flag for t2c/g2c types (bool, default false).\n\n"
       << "New-style schema (key 'integral_type' or 'recursion_type'; spellings are\n"
       << "case- and separator-insensitive, e.g. 'two_center' == 'TwoCenter'):\n"
       << "  integral_type  integral arity: two_center|2c, three_center|3c,\n"
       << "                 four_center|4c. Only two_center is wired in so far.\n"
       << "  recursion_type horizontal-recurrence transfer: hrr_bra_ket, hrr_bra,\n"
       << "                 hrr_ket. Not wired in yet.\n"
       << "                 (exactly one of integral_type / recursion_type required.)\n"
       << "  max_ang_mom    maximum angular momentum (int, required).\n"
       << "  min_ang_mom    minimum angular momentum (int, default 0).\n"
       << "  operator_type  integrand (default overlap): overlap, kinetic_energy,\n"
       << "                 nuclear_potential, electron_repulsion, dipole_momentum,\n"
       << "                 linear_momentum, local_ecp, projected_ecp,\n"
       << "                 three_center_overlap, three_center_r2,\n"
       << "                 three_center_r_dot_r2. two_center supports overlap,\n"
       << "                 kinetic_energy, electron_repulsion.\n"
       << "  hardware       target hardware (default cpu): cpu.\n"
       << "  language       target language (default C++): C++.\n"
       << "  storage_form   result container (default VeloxChemSparse).\n"
       << "  signature      kernel signature (default VeloxChemScreened).\n";
}

/// Reads the 'geom' key as a fixed-arity array, validating its length.
/// @param config The parsed configuration.
/// @param type The run-type family (for error messages).
/// @return The geometric-derivative orders (defaults to all zeros when absent).
template <std::size_t N>
std::array<int, N>
read_geom(const cfg::Config& config, const std::string& type)
{
    const auto values = config.get_int_array("geom", std::vector<int>(N, 0));

    if (values.size() != N)
    {
        throw cfg::ConfigError("config: run type '" + type + "' expects 'geom' with " +
                               std::to_string(N) + " entries, got " +
                               std::to_string(values.size()));
    }

    std::array<int, N> geom{};

    for (std::size_t i = 0; i < N; i++) geom[i] = values[i];

    return geom;
}

/// Reads the 'rec_form' key as a (bool, bool) pair.
/// @param config The parsed configuration.
/// @return The recursion-form pair (defaults to {true, false} when absent).
std::pair<bool, bool>
read_rec_form(const cfg::Config& config)
{
    const auto values = config.get_int_array("rec_form", std::vector<int>{1, 0});

    if (values.size() != 2)
    {
        throw cfg::ConfigError("config: 'rec_form' must have 2 entries, got " +
                               std::to_string(values.size()));
    }

    return {values[0] != 0, values[1] != 0};
}

/// True if every geometric-derivative order is zero (i.e. a plain integral run).
template <std::size_t N>
bool
is_plain(const std::array<int, N>& geom)
{
    for (const auto order : geom)
    {
        if (order != 0) return false;
    }

    return true;
}

/// Reports a validated run configuration to stdout.
/// @param run_config The validated run configuration.
void
describe(const cfg::RunConfiguration& run_config)
{
    std::cout << "Run configuration:\n";

    if (run_config.integral_type)
    {
        std::cout << "  integral_type = " << cfg::to_string(*run_config.integral_type) << "\n";
    }
    else
    {
        std::cout << "  recursion_type = " << cfg::to_string(*run_config.recursion_type) << "\n";
    }

    std::cout << "  operator_type = " << cfg::to_string(run_config.operator_type) << "\n"
              << "  hardware      = " << cfg::to_string(run_config.hardware) << "\n"
              << "  language      = " << cfg::to_string(run_config.language) << "\n"
              << "  ang_mom       = [" << run_config.min_ang_mom << ", "
              << run_config.max_ang_mom << "]\n"
              << "  storage_form  = " << cfg::to_string(run_config.storage_form) << "\n"
              << "  signature     = " << cfg::to_string(run_config.signature) << "\n";
}

/// Dispatches a parsed configuration to the matching code generator.
/// @param config The parsed configuration.
/// @return The process exit code (0 on success, 1 on an unknown run type).
int
run(const cfg::Config& config)
{
    // new-style configuration: an integral_type or recursion_type key selects the
    // orthogonal-field schema. Not every generator is wired in yet; validate and
    // report.

    if (config.has("integral_type") || config.has("recursion_type"))
    {
        const auto run_config = cfg::make_run_configuration(config);

        describe(run_config);

        if (run_config.recursion_type)
        {
            std::cout << "litmus: configuration is valid; "
                      << cfg::to_string(*run_config.recursion_type)
                      << " recursion generators are not wired in yet." << std::endl;
            return 0;
        }

        switch (*run_config.integral_type)
        {
            case cfg::IntegralType::two_center:
                TwoCenterGenerator().generate(run_config);
                return 0;

            case cfg::IntegralType::three_center:
            case cfg::IntegralType::four_center:
                std::cout << "litmus: configuration is valid; "
                          << cfg::to_string(*run_config.integral_type)
                          << " generators are not wired in yet." << std::endl;
                return 0;
        }
    }

    const auto type = config.get_string("type");

    const auto integral = config.get_string("integral", "none");

    const auto lmax = config.get_int("lmax", 0);

    // two-center integrals

    if (type == "t2c_cpu")
    {
        const auto geom = read_geom<3>(config, type);

        const auto rec_form = read_rec_form(config);

        const auto use_rs = config.get_bool("use_rs", false);

        if ((geom[0] + geom[2]) == 0)
        {
            T2CCPUGenerator().generate(integral, lmax, geom, rec_form, use_rs);
        }
        else
        {
            T2CGeomCPUGenerator().generate(integral, lmax, geom, rec_form, use_rs);
        }

        return 0;
    }

    if (type == "t2c_hrr_cpu")
    {
        T2CHRRCPUGenerator().generate(lmax);

        return 0;
    }

    if (type == "t2c_geom_cpu")
    {
        T2CGeomDerivCPUGenerator().generate(lmax, read_geom<3>(config, type));

        return 0;
    }

    // four-center integrals

    if (type == "t4c_cpu")
    {
        const auto geom = read_geom<5>(config, type);

        if (is_plain(geom))
        {
            T4CCPUGenerator().generate(integral, lmax);
        }
        else
        {
            T4CGeomCPUGenerator().generate(integral, lmax, geom);
        }

        return 0;
    }

    if (type == "t4c_geom_cpu")
    {
        T4CGeomDerivCPUGenerator().generate(lmax, read_geom<4>(config, type));

        return 0;
    }

    if (type == "t4c_geom_hrr_cpu")
    {
        T4CGeomHrrCPUGenerator().generate(integral, lmax, read_geom<4>(config, type));

        return 0;
    }

    if (type == "t4c_diag_cpu")
    {
        T4CDiagCPUGenerator().generate(integral, lmax);

        return 0;
    }

    if (type == "t4c_call_tree")
    {
        T4CCallTreeGenerator().generate(integral, lmax);

        return 0;
    }

    // three-center integrals

    if (type == "t3c_cpu")
    {
        const auto geom = read_geom<3>(config, type);

        const auto aux_lmax = config.get_int("aux_lmax", lmax + 2);

        if (is_plain(geom))
        {
            T3CCPUGenerator().generate(integral, lmax, aux_lmax);
        }
        else
        {
            T3CGeomCPUGenerator().generate(integral, lmax, aux_lmax, geom);
        }

        return 0;
    }

    if (type == "t3c_geom_hrr_cpu")
    {
        const auto aux_lmax = config.get_int("aux_lmax", lmax + 2);

        T3CGeomHrrCPUGenerator().generate(integral, aux_lmax, read_geom<3>(config, type));

        return 0;
    }

    // two-center integrals on a grid

    if (type == "g2c_cpu")
    {
        const auto geom = read_geom<3>(config, type);

        const auto use_rs = config.get_bool("use_rs", false);

        if (is_plain(geom))
        {
            G2CCPUGenerator().generate(integral, lmax, geom, use_rs);

            return 0;
        }

        throw cfg::ConfigError("config: geometric-derivative g2c integrals are not implemented");
    }

    // two-center ECP integrals

    if (type == "t2c_ecp_cpu")
    {
        const auto geom = read_geom<3>(config, type);

        if (is_plain(geom))
        {
            T2CECPCPUGenerator().generate(integral, lmax);
        }
        else
        {
            T2CECPGeomCPUGenerator().generate(integral, lmax, geom);
        }

        return 0;
    }

    if (type == "t2c_proj_ecp_cpu")
    {
        const auto geom = read_geom<3>(config, type);

        const auto proj_lmax = config.get_int("proj_lmax", 0);

        if (is_plain(geom))
        {
            T2CProjECPCPUGenerator().generate(integral, lmax, proj_lmax);
        }
        else
        {
            T2CGeomProjECPCPUGenerator().generate(integral, lmax, proj_lmax, geom);
        }

        return 0;
    }

    // Cartesian-to-spherical transformation factors (VeloxChem SphericalMomentum.hpp)

    if (type == "spherical_momentum")
    {
        SphericalMomentumGenerator().generate(lmax);

        return 0;
    }

    std::cerr << "litmus: unknown run type '" << type << "'.\n"
              << "valid types: " << valid_types << std::endl;

    return 1;
}

}  // namespace

int
main(int argc, char** argv)
{
    const std::vector<std::string> args(argv + 1, argv + argc);

    if (args.empty() || (args[0] == "-h") || (args[0] == "--help"))
    {
        print_usage(args.empty() ? std::cerr : std::cout);

        return args.empty() ? 1 : 0;
    }

    if ((args[0] != "run") || (args.size() != 2))
    {
        std::cerr << "litmus: expected 'litmus run <config-file>'.\n\n";

        print_usage(std::cerr);

        return 1;
    }

    try
    {
        const auto config = cfg::parse_file(args[1]);

        const auto stime = std::chrono::high_resolution_clock::now();

        const auto rc = run(config);

        const auto etime = std::chrono::high_resolution_clock::now();

        const auto dtime = std::chrono::duration_cast<std::chrono::seconds>(etime - stime);

        std::cout << "Elapsed time: " << dtime.count() << " seconds." << std::endl;

        return rc;
    }
    catch (const cfg::ConfigError& error)
    {
        std::cerr << "litmus: " << error.what() << std::endl;

        return 1;
    }
}
