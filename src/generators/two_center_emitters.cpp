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

#include "two_center_emitters.hpp"

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include "config.hpp"
#include "file_stream.hpp"
#include "operator.hpp"
#include "string_formater.hpp"
#include "tensor.hpp"

namespace {  // C++/CPU emitter helpers

/// The naming tags for an integrand operator: the (C++17 nested) namespace the
/// kernel lives in and the CamelCase label used in its file name.
struct OperatorTags
{
    /// The fully-qualified kernel namespace, e.g. "os2c::ovl".
    std::string ns;

    /// The CamelCase operator label used in file names, e.g. "Overlap".
    std::string file_label;

    /// The human-readable operator caption for documentation, e.g. "overlap".
    std::string caption;
};

/// Maps an integrand operator to its naming tags. Only the operators the
/// two-center generator currently supports are named; the rest fall through to a
/// thrown error. The switch carries no default, so a new OperatorType trips
/// -Wswitch here.
/// @param op The integrand operator type.
/// @return The naming tags for the operator.
OperatorTags
operator_tags(cfg::OperatorType op)
{
    switch (op)
    {
        case cfg::OperatorType::overlap:
            return {"os2c::ovl", "Overlap", "overlap"};

        case cfg::OperatorType::kinetic_energy:
            return {"os2c::kin", "KineticEnergy", "kinetic energy"};

        case cfg::OperatorType::electron_repulsion:
            return {"os2c::eri", "ElectronRepulsion", "electron repulsion"};

        // remaining operators have no two-center kernel yet
        case cfg::OperatorType::nuclear_potential:
        case cfg::OperatorType::dipole_momentum:
        case cfg::OperatorType::linear_momentum:
        case cfg::OperatorType::three_center_overlap:
        case cfg::OperatorType::three_center_r2:
        case cfg::OperatorType::three_center_r_dot_r2:
        case cfg::OperatorType::local_ecp:
        case cfg::OperatorType::projected_ecp:
            break;
    }

    throw cfg::ConfigError("two-center emitter: operator '" + cfg::to_string(op) +
                           "' has no two-center kernel naming");
}

/// The lowercase spectroscopic shell label (s, p, d, f, ...) of an angular
/// momentum, as used in the compute_<la>_<lb> function name.
/// @param ang_mom The angular momentum.
/// @return The lowercase shell label.
std::string
shell_label(int ang_mom)
{
    return fstr::lowercase(Tensor(ang_mom).label());
}

/// The number of Cartesian components of a shell of angular momentum l,
/// (l + 1)(l + 2)/2. This is the per-pair row count of a contracted Cartesian
/// integral block (s -> 1, p -> 3, d -> 6, ...).
/// @param l The angular momentum.
/// @return The Cartesian component count.
int
cartesian_count(int l)
{
    return (l + 1) * (l + 2) / 2;
}

/// The kernel function name for a target integral, e.g. "compute_p_p".
/// @param integral The target two-center integral.
/// @return The function name.
std::string
kernel_func_name(const I2CIntegral& integral)
{
    return "compute_" + shell_label(integral[0]) + "_" + shell_label(integral[1]);
}

/// The base file name (no extension) for a target integral's kernel, e.g.
/// "ObaraSaikaTwoCenterOverlapPP"; the .hpp/.cpp pair share it.
/// @param op The integrand operator type.
/// @param integral The target two-center integral.
/// @return The base file name.
std::string
kernel_file_name(cfg::OperatorType op, const I2CIntegral& integral)
{
    return "ObaraSaikaTwoCenter" + operator_tags(op).file_label + Tensor(integral[0]).label() +
           Tensor(integral[1]).label();
}

/// The return type of a kernel for a storage form. The storage form selects the
/// container the integral block is returned in. The switch carries no default,
/// so a new StorageForm trips -Wswitch here.
/// @param form The storage form.
/// @return The return type name.
std::string
return_type(cfg::StorageForm form)
{
    switch (form)
    {
        case cfg::StorageForm::veloxchem_sparse:
            return "osfunc::CArray<double>";
    }

    return std::string();  // unreachable: every StorageForm is handled above
}

/// The input parameters of a kernel for a signature, as "type name" fragments.
/// The signature selects the calling convention. The switch carries no default,
/// so a new Signature trips -Wswitch here.
/// @param signature The kernel signature convention.
/// @return The input parameter fragments.
std::vector<std::string>
input_params(cfg::Signature signature)
{
    switch (signature)
    {
        case cfg::Signature::veloxchem_screened:
            return {"const osfunc::CBasisFunctionPair& pair"};
    }

    return {};  // unreachable: every Signature is handled above
}

/// The header includes a kernel needs for its return and input types.
/// @param run_config The run configuration (selects storage form and signature).
/// @return The include lines (with quotes).
std::vector<std::string>
kernel_includes(const cfg::RunConfiguration& run_config)
{
    std::vector<std::string> vstr;

    switch (run_config.storage_form)
    {
        case cfg::StorageForm::veloxchem_sparse:
            vstr.push_back("#include \"Array.hpp\"");
            break;
    }

    switch (run_config.signature)
    {
        case cfg::Signature::veloxchem_screened:
            vstr.push_back("#include \"BasisFunctionPair.hpp\"");
            break;
    }

    return vstr;
}

/// Formats a two-center integral as a "(bra|operator|ket)" caption, appending a
/// non-zero recursion order as "^n", for use in generated comments.
/// @param integral The two-center integral.
/// @return The caption.
std::string
integral_caption(const I2CIntegral& integral)
{
    auto text = "(" + Tensor(integral[0]).label() + "|" + integral.integrand().name() +
                "|" + Tensor(integral[1]).label() + ")";

    if (const auto order = integral.order(); order != 0)
    {
        text += "^" + std::to_string(order);
    }

    return text;
}

/// The two-center emitter for the C++ language on CPU hardware. Produces a
/// header/definition pair whose body lays out the integral computation workflow
/// (VRR seeds -> VRR closure -> HRR transfer -> store), shaped by the configured
/// signature (inputs) and storage form (return).
class CppCpuTwoCenterEmitter : public TwoCenterEmitter
{
    /// Writes the kernel declaration header (.hpp).
    void _write_hpp(const cfg::RunConfiguration& run_config,
                    const I2CIntegral&           integral) const;

    /// Writes the kernel definition (.cpp) carrying the computation workflow.
    void _write_cpp(const cfg::RunConfiguration& run_config,
                    const I2CIntegral&           integral,
                    const SI2CIntegrals&         hrr_ints,
                    const SI2CIntegrals&         vrr_base_ints,
                    const SI2CIntegrals&         vrr_rest_ints) const;

    /// The kernel signature, as code lines: "compute_<la>_<lb>(<inputs>) ->
    /// <return>", broken across lines and aligned under the function name when
    /// there is more than one input parameter.
    /// @param run_config The run configuration (selects inputs and return type).
    /// @param integral The target integral (names the function).
    /// @param terminus Whether to terminate with a ';' (declaration).
    /// @return The signature lines.
    std::vector<std::string> _signature_lines(const cfg::RunConfiguration& run_config,
                                               const I2CIntegral&           integral,
                                               const bool                   terminus) const;

public:
    void emit(const cfg::RunConfiguration& run_config,
              const I2CIntegral&           integral,
              const SI2CIntegrals&         hrr_ints,
              const SI2CIntegrals&         vrr_base_ints,
              const SI2CIntegrals&         vrr_rest_ints) const override;
};

std::vector<std::string>
CppCpuTwoCenterEmitter::_signature_lines(const cfg::RunConfiguration& run_config,
                                         const I2CIntegral&           integral,
                                         const bool                   terminus) const
{
    std::vector<std::string> vstr;

    const auto name = kernel_func_name(integral) + "(";

    const auto spacer = std::string(name.size(), ' ');

    const auto params = input_params(run_config.signature);

    const auto tail = ") -> " + return_type(run_config.storage_form) + (terminus ? ";" : "");

    for (std::size_t i = 0; i < params.size(); i++)
    {
        const auto head = (i == 0) ? name : spacer;

        const auto last = (i + 1 == params.size());

        vstr.push_back(head + params[i] + (last ? tail : ","));
    }

    return vstr;
}

void
CppCpuTwoCenterEmitter::_write_hpp(const cfg::RunConfiguration& run_config,
                                   const I2CIntegral&           integral) const
{
    const auto tags = operator_tags(run_config.operator_type);

    const auto base = kernel_file_name(run_config.operator_type, integral);

    const auto guard = base + "_hpp";

    std::ofstream fstream;

    fstream.open((base + ".hpp").c_str(), std::ios_base::trunc);

    auto lines = VCodeLines();

    lines.push_back({0, 0, 1, "#ifndef " + guard});
    lines.push_back({0, 0, 2, "#define " + guard});

    for (const auto& include : kernel_includes(run_config))
    {
        lines.push_back({0, 0, 1, include});
    }
    lines.push_back({0, 0, 1, ""});

    lines.push_back({0, 0, 2, "namespace " + tags.ns + " {  // " + tags.caption +
                                  " two-center integrals"});

    lines.push_back({0, 0, 1, "/// @brief Computes " + integral_caption(integral) +
                                  " integrals for a screened pair of basis functions."});
    lines.push_back({0, 0, 1, "/// @param pair The screened pair of basis functions."});
    lines.push_back({0, 0, 1, "/// @return The matrix of computed integrals."});

    lines.push_back({0, 0, 1, "auto"});

    for (const auto& label : _signature_lines(run_config, integral, true))
    {
        lines.push_back({0, 0, 1, label});
    }

    lines.push_back({0, 0, 2, "}  // namespace " + tags.ns});

    lines.push_back({0, 0, 1, "#endif /* " + guard + " */"});

    ost::write_code_lines(fstream, lines);

    fstream.close();
}

void
CppCpuTwoCenterEmitter::_write_cpp(const cfg::RunConfiguration& run_config,
                                   const I2CIntegral&           integral,
                                   const SI2CIntegrals&         hrr_ints,
                                   const SI2CIntegrals&         vrr_base_ints,
                                   const SI2CIntegrals&         vrr_rest_ints) const
{
    const auto tags = operator_tags(run_config.operator_type);

    const auto base = kernel_file_name(run_config.operator_type, integral);

    const int la = integral[0];

    const int lb = integral[1];

    const int nspher = (2 * la + 1) * (2 * lb + 1);

    const bool has_hrr = (la > 0) && (lb > 0);

    // the vertical recurrence builds the smaller side; PB grows the ket (la <= lb),
    // PA grows the bra. The generic Pc kernel argument is fed the matching array.
    const bool ket_built = (la <= lb);

    const auto dist = ket_built ? std::string("pb") : std::string("pa");

    // order the VRR ladder (all built integrals) and the HRR-consumed base by total
    // angular momentum, lowest first.
    const auto by_l = [](const I2CIntegral& a, const I2CIntegral& b) {
        if ((a[0] + a[1]) != (b[0] + b[1])) return (a[0] + a[1]) < (b[0] + b[1]);

        return a < b;
    };

    std::vector<I2CIntegral> vrr_ordered(vrr_base_ints.begin(), vrr_base_ints.end());

    vrr_ordered.insert(vrr_ordered.end(), vrr_rest_ints.begin(), vrr_rest_ints.end());

    std::sort(vrr_ordered.begin(), vrr_ordered.end(), by_l);

    std::vector<I2CIntegral> base_ordered(vrr_base_ints.begin(), vrr_base_ints.end());

    std::sort(base_ordered.begin(), base_ordered.end(), by_l);

    // build the body and the set of os2c kernel headers it calls.
    auto body = VCodeLines();

    std::set<std::string> headers;

    body.push_back({1, 0, 1, "// number of screened atom pairs"});
    body.push_back({1, 0, 1, "const auto npairs = pair.number_of_pairs();"});
    body.push_back({0, 0, 1, ""});

    body.push_back({1, 0, 1, "// spherical (2*la+1) x (2*lb+1) result, one atom pair per column"});
    body.push_back({1, 0, 1, "osfunc::CArray<double> buffer(" + std::to_string(nspher) +
                                 ", npairs);"});
    body.push_back({0, 0, 1, ""});

    if (run_config.operator_type != cfg::OperatorType::overlap)
    {
        // only the overlap VRR/HRR kernels are generated so far; other operators
        // return the zeroed buffer until their kernels exist.
        body.push_back({1, 0, 1, "// TODO: " + tags.caption + " two-center kernels not generated yet"});
    }
    else if (la == 0 && lb == 0)
    {
        // (s|s): the contracted primitive overlaps are the spherical result.
        body.push_back({1, 0, 1, "// (s|s): the contracted primitive overlaps are the result"});
        body.push_back({1, 0, 1, "osfunc::contract(buffer, osfunc::compute_overlap(pair));"});
    }
    else
    {
        body.push_back({1, 0, 1, "// primitive (s|s) seed and Pc (" + dist + ") distances"});
        body.push_back({1, 0, 1, "const auto ss = osfunc::compute_overlap(pair);"});
        body.push_back({1, 0, 1, "const auto " + dist + " = osfunc::compute_" + dist + "(pair);"});
        body.push_back({0, 0, 1, ""});

        if (!has_hrr)
        {
            // (s|lb) or (la|s): the spherical VRR fuses the recurrence, contraction
            // and Cartesian-to-spherical transform.
            const auto lval = la + lb;

            body.push_back({1, 0, 1, "// vertical recurrence with the Cartesian-to-spherical transform"});
            body.push_back({1, 0, 1, tags.ns + "::compute_" + shell_label(lval) + "_sph(pair, ss, " +
                                         dist + ", buffer);"});

            headers.insert("ObaraSaikaTwoCenterOverlapVrrSph" + Tensor(lval).label() + ".hpp");
        }
        else
        {
            body.push_back({1, 0, 1, "// number of primitive pairs"});
            body.push_back({1, 0, 1, "const auto nprims = pair.number_of_primitive_pairs();"});
            body.push_back({0, 0, 1, ""});

            // vertical recurrence: build the primitive Cartesian ladder up to the
            // HRR base, each step from the two below it.
            body.push_back({1, 0, 1, "// vertical recurrence: primitive Cartesian base integrals"});

            for (const auto& tint : vrr_ordered)
            {
                const auto lval = tint[0] + tint[1];

                if (lval == 0) continue;  // the (s|s) seed is "ss"

                const auto kb = (tint[0] == 0);

                const auto name = shell_label(tint[0]) + shell_label(tint[1]);

                const auto rows = cartesian_count(lval);

                std::string args = "pair, ";

                for (int m = lval - 1; (m >= 0) && (m >= lval - 2); m--)
                {
                    args += (kb ? ("s" + shell_label(m)) : (shell_label(m) + "s")) + ", ";
                }

                body.push_back({1, 0, 1, "osfunc::CArray<double> " + name + "(nprims * " +
                                             std::to_string(rows) + ", npairs);"});
                body.push_back({1, 0, 1, "os2c::vrr::ovl::compute_" + shell_label(lval) + "(" + args +
                                             dist + ", " + name + ");"});

                headers.insert("ObaraSaikaTwoCenterOverlapVrrCart" + Tensor(lval).label() + ".hpp");
            }

            body.push_back({0, 0, 1, ""});

            // contract the base integrals the HRR consumes.
            body.push_back({1, 0, 1, "// contract the base integrals consumed by the horizontal recurrence"});

            for (const auto& tint : base_ordered)
            {
                const auto name = shell_label(tint[0]) + shell_label(tint[1]);

                const auto rows = cartesian_count(tint[0]) * cartesian_count(tint[1]);

                body.push_back({1, 0, 1, "osfunc::CArray<double> c" + name + "(" +
                                             std::to_string(rows) + ", npairs);"});
                body.push_back({1, 0, 1, "osfunc::contract(c" + name + ", " + name + ");"});
            }

            body.push_back({0, 0, 1, ""});

            // horizontal recurrence with the transform fused; it transfers momentum
            // between the centers using the AB distances.
            std::string args = "";

            for (const auto& tint : base_ordered)
            {
                args += "c" + shell_label(tint[0]) + shell_label(tint[1]) + ", ";
            }

            body.push_back({1, 0, 1, "// horizontal recurrence (Cartesian-to-spherical transform fused)"});
            body.push_back({1, 0, 1, "const auto ab = osfunc::compute_ab(pair);"});
            body.push_back({1, 0, 1, "os2c::hrr::compute_" + shell_label(la) + "_" + shell_label(lb) +
                                         "(" + args + "ab, buffer);"});

            headers.insert("ObaraSaikaTwoCenterHrr" + Tensor(la).label() + Tensor(lb).label() + ".hpp");
        }
    }

    body.push_back({0, 0, 1, ""});
    body.push_back({1, 0, 1, "return buffer;"});

    // assemble the file: includes, namespace, signature, body.
    std::ofstream fstream;

    fstream.open((base + ".cpp").c_str(), std::ios_base::trunc);

    auto lines = VCodeLines();

    lines.push_back({0, 0, 2, "#include \"" + base + ".hpp\""});

    lines.push_back({0, 0, 1, "#include \"ObaraSaikaFunc.hpp\""});

    for (const auto& header : headers) lines.push_back({0, 0, 1, "#include \"" + header + "\""});

    lines.push_back({0, 0, 1, ""});

    lines.push_back({0, 0, 2, "namespace " + tags.ns + " {  // " + tags.caption +
                                  " two-center integrals"});

    lines.push_back({0, 0, 1, "auto"});

    for (const auto& label : _signature_lines(run_config, integral, false))
    {
        lines.push_back({0, 0, 1, label});
    }

    lines.push_back({0, 0, 1, "{"});

    for (const auto& line : body) lines.push_back(line);

    lines.push_back({0, 0, 2, "}"});

    lines.push_back({0, 0, 1, "}  // namespace " + tags.ns});

    ost::write_code_lines(fstream, lines);

    fstream.close();
}

void
CppCpuTwoCenterEmitter::emit(const cfg::RunConfiguration& run_config,
                             const I2CIntegral&           integral,
                             const SI2CIntegrals&         hrr_ints,
                             const SI2CIntegrals&         vrr_base_ints,
                             const SI2CIntegrals&         vrr_rest_ints) const
{
    _write_hpp(run_config, integral);

    _write_cpp(run_config, integral, hrr_ints, vrr_base_ints, vrr_rest_ints);
}

}  // namespace

std::unique_ptr<TwoCenterEmitter>
make_two_center_emitter(const cfg::RunConfiguration& run_config)
{
    switch (run_config.hardware)
    {
        case cfg::Hardware::cpu:
        {
            switch (run_config.language)
            {
                case cfg::Language::cpp:
                    return std::make_unique<CppCpuTwoCenterEmitter>();
            }

            break;
        }
    }

    throw cfg::ConfigError("two-center generator: no emitter for hardware '" +
                           cfg::to_string(run_config.hardware) + "' and language '" +
                           cfg::to_string(run_config.language) + "'");
}
