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

/// The namespace abbreviation for a VRR integral's integrand, as used in the
/// os2c::vrr::<abbrev> kernel namespace and the p<abbrev>_<shells> result name.
/// Matches the operator_tags abbreviations but keys off the algebra integrand of
/// an individual VRR integral (kinetic recurrences also produce overlap
/// auxiliaries, so the integrand can differ from the kernel's operator).
/// @param integrand The integrand operator of the VRR integral.
/// @return The abbreviation ("ovl", "kin", "eri").
std::string
integrand_abbrev(const Operator& integrand)
{
    const auto name = integrand.name();

    if (name == "1") return "ovl";

    if (name == "T") return "kin";

    if (name == "1/|r-r'|") return "eri";

    throw cfg::ConfigError("two-center emitter: no VRR namespace for integrand '" + name + "'");
}

/// Whether a VRR integral is the elementary primitive-overlap seed, i.e. the
/// (s|1|s) overlap already produced by osfunc::compute_overlap (povl_ss). Such an
/// integral is not re-computed by a VRR kernel.
/// @param integral The VRR integral.
/// @return True if it is the primitive-overlap seed.
bool
is_primitive_overlap_seed(const I2CIntegral& integral)
{
    return (integral.integrand().name() == "1") && (integral[0] == 0) && (integral[1] == 0);
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
            return "CDenseMatrix";
    }

    return std::string();  // unreachable: every StorageForm is handled above
}

/// The body statements that allocate and zero the result buffer for a storage
/// form. The buffer is sized (2*la+1)*(2*lb+1) rows (spherical components of the
/// bra/ket shells) by one column per screened atom pair (the local `npairs`).
/// The switch carries no default, so a new StorageForm trips -Wswitch here.
/// @param form The storage form.
/// @param nrows The number of result rows, (2*la+1)*(2*lb+1).
/// @return The buffer-allocation statements.
std::vector<std::string>
result_buffer_lines(cfg::StorageForm form, int nrows)
{
    switch (form)
    {
        case cfg::StorageForm::veloxchem_sparse:
            return {"// result buffer: (2*la+1)*(2*lb+1) spherical components per row,",
                    "// one screened atom pair per column (zero-initialized by construction)",
                    "auto buffer = CDenseMatrix(" + std::to_string(nrows) +
                        ", static_cast<int>(npairs));"};
    }

    return {};  // unreachable: every StorageForm is handled above
}

/// Whether an operator's vertical recurrence is seeded by the elementary
/// primitive overlaps (osfunc::compute_overlap). True for overlap and kinetic
/// energy; electron repulsion seeds from the Boys function instead. The switch
/// carries no default, so a new OperatorType trips -Wswitch here.
/// @param op The integrand operator type.
/// @return True if the recurrence seeds from the primitive overlaps.
bool
seeds_with_primitive_overlap(cfg::OperatorType op)
{
    switch (op)
    {
        case cfg::OperatorType::overlap:
        case cfg::OperatorType::kinetic_energy:
            return true;

        case cfg::OperatorType::electron_repulsion:
        case cfg::OperatorType::nuclear_potential:
        case cfg::OperatorType::dipole_momentum:
        case cfg::OperatorType::linear_momentum:
        case cfg::OperatorType::three_center_overlap:
        case cfg::OperatorType::three_center_r2:
        case cfg::OperatorType::three_center_r_dot_r2:
        case cfg::OperatorType::local_ecp:
        case cfg::OperatorType::projected_ecp:
            return false;
    }

    return false;  // unreachable: every OperatorType is handled above
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
            return {"const CScreenedBasisFunctionPair& pair"};
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
            vstr.push_back("#include \"DenseMatrix.hpp\"");
            break;
    }

    switch (run_config.signature)
    {
        case cfg::Signature::veloxchem_screened:
            vstr.push_back("#include \"ScreenedBasisFunctionPair.hpp\"");
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

/// Which side the vertical recurrence builds up. The VRR seeds keep the smaller
/// side at zero and ladder the other: (0|o|b) seeds grow the ket, (b|o|0) seeds
/// grow the bra. A pure (0|o|0) target grows neither.
enum class VrrBuildSide
{
    none,
    bra,
    ket
};

/// Reads the VRR build side off a set of VRR integrals: the ket side if any
/// integral is (0|o|b) with b > 0, the bra side if any is (b|o|0) with b > 0, and
/// neither otherwise (a (0|o|0)-only target).
/// @param vrr_ints The VRR integrals (base seeds and their closure).
/// @return The side the recurrence builds up.
VrrBuildSide
vrr_build_side(const SI2CIntegrals& vrr_ints)
{
    for (const auto& tint : vrr_ints)
    {
        if ((tint[0] == 0) && (tint[1] > 0)) return VrrBuildSide::ket;

        if ((tint[0] > 0) && (tint[1] == 0)) return VrrBuildSide::bra;
    }

    return VrrBuildSide::none;
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

    std::ofstream fstream;

    fstream.open((base + ".cpp").c_str(), std::ios_base::trunc);

    auto lines = VCodeLines();

    lines.push_back({0, 0, 1, "#include \"" + base + ".hpp\""});

    // the Cartesian-to-spherical transform is only emitted off the (s|s) special
    // case, so its header is only included when a transform call is generated.
    const auto emits_transform = seeds_with_primitive_overlap(run_config.operator_type) &&
                                 !is_primitive_overlap_seed(integral);

    lines.push_back({0, 0, emits_transform ? 1 : 2, "#include \"ObaraSaikaFunc.hpp\""});

    if (emits_transform)
    {
        lines.push_back({0, 0, 2, "#include \"CartesianToSphericalFunc.hpp\""});
    }

    lines.push_back({0, 0, 2, "namespace " + tags.ns + " {  // " + tags.caption +
                                  " two-center integrals"});

    lines.push_back({0, 0, 1, "auto"});

    for (const auto& label : _signature_lines(run_config, integral, false))
    {
        lines.push_back({0, 0, 1, label});
    }

    lines.push_back({0, 0, 1, "{"});

    // the screened pair carries the bra/ket basis functions and the surviving
    // atom pairs.
    lines.push_back({1, 0, 1, "// number of screened atom pairs to evaluate"});
    lines.push_back({1, 0, 1, "const auto npairs = pair.number_of_pairs();"});
    lines.push_back({0, 0, 1, ""});

    const auto nrows = (2 * integral[0] + 1) * (2 * integral[1] + 1);

    for (const auto& label : result_buffer_lines(run_config.storage_form, nrows))
    {
        lines.push_back({1, 0, 1, label});
    }
    lines.push_back({0, 0, 1, ""});

    // Obara-Saika PA/PB distances seed the vertical recurrence: PB when it builds
    // up the ket ((0|o|b) seeds), PA when it builds up the bra ((b|o|0) seeds).
    // A pure (0|o|0) target needs neither.
    SI2CIntegrals vrr_ints = vrr_base_ints;

    vrr_ints.insert(vrr_rest_ints.begin(), vrr_rest_ints.end());

    if (const auto side = vrr_build_side(vrr_ints); side == VrrBuildSide::ket)
    {
        lines.push_back({1, 0, 1, "// Obara-Saika PB distances (vertical recurrence builds the ket)"});
        lines.push_back({1, 0, 1, "const auto pb = osfunc::compute_pb(pair);"});
        lines.push_back({0, 0, 1, ""});
    }
    else if (side == VrrBuildSide::bra)
    {
        lines.push_back({1, 0, 1, "// Obara-Saika PA distances (vertical recurrence builds the bra)"});
        lines.push_back({1, 0, 1, "const auto pa = osfunc::compute_pa(pair);"});
        lines.push_back({0, 0, 1, ""});
    }

    // elementary primitive overlaps seed the vertical recurrence for overlap and
    // kinetic-energy integrals.
    if (seeds_with_primitive_overlap(run_config.operator_type))
    {
        if (is_primitive_overlap_seed(integral))
        {
            // (s|s) overlap: the contracted primitive overlaps are already the
            // spherical result (one component per side), so contract them straight
            // into the buffer and skip the Cartesian-to-spherical transform. The
            // seed is folded into the contract since nothing else consumes it here.
            lines.push_back({1, 0, 1, "// (s|s): contract the primitive overlaps into the buffer"});
            lines.push_back({1, 0, 1, "osfunc::contract(buffer, osfunc::compute_overlap(pair));"});
            lines.push_back({0, 0, 1, ""});

            lines.push_back({1, 0, 1, "return buffer;"});
            lines.push_back({0, 0, 2, "}"});
            lines.push_back({0, 0, 1, "}  // namespace " + tags.ns});

            ost::write_code_lines(fstream, lines);

            fstream.close();

            return;
        }

        // every other target consumes the seed in its VRR steps, so name it.
        lines.push_back({1, 0, 1, "// elementary primitive overlaps seed the vertical recurrence"});
        lines.push_back({1, 0, 1, "const auto povl_ss = osfunc::compute_overlap(pair);"});
        lines.push_back({0, 0, 1, ""});

        // each remaining VRR integral is produced by its own single-step VRR kernel
        // (os2c::vrr::<abbrev>::compute_<L>, L the built-side angular momentum);
        // those kernels are placeholders added separately. Emit lowest angular
        // momentum first so the lower integrals a step consumes precede it.
        std::vector<I2CIntegral> ordered(vrr_ints.begin(), vrr_ints.end());

        std::sort(ordered.begin(), ordered.end(), [](const I2CIntegral& a, const I2CIntegral& b) {
            if ((a[0] + a[1]) != (b[0] + b[1])) return (a[0] + a[1]) < (b[0] + b[1]);

            return a < b;
        });

        for (const auto& tint : ordered)
        {
            if (is_primitive_overlap_seed(tint)) continue;

            const auto abbr = integrand_abbrev(tint.integrand());

            // exactly one side is zero in a VRR integral; the other carries the
            // built angular momentum L. The kernel is labelled by L alone.
            const auto lval = tint[0] + tint[1];

            const auto ket_built = (tint[0] == 0);

            const auto var = "p" + abbr + "_" + shell_label(tint[0]) + shell_label(tint[1]);

            // a single VRR step for built momentum L takes the pair, the PA/PB
            // distances of the side it builds, and the two lower integrals it
            // recurs from (built momentum L-1 and L-2, highest first).
            std::string args = "pair";

            if (lval >= 1)
            {
                args += ket_built ? ", pb" : ", pa";

                for (int m = lval - 1; (m >= 0) && (m >= lval - 2); m--)
                {
                    const auto lower = ket_built ? ("s" + shell_label(m)) : (shell_label(m) + "s");

                    args += ", p" + abbr + "_" + lower;
                }
            }

            lines.push_back({1, 0, 1, "const auto " + var + " = os2c::vrr::" + abbr + "::compute_" +
                                          shell_label(lval) + "(" + args + ");"});
        }

        lines.push_back({0, 0, 1, ""});

        // contract the VRR base integrals (those HRR consumes) from primitive to
        // contracted form via osfunc::contract: each contracted block is sized to
        // its Cartesian component count and accumulated into from the primitives.
        std::vector<I2CIntegral> base_ordered(vrr_base_ints.begin(), vrr_base_ints.end());

        std::sort(base_ordered.begin(), base_ordered.end(),
                  [](const I2CIntegral& a, const I2CIntegral& b) {
                      if ((a[0] + a[1]) != (b[0] + b[1])) return (a[0] + a[1]) < (b[0] + b[1]);

                      return a < b;
                  });

        lines.push_back({1, 0, 1, "// contract the VRR base integrals consumed by HRR"});

        for (const auto& tint : base_ordered)
        {
            const auto abbr = integrand_abbrev(tint.integrand());

            const auto shells = shell_label(tint[0]) + shell_label(tint[1]);

            const auto rows = cartesian_count(tint[0]) * cartesian_count(tint[1]);

            lines.push_back({1, 0, 1, "auto c" + abbr + "_" + shells + " = CDenseMatrix(" +
                                          std::to_string(rows) + ", static_cast<int>(npairs));"});

            lines.push_back({1, 0, 1, "osfunc::contract(c" + abbr + "_" + shells + ", p" + abbr +
                                          "_" + shells + ");"});
        }

        lines.push_back({0, 0, 1, ""});

        // horizontal recurrence: when both sides carry angular momentum, transfer
        // momentum from one center to the other to reach the target, consuming the
        // contracted base integrals. HRR is operator-independent (a pure A-B
        // transfer), so its namespace carries no operator abbreviation.
        if (!hrr_ints.empty())
        {
            const auto abbr = integrand_abbrev(integral.integrand());

            const auto shells = shell_label(integral[0]) + shell_label(integral[1]);

            std::string args = "pair";

            for (const auto& tint : base_ordered)
            {
                args += ", c" + integrand_abbrev(tint.integrand()) + "_" + shell_label(tint[0]) +
                        shell_label(tint[1]);
            }

            lines.push_back({1, 0, 1, "// horizontal recurrence to the target integral"});

            lines.push_back({1, 0, 1, "const auto c" + abbr + "_" + shells + " = os2c::hrr::compute_" +
                                          shell_label(integral[0]) + "_" + shell_label(integral[1]) +
                                          "(" + args + ");"});

            lines.push_back({0, 0, 1, ""});
        }

        // transform the contracted Cartesian target to the spherical (real solid
        // harmonic) result, accumulating it into the return buffer.
        const auto target = "c" + integrand_abbrev(integral.integrand()) + "_" +
                            shell_label(integral[0]) + shell_label(integral[1]);

        lines.push_back({1, 0, 1, "// transform the Cartesian target to spherical and distribute"});

        lines.push_back({1, 0, 1, "osfunc::transform<" + std::to_string(integral[0]) + ", " +
                                      std::to_string(integral[1]) + ">(buffer, " + target + ");"});

        lines.push_back({0, 0, 1, ""});
    }

    lines.push_back({1, 0, 1, "return buffer;"});

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
