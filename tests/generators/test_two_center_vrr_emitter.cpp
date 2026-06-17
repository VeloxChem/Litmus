// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "run_configuration.hpp"
#include "two_center_vrr_emitter.hpp"
#include "two_center_vrr_generators.hpp"

namespace {

/// True if haystack contains needle.
bool
contains(const std::string& haystack, const std::string& needle)
{
    return haystack.find(needle) != std::string::npos;
}

/// Runs TwoCenterVrrGenerator::generate inside a private temporary directory (the
/// generator writes relative to the working directory) and returns that directory.
std::filesystem::path
generate_in_temp_dir(const cfg::RunConfiguration& run_config, const std::string& tag)
{
    const auto dir = std::filesystem::path(testing::TempDir()) / ("litmus_vrr_" + tag);

    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);

    const auto cwd = std::filesystem::current_path();
    std::filesystem::current_path(dir);

    TwoCenterVrrGenerator().generate(run_config);

    std::filesystem::current_path(cwd);

    return dir;
}

std::string
read_file(const std::filesystem::path& path)
{
    std::ifstream in(path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

cfg::RunConfiguration
vrr_config(cfg::RecursionType type, int min_ang_mom, int max_ang_mom)
{
    cfg::RunConfiguration run_config;
    run_config.recursion_type = type;
    run_config.min_ang_mom    = min_ang_mom;
    run_config.max_ang_mom    = max_ang_mom;
    return run_config;
}

}  // namespace

TEST(TwoCenterVrrEmitterTest, CartesianSignatureLowerShells)
{
    // (s|d) single step consumes the lower (s|p) and (s|s) Cartesian integrals.
    const auto sig = format_vrr_cartesian_signature(2);

    EXPECT_TRUE(contains(sig, "void compute_d("));
    EXPECT_TRUE(contains(sig, "const osfunc::CBasisFunctionPair& pair,"));
    EXPECT_TRUE(contains(sig, "const osfunc::CArray<double>& sp,"));
    EXPECT_TRUE(contains(sig, "const osfunc::CArray<double>& ss,"));
    EXPECT_TRUE(contains(sig, "const osfunc::CArray<double>& pc,"));
    EXPECT_TRUE(contains(sig, "osfunc::CArray<double>& sd)"));
}

TEST(TwoCenterVrrEmitterTest, CartesianSignaturePHasOnlyOneLowerShell)
{
    // (s|p) has no (s|lb-2) term, so it takes only the (s|s) seed.
    const auto sig = format_vrr_cartesian_signature(1);

    EXPECT_TRUE(contains(sig, "void compute_p("));
    EXPECT_TRUE(contains(sig, "const osfunc::CArray<double>& ss,"));
    EXPECT_TRUE(contains(sig, "osfunc::CArray<double>& sp)"));
}

TEST(TwoCenterVrrEmitterTest, SphericalSignatureSeedOnly)
{
    // the full spherical (s|d) builds from the (s|s) seed alone.
    const auto sig = format_vrr_spherical_signature(2);

    EXPECT_TRUE(contains(sig, "void compute_d_sph("));
    EXPECT_TRUE(contains(sig, "const osfunc::CArray<double>& ss,"));
    EXPECT_TRUE(contains(sig, "const osfunc::CArray<double>& pc,"));
    EXPECT_TRUE(contains(sig, "osfunc::CArray<double>& sd)"));
    EXPECT_FALSE(contains(sig, "& sp,"));  // no lower (s|p) input
}

TEST(TwoCenterVrrEmitterTest, CartesianKernelRetainsFe)
{
    // single-step Cartesian VRR is per-primitive and keeps the fe = 1/(2 eta) term.
    const auto src = format_vrr_cartesian_kernel(2);

    EXPECT_TRUE(contains(src, "const auto fe = 1.0 / (2.0 * (bra_exps[p] + ket_exps[q]));"));
    EXPECT_TRUE(contains(src, "const auto ip = p * ket.number_of_primitive_functions() + q;"));

    // per-integral primitive offsets (component count of each integral).
    EXPECT_TRUE(contains(src, "sd.row(ip * 6 + 0)"));   // (s|d): 6 Cartesian
    EXPECT_TRUE(contains(src, "sp.row(ip * 3 + 0)"));   // (s|p): 3 Cartesian
    EXPECT_TRUE(contains(src, "ss.row(ip * 1 + 0)"));   // (s|s): 1

    // the recurrence: diagonal carries the fe lowering term, off-diagonal does not.
    EXPECT_TRUE(contains(src, "sd_0[i] = pc_x[i] * sp_0[i] + fe * ss_0[i];"));  // xx
    EXPECT_TRUE(contains(src, "sd_1[i] = pc_x[i] * sp_1[i];"));                 // xy
}

TEST(TwoCenterVrrEmitterTest, CartesianKernelCountsLoweringMultiplicity)
{
    // (s|f) lowers a d-shell exponent of 2, so the fe term carries the count 2.
    const auto src = format_vrr_cartesian_kernel(3);

    EXPECT_TRUE(contains(src, "sf_0[i] = pc_x[i] * sd_0[i] + 2.0 * fe * sp_0[i];"));  // xxx
}

TEST(TwoCenterVrrEmitterTest, SphericalKernelCancelsFe)
{
    // the traceless spherical transform cancels every fe term, so the spherical
    // kernel never computes fe and never touches the exponents.
    const auto src = format_vrr_spherical_kernel(2);

    EXPECT_FALSE(contains(src, "fe"));
    EXPECT_FALSE(contains(src, "bra_exps"));

    // it accumulates the pure pc polynomial times the seed into the result.
    EXPECT_TRUE(contains(src, "const double f3 = std::sqrt(3.0);"));
    EXPECT_TRUE(contains(src, "sd_0[i] += f3 * pc_x[i] * pc_y[i] * t_ss[i];"));

    // the m=0 component factors out -0.5; the long inner sum wraps one term per line.
    EXPECT_TRUE(contains(src, "sd_2[i] += -0.5 * (pc_x[i] * pc_x[i]"));
    EXPECT_TRUE(contains(src, "- 2.0 * pc_z[i] * pc_z[i]) * t_ss[i];"));
}

TEST(TwoCenterVrrEmitterTest, SphericalKernelGcdFactoring)
{
    // the (s|f) spherical components factor the common rational GCD (1/4) and the
    // radical (f10), leaving integer inner coefficients.
    const auto src = format_vrr_spherical_kernel(3);

    EXPECT_TRUE(contains(src, "sf_0[i] += 0.25 * f10 * (3.0 * pc_x[i] * pc_x[i] * pc_y[i]"));
    EXPECT_TRUE(contains(src, "- pc_y[i] * pc_y[i] * pc_y[i]) * t_ss[i];"));
}

TEST(TwoCenterVrrGeneratorTest, CartesianWritesKernelPairs)
{
    const auto dir = generate_in_temp_dir(vrr_config(cfg::RecursionType::vrr_cartesian, 1, 3), "cart");

    for (const std::string base : {"ObaraSaikaTwoCenterOverlapVrrCartP",
                                   "ObaraSaikaTwoCenterOverlapVrrCartD",
                                   "ObaraSaikaTwoCenterOverlapVrrCartF"})
    {
        EXPECT_TRUE(std::filesystem::exists(dir / (base + ".hpp"))) << base;
        EXPECT_TRUE(std::filesystem::exists(dir / (base + ".cpp"))) << base;
    }

    // the Cartesian VRR lives in the os2c::vrr::ovl namespace.
    const auto hpp = read_file(dir / "ObaraSaikaTwoCenterOverlapVrrCartD.hpp");
    EXPECT_TRUE(contains(hpp, "#include \"BasisFunctionPair.hpp\""));
    EXPECT_TRUE(contains(hpp, "namespace os2c::vrr::ovl {"));
    EXPECT_TRUE(contains(hpp, "void compute_d("));

    // the source carries the definition; it has no transcendental factors.
    const auto cpp = read_file(dir / "ObaraSaikaTwoCenterOverlapVrrCartD.cpp");
    EXPECT_TRUE(contains(cpp, "#include \"ObaraSaikaTwoCenterOverlapVrrCartD.hpp\""));
    EXPECT_FALSE(contains(cpp, "#include <cmath>"));
    EXPECT_TRUE(contains(cpp, "#pragma omp simd aligned("));
}

TEST(TwoCenterVrrGeneratorTest, SphericalWritesKernelPairs)
{
    const auto dir = generate_in_temp_dir(vrr_config(cfg::RecursionType::vrr_spherical, 1, 3), "sph");

    EXPECT_TRUE(std::filesystem::exists(dir / "ObaraSaikaTwoCenterOverlapVrrSphD.hpp"));
    EXPECT_TRUE(std::filesystem::exists(dir / "ObaraSaikaTwoCenterOverlapVrrSphD.cpp"));

    // the spherical VRR lives in the os2c::ovl namespace and needs <cmath>.
    const auto hpp = read_file(dir / "ObaraSaikaTwoCenterOverlapVrrSphD.hpp");
    EXPECT_TRUE(contains(hpp, "namespace os2c::ovl {"));
    EXPECT_TRUE(contains(hpp, "void compute_d_sph("));

    const auto cpp = read_file(dir / "ObaraSaikaTwoCenterOverlapVrrSphD.cpp");
    EXPECT_TRUE(contains(cpp, "#include <cmath>"));
}
