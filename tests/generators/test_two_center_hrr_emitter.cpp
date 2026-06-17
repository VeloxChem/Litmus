// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "run_configuration.hpp"
#include "two_center_hrr_emitter.hpp"
#include "two_center_hrr_generators.hpp"

namespace {

/// True if haystack contains needle.
bool
contains(const std::string& haystack, const std::string& needle)
{
    return haystack.find(needle) != std::string::npos;
}

/// Runs TwoCenterHrrGenerator::generate inside a private temporary directory (the
/// generator writes relative to the working directory) and returns that directory.
std::filesystem::path
generate_in_temp_dir(const cfg::RunConfiguration& run_config, const std::string& tag)
{
    const auto dir = std::filesystem::path(testing::TempDir()) / ("litmus_hrr_" + tag);

    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);

    const auto cwd = std::filesystem::current_path();
    std::filesystem::current_path(dir);

    TwoCenterHrrGenerator().generate(run_config);

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
hrr_config(cfg::RecursionType type, int min_ang_mom, int max_ang_mom)
{
    cfg::RunConfiguration run_config;
    run_config.recursion_type = type;
    run_config.min_ang_mom    = min_ang_mom;
    run_config.max_ang_mom    = max_ang_mom;
    return run_config;
}

}  // namespace

TEST(TwoCenterHrrEmitterTest, SignatureBraIncrement)
{
    // (p|p): bra-incremented, full reduction consumes (s|p) and (s|d).
    const auto sig = format_hrr_signature(1, 1);

    EXPECT_TRUE(contains(sig, "void compute_p_p("));
    EXPECT_TRUE(contains(sig, "const osfunc::CArray<double>& sp,"));
    EXPECT_TRUE(contains(sig, "const osfunc::CArray<double>& sd,"));
    EXPECT_TRUE(contains(sig, "const osfunc::CArray<double>& ab,"));
    EXPECT_TRUE(contains(sig, "osfunc::CArray<double>& pp)"));
}

TEST(TwoCenterHrrEmitterTest, SignatureFullReductionToS)
{
    // (d|d) full recursion reduces the bra to s, so the bases are (s|d), (s|f),
    // (s|g) -- not the single-step (p|d), (p|f).
    const auto sig = format_hrr_signature(2, 2);

    EXPECT_TRUE(contains(sig, "void compute_d_d("));
    EXPECT_TRUE(contains(sig, "& sd,"));
    EXPECT_TRUE(contains(sig, "& sf,"));
    EXPECT_TRUE(contains(sig, "& sg,"));
    EXPECT_TRUE(contains(sig, "& dd)"));
    EXPECT_FALSE(contains(sig, "& pd,"));
    EXPECT_FALSE(contains(sig, "& pf,"));
}

TEST(TwoCenterHrrEmitterTest, SignatureKetIncrement)
{
    // (d|p): ket-incremented (la > lb), bases (d|s) and (f|s).
    const auto sig = format_hrr_signature(2, 1);

    EXPECT_TRUE(contains(sig, "void compute_d_p("));
    EXPECT_TRUE(contains(sig, "& ds,"));
    EXPECT_TRUE(contains(sig, "& fs,"));
    EXPECT_TRUE(contains(sig, "osfunc::CArray<double>& dp)"));
}

TEST(TwoCenterHrrEmitterTest, KernelScaffolding)
{
    const auto src = format_hrr_kernel(1, 1);

    EXPECT_TRUE(contains(src, "const std::size_t ncomps = 9;"));            // (2*1+1)^2
    EXPECT_TRUE(contains(src, "const auto nblocks = pp.nrows() / ncomps;"));
    EXPECT_TRUE(contains(src, "const auto npairs = pp.ncols();"));
    EXPECT_TRUE(contains(src, "for (std::size_t iblock = 0; iblock < nblocks; iblock++)"));
    EXPECT_TRUE(contains(src, "#pragma omp simd aligned("));
    EXPECT_TRUE(contains(src, "for (std::size_t i = 0; i < npairs; i++)"));
}

TEST(TwoCenterHrrEmitterTest, KernelPerIntegralOffsets)
{
    // each integral is offset by iblock times its own component count.
    const auto src = format_hrr_kernel(1, 1);

    EXPECT_TRUE(contains(src, "const auto sp_off = iblock * 3;"));   // (s|p): 3 Cartesian
    EXPECT_TRUE(contains(src, "const auto sd_off = iblock * 6;"));   // (s|d): 6 Cartesian
    EXPECT_TRUE(contains(src, "const auto pp_off = iblock * 9;"));   // (p|p): 9 spherical
}

TEST(TwoCenterHrrEmitterTest, KernelPPRecurrenceValues)
{
    // bra-side HRR (p_i|p_j) = (s|d_{i+j}) - AB_i (s|p_j), with the p-shell
    // spherical ordering m=-1->y, 0->z, +1->x. Terms are grouped by AB product
    // (the AB-degree term leads), one summand per line.
    const auto src = format_hrr_kernel(1, 1);

    EXPECT_TRUE(contains(src, "pp_0[i] = -sp_1[i] * ab_y[i]\n                    + sd_3[i];"));  // (p_y|p_y)
    EXPECT_TRUE(contains(src, "pp_8[i] = -sp_0[i] * ab_x[i]\n                    + sd_0[i];"));  // (p_x|p_x)
}

TEST(TwoCenterHrrEmitterTest, KernelHoistsRadicalConstant)
{
    // the d-shell ket transform introduces sqrt(3), hoisted once outside the loops.
    const auto src = format_hrr_kernel(1, 2);

    EXPECT_TRUE(contains(src, "const double f3 = std::sqrt(3.0);"));
    EXPECT_TRUE(contains(src, "* f3 *"));
}

TEST(TwoCenterHrrEmitterTest, KernelFullReductionConsumesSShells)
{
    // (d|d) reduces all the way to s: the body reads sd/sf/sg, never pd/pf.
    const auto src = format_hrr_kernel(2, 2);

    EXPECT_TRUE(contains(src, "const std::size_t ncomps = 25;"));
    EXPECT_TRUE(contains(src, "sd.row("));
    EXPECT_TRUE(contains(src, "sf.row("));
    EXPECT_TRUE(contains(src, "sg.row("));
    EXPECT_FALSE(contains(src, "pd.row("));
    EXPECT_FALSE(contains(src, "pf.row("));
}

TEST(TwoCenterHrrEmitterTest, KernelGroupsAndFactorsCoefficients)
{
    // terms sharing an AB product are grouped and the common rational factored out.
    const auto src = format_hrr_kernel(2, 2);

    EXPECT_TRUE(contains(
        src, "0.25 * (sd_0[i] + sd_3[i] - 2.0 * sd_5[i]) * ab_x[i] * ab_x[i]"));
}

TEST(TwoCenterHrrEmitterTest, KetIncrementGroupsByKet)
{
    // for la > lb the single-component loops are grouped by the ket (the
    // incremented side), so each loop pins one AB axis.
    const auto src = format_hrr_kernel(2, 1);

    EXPECT_TRUE(contains(src, "// ket spherical component 0"));
    EXPECT_FALSE(contains(src, "// bra spherical component 0"));
}

TEST(TwoCenterHrrGeneratorTest, WritesSelectedKernelPairs)
{
    const auto dir = generate_in_temp_dir(hrr_config(cfg::RecursionType::hrr_bra_ket, 1, 2), "all");

    // both sides >= 1, so the s-involving (SP, PS, ...) targets are skipped.
    for (const std::string base : {"ObaraSaikaTwoCenterHrrPP",
                                   "ObaraSaikaTwoCenterHrrPD",
                                   "ObaraSaikaTwoCenterHrrDP",
                                   "ObaraSaikaTwoCenterHrrDD"})
    {
        EXPECT_TRUE(std::filesystem::exists(dir / (base + ".hpp"))) << base;
        EXPECT_TRUE(std::filesystem::exists(dir / (base + ".cpp"))) << base;
    }

    // the header guards the declaration in the os2c::hrr namespace.
    const auto hpp = read_file(dir / "ObaraSaikaTwoCenterHrrPP.hpp");
    EXPECT_TRUE(contains(hpp, "#ifndef ObaraSaikaTwoCenterHrrPP_hpp"));
    EXPECT_TRUE(contains(hpp, "#include \"Array.hpp\""));
    EXPECT_TRUE(contains(hpp, "namespace os2c::hrr {"));
    EXPECT_TRUE(contains(hpp, "void compute_p_p("));

    // the source includes its header and <cmath> and carries the definition.
    const auto cpp = read_file(dir / "ObaraSaikaTwoCenterHrrPP.cpp");
    EXPECT_TRUE(contains(cpp, "#include \"ObaraSaikaTwoCenterHrrPP.hpp\""));
    EXPECT_TRUE(contains(cpp, "#include <cmath>"));
    EXPECT_TRUE(contains(cpp, "#pragma omp simd aligned("));
}

TEST(TwoCenterHrrGeneratorTest, RecursionTypeSelectsTransferSide)
{
    // hrr_bra keeps only la <= lb; hrr_ket only la > lb.
    const auto bra = generate_in_temp_dir(hrr_config(cfg::RecursionType::hrr_bra, 1, 2), "bra");
    EXPECT_TRUE(std::filesystem::exists(bra / "ObaraSaikaTwoCenterHrrPD.hpp"));
    EXPECT_FALSE(std::filesystem::exists(bra / "ObaraSaikaTwoCenterHrrDP.hpp"));

    const auto ket = generate_in_temp_dir(hrr_config(cfg::RecursionType::hrr_ket, 1, 2), "ket");
    EXPECT_TRUE(std::filesystem::exists(ket / "ObaraSaikaTwoCenterHrrDP.hpp"));
    EXPECT_FALSE(std::filesystem::exists(ket / "ObaraSaikaTwoCenterHrrPD.hpp"));
}
