// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>

#include "run_configuration.hpp"
#include "two_center_generators.hpp"

namespace {

// Runs TwoCenterGenerator::generate inside a private temporary directory (the
// emitter writes its files relative to the working directory) and returns the
// path that directory, so a test can read the emitted files back.
std::filesystem::path
generate_in_temp_dir(const cfg::RunConfiguration& run_config, const std::string& tag)
{
    const auto dir = std::filesystem::path(testing::TempDir()) / ("litmus_emit_" + tag);

    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);

    const auto cwd = std::filesystem::current_path();
    std::filesystem::current_path(dir);

    TwoCenterGenerator().generate(run_config);

    std::filesystem::current_path(cwd);

    return dir;
}

// Reads a whole file into a string.
std::string
read_file(const std::filesystem::path& path)
{
    std::ifstream in(path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

cfg::RunConfiguration
overlap_config(int max_ang_mom)
{
    cfg::RunConfiguration run_config;
    run_config.integral_type = cfg::IntegralType::two_center;
    run_config.operator_type = cfg::OperatorType::overlap;
    run_config.max_ang_mom   = max_ang_mom;
    return run_config;
}

}  // namespace

TEST(TwoCenterEmittersTest, EmitsHeaderAndSourcePairPerTarget)
{
    const auto dir = generate_in_temp_dir(overlap_config(1), "pair");

    // bra/ket each span [0, 1] -> SS, SP, PS, PP, each an .hpp/.cpp pair.
    for (const auto& base : {"ObaraSaikaTwoCenterOverlapSS",
                             "ObaraSaikaTwoCenterOverlapSP",
                             "ObaraSaikaTwoCenterOverlapPS",
                             "ObaraSaikaTwoCenterOverlapPP"})
    {
        EXPECT_TRUE(std::filesystem::exists(dir / (std::string(base) + ".hpp"))) << base;
        EXPECT_TRUE(std::filesystem::exists(dir / (std::string(base) + ".cpp"))) << base;
    }
}

TEST(TwoCenterEmittersTest, HeaderDeclaresGuardedNamespacedKernel)
{
    const auto dir = generate_in_temp_dir(overlap_config(0), "header");

    const auto hpp = read_file(dir / "ObaraSaikaTwoCenterOverlapSS.hpp");

    EXPECT_NE(hpp.find("#ifndef ObaraSaikaTwoCenterOverlapSS_hpp"), std::string::npos);
    EXPECT_NE(hpp.find("#endif /* ObaraSaikaTwoCenterOverlapSS_hpp */"), std::string::npos);
    EXPECT_NE(hpp.find("namespace os2c::ovl"), std::string::npos);
    EXPECT_NE(hpp.find("auto\ncompute_s_s("), std::string::npos);

    // input from the signature, return type from the storage form
    EXPECT_NE(hpp.find("const osfunc::CBasisFunctionPair& pair"), std::string::npos);
    EXPECT_NE(hpp.find("-> osfunc::CArray<double>;"), std::string::npos);

    // the kernel pulls in the osfunc data-structure headers of its input/return
    EXPECT_NE(hpp.find("#include \"Array.hpp\""), std::string::npos);
    EXPECT_NE(hpp.find("#include \"BasisFunctionPair.hpp\""), std::string::npos);
}

TEST(TwoCenterEmittersTest, SourceDefinesBodyWithBufferAndReturn)
{
    const auto dir = generate_in_temp_dir(overlap_config(0), "body");

    const auto cpp = read_file(dir / "ObaraSaikaTwoCenterOverlapSS.cpp");

    // includes its own header, the Obara-Saika helpers, and reopens the namespace
    EXPECT_NE(cpp.find("#include \"ObaraSaikaTwoCenterOverlapSS.hpp\""), std::string::npos);
    EXPECT_NE(cpp.find("#include \"ObaraSaikaFunc.hpp\""), std::string::npos);
    EXPECT_NE(cpp.find("auto\ncompute_s_s("), std::string::npos);
    EXPECT_NE(cpp.find("-> osfunc::CArray<double>\n{"), std::string::npos);

    // signature-aware input and storage-form-aware return buffer
    EXPECT_NE(cpp.find("const osfunc::CBasisFunctionPair& pair"), std::string::npos);
    EXPECT_NE(cpp.find("pair.number_of_pairs()"), std::string::npos);
    // (S|1|S): (2*0+1)*(2*0+1) = 1 row, one column per pair, zero-initialized
    EXPECT_NE(cpp.find("osfunc::CArray<double> buffer(1, npairs);"), std::string::npos);
    EXPECT_EQ(cpp.find("buffer.zero();"), std::string::npos);  // ctor already zeroes
    EXPECT_NE(cpp.find("return buffer;"), std::string::npos);
}

TEST(TwoCenterEmittersTest, VrrSeedsSelectPaPbOrNeither)
{
    const auto dir = generate_in_temp_dir(overlap_config(1), "papb");

    // (S|1|S): la == lb == 0, neither PA nor PB is computed.
    const auto ss = read_file(dir / "ObaraSaikaTwoCenterOverlapSS.cpp");
    EXPECT_EQ(ss.find("osfunc::compute_pa(pair)"), std::string::npos);
    EXPECT_EQ(ss.find("osfunc::compute_pb(pair)"), std::string::npos);

    // (S|1|P): VRR has (0|o|b) seeds -> the ket is built -> PB.
    const auto sp = read_file(dir / "ObaraSaikaTwoCenterOverlapSP.cpp");
    EXPECT_NE(sp.find("const auto pb = osfunc::compute_pb(pair);"), std::string::npos);
    EXPECT_EQ(sp.find("osfunc::compute_pa(pair)"), std::string::npos);

    // (P|1|S): VRR has (b|o|0) seeds -> the bra is built -> PA.
    const auto ps = read_file(dir / "ObaraSaikaTwoCenterOverlapPS.cpp");
    EXPECT_NE(ps.find("const auto pa = osfunc::compute_pa(pair);"), std::string::npos);
    EXPECT_EQ(ps.find("osfunc::compute_pb(pair)"), std::string::npos);

    // (P|1|P): a <= b keeps the bra at 0 and builds the ket -> PB.
    const auto pp = read_file(dir / "ObaraSaikaTwoCenterOverlapPP.cpp");
    EXPECT_NE(pp.find("const auto pb = osfunc::compute_pb(pair);"), std::string::npos);
    EXPECT_EQ(pp.find("osfunc::compute_pa(pair)"), std::string::npos);
}

cfg::RunConfiguration
config_for(cfg::OperatorType op, int max_ang_mom)
{
    cfg::RunConfiguration run_config;
    run_config.integral_type = cfg::IntegralType::two_center;
    run_config.operator_type = op;
    run_config.max_ang_mom   = max_ang_mom;
    return run_config;
}

TEST(TwoCenterEmittersTest, OverlapEmitsWorkflowOtherOperatorsStub)
{
    // overlap (s|s) contracts the primitive overlaps straight into the result.
    const auto ovl_dir = generate_in_temp_dir(config_for(cfg::OperatorType::overlap, 0), "ovl_seed");
    const auto ovl = read_file(ovl_dir / "ObaraSaikaTwoCenterOverlapSS.cpp");
    EXPECT_NE(ovl.find("osfunc::contract(buffer, osfunc::compute_overlap(pair));"), std::string::npos);

    // only the overlap kernels are generated so far; other operators return the
    // zeroed buffer with a TODO and call no kernels.
    const auto kin_dir =
        generate_in_temp_dir(config_for(cfg::OperatorType::kinetic_energy, 0), "kin_seed");
    const auto kin = read_file(kin_dir / "ObaraSaikaTwoCenterKineticEnergySS.cpp");
    EXPECT_NE(kin.find("not generated yet"), std::string::npos);
    EXPECT_EQ(kin.find("osfunc::compute_overlap"), std::string::npos);
    EXPECT_NE(kin.find("return buffer;"), std::string::npos);

    const auto eri_dir =
        generate_in_temp_dir(config_for(cfg::OperatorType::electron_repulsion, 0), "eri_seed");
    const auto eri = read_file(eri_dir / "ObaraSaikaTwoCenterElectronRepulsionSS.cpp");
    EXPECT_NE(eri.find("not generated yet"), std::string::npos);
    EXPECT_EQ(eri.find("osfunc::compute_overlap"), std::string::npos);
}

TEST(TwoCenterEmittersTest, HrrTargetsBuildCartesianVrrLadderThenContract)
{
    const auto dir = generate_in_temp_dir(overlap_config(1), "vrr");

    // (P|1|P): both sides carry momentum, so the ket is built by the Cartesian VRR
    // (os2c::vrr::ovl::compute_<L>(pair, <lower>, pc, <out>)), lowest first, then
    // contracted into Cartesian blocks sized ncart(la)*ncart(lb).
    const auto pp = read_file(dir / "ObaraSaikaTwoCenterOverlapPP.cpp");
    EXPECT_NE(pp.find("osfunc::CArray<double> sp(nprims * 3, npairs);"), std::string::npos);
    EXPECT_NE(pp.find("os2c::vrr::ovl::compute_p(pair, ss, pb, sp);"), std::string::npos);
    EXPECT_NE(pp.find("os2c::vrr::ovl::compute_d(pair, sp, ss, pb, sd);"), std::string::npos);
    EXPECT_LT(pp.find("compute_p(pair, ss, pb, sp)"), pp.find("compute_d(pair, sp, ss, pb, sd)"));

    EXPECT_NE(pp.find("osfunc::CArray<double> csp(3, npairs);"), std::string::npos);
    EXPECT_NE(pp.find("osfunc::contract(csp, sp);"), std::string::npos);
    EXPECT_NE(pp.find("osfunc::CArray<double> csd(6, npairs);"), std::string::npos);
    EXPECT_NE(pp.find("osfunc::contract(csd, sd);"), std::string::npos);
}

TEST(TwoCenterEmittersTest, NoHrrTargetsCallTheFusedSphericalVrr)
{
    const auto dir = generate_in_temp_dir(overlap_config(1), "sph");

    // (S|1|P): one side is zero, so the spherical VRR builds the ket with PB and
    // fuses contraction + transform; no Cartesian VRR or HRR is needed.
    const auto sp = read_file(dir / "ObaraSaikaTwoCenterOverlapSP.cpp");
    EXPECT_NE(sp.find("os2c::ovl::compute_p_sph(pair, ss, pb, buffer);"), std::string::npos);
    EXPECT_EQ(sp.find("os2c::vrr::"), std::string::npos);
    EXPECT_EQ(sp.find("os2c::hrr::"), std::string::npos);

    // (P|1|S): the bra is built instead, with the PA distances.
    const auto ps = read_file(dir / "ObaraSaikaTwoCenterOverlapPS.cpp");
    EXPECT_NE(ps.find("os2c::ovl::compute_p_sph(pair, ss, pa, buffer);"), std::string::npos);

    // (S|1|S): only the contracted seed, no recurrence kernels at all.
    const auto ss = read_file(dir / "ObaraSaikaTwoCenterOverlapSS.cpp");
    EXPECT_EQ(ss.find("os2c::vrr::"), std::string::npos);
    EXPECT_EQ(ss.find("os2c::ovl::compute_"), std::string::npos);
}

TEST(TwoCenterEmittersTest, HrrFusesTransformNoStandaloneTransform)
{
    const auto dir = generate_in_temp_dir(overlap_config(1), "hrr");

    // (P|1|P): the HRR consumes the contracted base integrals and writes the
    // spherical buffer directly (the transform is fused), so no osfunc::transform.
    const auto pp = read_file(dir / "ObaraSaikaTwoCenterOverlapPP.cpp");
    EXPECT_NE(pp.find("const auto ab = osfunc::compute_ab(pair);"), std::string::npos);
    EXPECT_NE(pp.find("os2c::hrr::compute_p_p(csp, csd, ab, buffer);"), std::string::npos);
    EXPECT_LT(pp.find("os2c::hrr::compute_p_p"), pp.find("return buffer;"));

    // the orchestration includes the kernel headers it calls.
    EXPECT_NE(pp.find("#include \"ObaraSaikaTwoCenterHrrPP.hpp\""), std::string::npos);
    EXPECT_NE(pp.find("#include \"ObaraSaikaTwoCenterOverlapVrrCartD.hpp\""), std::string::npos);

    // no standalone Cartesian-to-spherical transform anywhere.
    for (const auto* base : {"ObaraSaikaTwoCenterOverlapPP",
                             "ObaraSaikaTwoCenterOverlapSP",
                             "ObaraSaikaTwoCenterOverlapSS"})
    {
        const auto cpp = read_file(dir / (std::string(base) + ".cpp"));
        EXPECT_EQ(cpp.find("osfunc::transform"), std::string::npos) << base;
        EXPECT_EQ(cpp.find("CartesianToSphericalFunc.hpp"), std::string::npos) << base;
    }
}

TEST(TwoCenterEmittersTest, MakeEmitterRejectsNothingForSupportedCpuCpp)
{
    // the default configuration targets cpu/C++, for which an emitter exists.
    EXPECT_NO_THROW({ generate_in_temp_dir(overlap_config(0), "supported"); });
}
