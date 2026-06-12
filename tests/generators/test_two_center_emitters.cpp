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
    EXPECT_NE(hpp.find("const CScreenedBasisFunctionPair& pair"), std::string::npos);
    EXPECT_NE(hpp.find("-> CDenseMatrix;"), std::string::npos);

    // the kernel pulls in the headers of its input and return types
    EXPECT_NE(hpp.find("#include \"DenseMatrix.hpp\""), std::string::npos);
    EXPECT_NE(hpp.find("#include \"ScreenedBasisFunctionPair.hpp\""), std::string::npos);
}

TEST(TwoCenterEmittersTest, SourceDefinesBodyWithBufferAndReturn)
{
    const auto dir = generate_in_temp_dir(overlap_config(0), "body");

    const auto cpp = read_file(dir / "ObaraSaikaTwoCenterOverlapSS.cpp");

    // includes its own header, the Obara-Saika helpers, and reopens the namespace
    EXPECT_NE(cpp.find("#include \"ObaraSaikaTwoCenterOverlapSS.hpp\""), std::string::npos);
    EXPECT_NE(cpp.find("#include \"ObaraSaikaFunc.hpp\""), std::string::npos);
    EXPECT_NE(cpp.find("auto\ncompute_s_s("), std::string::npos);
    EXPECT_NE(cpp.find("-> CDenseMatrix\n{"), std::string::npos);

    // signature-aware input and storage-form-aware return buffer
    EXPECT_NE(cpp.find("const CScreenedBasisFunctionPair& pair"), std::string::npos);
    EXPECT_NE(cpp.find("pair.number_of_pairs()"), std::string::npos);
    // (S|1|S): (2*0+1)*(2*0+1) = 1 row, one column per pair, zero-initialized
    EXPECT_NE(cpp.find("auto buffer = CDenseMatrix(1, static_cast<int>(npairs));"),
              std::string::npos);
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

TEST(TwoCenterEmittersTest, PrimitiveOverlapSeedsOverlapAndKineticNotEri)
{
    // overlap seeds from the primitive overlaps, even at (S|1|S) (where the seed
    // is folded straight into the contract rather than named).
    const auto ovl_dir = generate_in_temp_dir(config_for(cfg::OperatorType::overlap, 0), "ovl_seed");
    const auto ovl = read_file(ovl_dir / "ObaraSaikaTwoCenterOverlapSS.cpp");
    EXPECT_NE(ovl.find("osfunc::compute_overlap(pair)"), std::string::npos);

    // kinetic energy also seeds from the primitive overlaps.
    const auto kin_dir =
        generate_in_temp_dir(config_for(cfg::OperatorType::kinetic_energy, 0), "kin_seed");
    const auto kin = read_file(kin_dir / "ObaraSaikaTwoCenterKineticEnergySS.cpp");
    EXPECT_NE(kin.find("osfunc::compute_overlap(pair)"), std::string::npos);

    // electron repulsion does not seed from the primitive overlaps.
    const auto eri_dir =
        generate_in_temp_dir(config_for(cfg::OperatorType::electron_repulsion, 0), "eri_seed");
    const auto eri = read_file(eri_dir / "ObaraSaikaTwoCenterElectronRepulsionSS.cpp");
    EXPECT_EQ(eri.find("osfunc::compute_overlap(pair)"), std::string::npos);
}

TEST(TwoCenterEmittersTest, VrrIntegralsCallNamespacedKernels)
{
    const auto dir = generate_in_temp_dir(overlap_config(1), "vrr");

    // (P|1|P): VRR = {(s|s), (s|p), (s|d)}. The (s|s) seed is povl_ss; the others
    // call os2c::vrr::ovl::compute_<L> with the PB distances and the lower
    // integrals the single VRR step consumes (L-1, L-2; highest first).
    const auto pp = read_file(dir / "ObaraSaikaTwoCenterOverlapPP.cpp");
    EXPECT_NE(pp.find("const auto povl_sp = os2c::vrr::ovl::compute_p(pair, pb, povl_ss);"),
              std::string::npos);
    EXPECT_NE(pp.find("const auto povl_sd = os2c::vrr::ovl::compute_d(pair, pb, povl_sp, povl_ss);"),
              std::string::npos);
    EXPECT_LT(pp.find("povl_sp ="), pp.find("povl_sd ="));  // dependency order

    // (P|1|S): VRR builds the bra -> compute_p with the PA distances.
    const auto ps = read_file(dir / "ObaraSaikaTwoCenterOverlapPS.cpp");
    EXPECT_NE(ps.find("const auto povl_ps = os2c::vrr::ovl::compute_p(pair, pa, povl_ss);"),
              std::string::npos);

    // (S|1|S): only the seed, no VRR kernel calls.
    const auto ss = read_file(dir / "ObaraSaikaTwoCenterOverlapSS.cpp");
    EXPECT_EQ(ss.find("os2c::vrr::"), std::string::npos);
}

TEST(TwoCenterEmittersTest, VrrBaseIntegralsAreContractedForHrr)
{
    const auto dir = generate_in_temp_dir(overlap_config(1), "contract");

    // (P|1|P): VRR base = {(s|p), (s|d)} -> contract into Cartesian blocks sized
    // ncart(la)*ncart(lb): (s|p) -> 1*3 = 3 rows, (s|d) -> 1*6 = 6 rows.
    const auto pp = read_file(dir / "ObaraSaikaTwoCenterOverlapPP.cpp");
    EXPECT_NE(pp.find("auto covl_sp = CDenseMatrix(3, static_cast<int>(npairs));"),
              std::string::npos);
    EXPECT_NE(pp.find("osfunc::contract(covl_sp, povl_sp);"), std::string::npos);
    EXPECT_NE(pp.find("auto covl_sd = CDenseMatrix(6, static_cast<int>(npairs));"),
              std::string::npos);
    EXPECT_NE(pp.find("osfunc::contract(covl_sd, povl_sd);"), std::string::npos);

    // (S|1|S) overlap: the primitive overlaps are contracted directly into the
    // buffer with the seed folded into the contract (no named povl_ss, no
    // intermediate covl_ss block, no transform, no transform header).
    const auto ss = read_file(dir / "ObaraSaikaTwoCenterOverlapSS.cpp");
    EXPECT_NE(ss.find("osfunc::contract(buffer, osfunc::compute_overlap(pair));"),
              std::string::npos);
    EXPECT_EQ(ss.find("povl_ss"), std::string::npos);
    EXPECT_EQ(ss.find("covl_ss"), std::string::npos);
    EXPECT_EQ(ss.find("osfunc::transform"), std::string::npos);
    EXPECT_EQ(ss.find("CartesianToSphericalFunc.hpp"), std::string::npos);
}

TEST(TwoCenterEmittersTest, HrrTransfersToTargetWhenBothSidesCarryMomentum)
{
    const auto dir = generate_in_temp_dir(overlap_config(1), "hrr");

    // (P|1|P): both sides carry momentum -> HRR consumes the contracted base
    // integrals (covl_sp, covl_sd) to build the target. The HRR namespace has no
    // operator abbreviation; the result keeps the operator-tagged covl name.
    const auto pp = read_file(dir / "ObaraSaikaTwoCenterOverlapPP.cpp");
    EXPECT_NE(pp.find("const auto covl_pp = os2c::hrr::compute_p_p(pair, covl_sp, covl_sd);"),
              std::string::npos);

    // (S|1|P), (P|1|S), (S|1|S): one side is zero -> no HRR.
    for (const auto* base : {"ObaraSaikaTwoCenterOverlapSP",
                             "ObaraSaikaTwoCenterOverlapPS",
                             "ObaraSaikaTwoCenterOverlapSS"})
    {
        const auto cpp = read_file(dir / (std::string(base) + ".cpp"));
        EXPECT_EQ(cpp.find("os2c::hrr::"), std::string::npos) << base;
    }
}

TEST(TwoCenterEmittersTest, CartesianTargetTransformedIntoSphericalBuffer)
{
    const auto dir = generate_in_temp_dir(overlap_config(1), "transform");

    // (P|1|P): the HRR target covl_pp is transformed (la=lb=1) into the buffer.
    const auto pp = read_file(dir / "ObaraSaikaTwoCenterOverlapPP.cpp");
    EXPECT_NE(pp.find("#include \"CartesianToSphericalFunc.hpp\""), std::string::npos);
    EXPECT_NE(pp.find("osfunc::transform<1, 1>(buffer, covl_pp);"), std::string::npos);
    // the transform precedes the return of the buffer
    EXPECT_LT(pp.find("osfunc::transform<1, 1>"), pp.find("return buffer;"));

    // (S|1|P): no HRR; the contracted base (s|p) is transformed directly.
    const auto sp = read_file(dir / "ObaraSaikaTwoCenterOverlapSP.cpp");
    EXPECT_NE(sp.find("osfunc::transform<0, 1>(buffer, covl_sp);"), std::string::npos);

    // (S|1|S) overlap is the special-cased path: no transform at all.
    const auto ss = read_file(dir / "ObaraSaikaTwoCenterOverlapSS.cpp");
    EXPECT_EQ(ss.find("osfunc::transform"), std::string::npos);
}

TEST(TwoCenterEmittersTest, MakeEmitterRejectsNothingForSupportedCpuCpp)
{
    // the default configuration targets cpu/C++, for which an emitter exists.
    EXPECT_NO_THROW({ generate_in_temp_dir(overlap_config(0), "supported"); });
}
