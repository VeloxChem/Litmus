// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <string>

#include "provenance.hpp"

namespace {

bool contains(const std::string& text, const std::string& needle)
{
    return text.find(needle) != std::string::npos;
}

}  // namespace

TEST(ProvenanceTest, BannerStartsWithVeloxChemLicense)
{
    ost::set_generation_config("");

    const auto banner = ost::generated_file_banner();

    EXPECT_EQ(banner.rfind("//\n//                                   VELOXCHEM\n", 0), 0u);
    EXPECT_TRUE(contains(banner, "SPDX-License-Identifier: BSD-3-Clause"));
    EXPECT_TRUE(contains(banner, "Copyright 2018-2025 VeloxChem developers"));
}

TEST(ProvenanceTest, BannerMarksFileAsGenerated)
{
    ost::set_generation_config("");

    const auto banner = ost::generated_file_banner();

    EXPECT_TRUE(contains(banner, "generated automatically by Litmus"));
    EXPECT_TRUE(contains(banner, "DO NOT EDIT"));
    EXPECT_TRUE(contains(banner, "// Litmus revision: " + ost::litmus_revision() + "\n"));
    EXPECT_FALSE(ost::litmus_revision().empty());
}

TEST(ProvenanceTest, BannerWithoutConfigOmitsRegenerationRecipe)
{
    ost::set_generation_config("");

    EXPECT_FALSE(contains(ost::generated_file_banner(), "To regenerate"));
}

TEST(ProvenanceTest, BannerCarriesConfigWithoutCommentsOrBlankLines)
{
    ost::set_generation_config("# A comment line\n\ntype = \"t2c_cpu\"   # family\n  lmax = 2\n");

    const auto banner = ost::generated_file_banner();

    EXPECT_TRUE(contains(banner, "//     litmus.x run <config-file>\n"));
    EXPECT_TRUE(contains(banner, "//     type = \"t2c_cpu\"   # family\n"));
    EXPECT_TRUE(contains(banner, "//     lmax = 2\n"));
    EXPECT_FALSE(contains(banner, "A comment line"));

    ost::set_generation_config("");
}

TEST(ProvenanceTest, OpenGeneratedFileWritesBannerFirst)
{
    ost::set_generation_config("lmax = 1\n");

    const std::string path = testing::TempDir() + "/litmus_provenance_test.hpp";

    {
        std::ofstream fstream;
        ost::open_generated_file(fstream, path);
        fstream << "#ifndef guard\n";
    }

    std::ifstream in(path);
    std::stringstream buffer;
    buffer << in.rdbuf();

    EXPECT_EQ(buffer.str(), ost::generated_file_banner() + "#ifndef guard\n");

    ost::set_generation_config("");
}
