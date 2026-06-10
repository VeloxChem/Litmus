// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include <fstream>
#include <sstream>
#include <string>

#include "file_stream.hpp"

namespace {

// Writes the given code lines through ost::write_code_lines into a temporary
// file and returns the resulting file contents verbatim.
std::string render(const VCodeLines& lines)
{
    const std::string path = testing::TempDir() + "/litmus_file_stream_test.txt";

    {
        std::ofstream fstream(path);
        ost::write_code_lines(fstream, lines);
    }

    std::ifstream in(path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

}  // namespace

TEST(FileStreamTest, EmptyLinesProduceEmptyOutput)
{
    EXPECT_EQ(render(VCodeLines({})), "");
}

TEST(FileStreamTest, IndentationIsFourSpacesPerSpacerPlusOffset)
{
    // {nspacers, offset, nends, text}: 1 spacer (4 spaces) + 2 offset = 6 spaces.
    EXPECT_EQ(render(VCodeLines({{1, 2, 1, "code"}})), "      code\n");
}

TEST(FileStreamTest, ZeroSpacersAndOffsetGivesNoIndent)
{
    EXPECT_EQ(render(VCodeLines({{0, 0, 1, "x"}})), "x\n");
}

TEST(FileStreamTest, EndsControlNewlineCount)
{
    EXPECT_EQ(render(VCodeLines({{0, 0, 0, "noeol"}})), "noeol");
    EXPECT_EQ(render(VCodeLines({{0, 0, 2, "blank"}})), "blank\n\n");
}

TEST(FileStreamTest, MultipleLinesAreConcatenatedInOrder)
{
    const VCodeLines lines({{1, 0, 1, "first"}, {0, 0, 2, "second"}});

    EXPECT_EQ(render(lines), "    first\nsecond\n\n");
}
