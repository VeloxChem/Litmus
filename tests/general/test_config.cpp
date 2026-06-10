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

#include <gtest/gtest.h>

#include <vector>

#include "config.hpp"

using cfg::ConfigError;

TEST(ConfigTest, ParsesScalarTypes)
{
    const auto config = cfg::parse_string(R"(
        type      = "t2c_geom_cpu"
        lmax      = 2
        proj_lmax = -1
        use_rs    = true
        diagonal  = false
    )");

    EXPECT_EQ(config.get_string("type"), "t2c_geom_cpu");
    EXPECT_EQ(config.get_int("lmax"), 2);
    EXPECT_EQ(config.get_int("proj_lmax"), -1);
    EXPECT_TRUE(config.get_bool("use_rs"));
    EXPECT_FALSE(config.get_bool("diagonal"));
}

TEST(ConfigTest, ParsesBareString)
{
    // an unquoted, non-numeric, non-boolean value is taken as a string
    const auto config = cfg::parse_string("type = t4c_cpu");

    EXPECT_EQ(config.get_string("type"), "t4c_cpu");
}

TEST(ConfigTest, ParsesIntegerArrays)
{
    const auto config = cfg::parse_string(R"(
        geom = [1, 1, 0, 0, 0]
        one  = [3]
        none = []
    )");

    EXPECT_EQ(config.get_int_array("geom"), (std::vector<int>{1, 1, 0, 0, 0}));
    EXPECT_EQ(config.get_int_array("one"), (std::vector<int>{3}));
    EXPECT_TRUE(config.get_int_array("none").empty());
}

TEST(ConfigTest, IgnoresCommentsAndBlankLines)
{
    const auto config = cfg::parse_string(R"(
        # a leading comment
        type = "t3c_cpu"   # trailing comment after a value

        lmax = 3
    )");

    EXPECT_EQ(config.get_string("type"), "t3c_cpu");
    EXPECT_EQ(config.get_int("lmax"), 3);
}

TEST(ConfigTest, KeepsHashInsideQuotedString)
{
    const auto config = cfg::parse_string(R"(integral = "1/|r-r'| # eri")");

    EXPECT_EQ(config.get_string("integral"), "1/|r-r'| # eri");
}

TEST(ConfigTest, HasReportsPresence)
{
    const auto config = cfg::parse_string("lmax = 1");

    EXPECT_TRUE(config.has("lmax"));
    EXPECT_FALSE(config.has("type"));
}

TEST(ConfigTest, FallbackOverloadsUseDefaultWhenAbsent)
{
    const auto config = cfg::parse_string("lmax = 4");

    EXPECT_EQ(config.get_string("type", "none"), "none");
    EXPECT_EQ(config.get_int("proj_lmax", 7), 7);
    EXPECT_TRUE(config.get_bool("use_rs", true));
    EXPECT_EQ(config.get_int_array("geom", {0, 0, 0}), (std::vector<int>{0, 0, 0}));

    // present keys still read through the fallback overloads
    EXPECT_EQ(config.get_int("lmax", 99), 4);
}

TEST(ConfigTest, MissingKeyThrows)
{
    const auto config = cfg::parse_string("lmax = 1");

    EXPECT_THROW(config.get_string("type"), ConfigError);
    EXPECT_THROW(config.get_int("absent"), ConfigError);
}

TEST(ConfigTest, TypeMismatchThrows)
{
    const auto config = cfg::parse_string(R"(
        type = "t2c_cpu"
        lmax = 2
        geom = [0, 0, 0]
    )");

    EXPECT_THROW(config.get_int("type"), ConfigError);       // string read as int
    EXPECT_THROW(config.get_string("lmax"), ConfigError);    // int read as string
    EXPECT_THROW(config.get_bool("lmax"), ConfigError);      // int read as bool
    EXPECT_THROW(config.get_int_array("lmax"), ConfigError); // int read as array
}

TEST(ConfigTest, MalformedLineThrows)
{
    EXPECT_THROW(cfg::parse_string("type t2c_cpu"), ConfigError);   // no '='
    EXPECT_THROW(cfg::parse_string("= value"), ConfigError);        // empty key
    EXPECT_THROW(cfg::parse_string("lmax ="), ConfigError);         // missing value
    EXPECT_THROW(cfg::parse_string("geom = [1, 2"), ConfigError);   // unterminated array
    EXPECT_THROW(cfg::parse_string("geom = [1, x]"), ConfigError);  // non-integer element
    EXPECT_THROW(cfg::parse_string(R"(s = "abc)"), ConfigError);    // unterminated string
}

TEST(ConfigTest, MissingFileThrows)
{
    EXPECT_THROW(cfg::parse_file("/nonexistent/litmus/run.toml"), ConfigError);
}
