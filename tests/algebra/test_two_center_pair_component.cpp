// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "two_center_pair_component.hpp"

namespace {

TwoCenterPairComponent make_pair(const TensorComponent& f, const TensorComponent& s)
{
    return TwoCenterPairComponent({"a", "b"}, {f, s});
}

}  // namespace

TEST(TwoCenterPairComponentTest, DefaultConstruction)
{
    const TwoCenterPairComponent pair;

    EXPECT_EQ(pair.names(), (std::array<std::string, 2>({"", ""})));
    EXPECT_EQ(pair.centers(), 2);
    EXPECT_EQ(pair[0], TensorComponent(0, 0, 0));
    EXPECT_EQ(pair[1], TensorComponent(0, 0, 0));
}

TEST(TwoCenterPairComponentTest, AccessorsAndIndexing)
{
    const auto pair = make_pair(TensorComponent(1, 0, 0), TensorComponent(0, 1, 0));

    EXPECT_EQ(pair.names(), (std::array<std::string, 2>({"a", "b"})));
    EXPECT_EQ(pair[0], TensorComponent(1, 0, 0));
    EXPECT_EQ(pair[1], TensorComponent(0, 1, 0));
}

TEST(TwoCenterPairComponentTest, Equality)
{
    const auto pair = make_pair(TensorComponent(1, 0, 0), TensorComponent(0, 1, 0));

    EXPECT_EQ(pair, make_pair(TensorComponent(1, 0, 0), TensorComponent(0, 1, 0)));
    EXPECT_NE(pair, make_pair(TensorComponent(1, 0, 0), TensorComponent(0, 0, 1)));
    EXPECT_NE(pair, TwoCenterPairComponent({"a", "c"},
                                           {TensorComponent(1, 0, 0), TensorComponent(0, 1, 0)}));
}

TEST(TwoCenterPairComponentTest, OrderingByNamesThenShapes)
{
    EXPECT_TRUE(TwoCenterPairComponent({"a", "b"}, {TensorComponent(2, 0, 0), TensorComponent(0, 0, 0)}) <
                TwoCenterPairComponent({"a", "c"}, {TensorComponent(0, 0, 0), TensorComponent(0, 0, 0)}));
    EXPECT_TRUE(make_pair(TensorComponent(0, 1, 0), TensorComponent(0, 0, 0)) <
                make_pair(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)));
}

TEST(TwoCenterPairComponentTest, SimilarityPerCenter)
{
    const auto pair = make_pair(TensorComponent(1, 0, 0), TensorComponent(2, 0, 0));

    EXPECT_TRUE(pair.similar(make_pair(TensorComponent(0, 1, 0), TensorComponent(0, 2, 0))));
    EXPECT_FALSE(pair.similar(make_pair(TensorComponent(2, 0, 0), TensorComponent(2, 0, 0))));
    EXPECT_FALSE(pair.similar(make_pair(TensorComponent(1, 0, 0), TensorComponent(1, 0, 0))));
}

TEST(TwoCenterPairComponentTest, Label)
{
    EXPECT_EQ(make_pair(TensorComponent(1, 0, 0), TensorComponent(0, 0, 0)).label(), "x_0");
    EXPECT_EQ(make_pair(TensorComponent(1, 1, 0), TensorComponent(0, 0, 1)).label(), "xy_z");
}

TEST(TwoCenterPairComponentTest, ToString)
{
    EXPECT_EQ(make_pair(TensorComponent(1, 0, 0), TensorComponent(0, 1, 0)).to_string(),
              "{a:(1,0,0);b:(0,1,0)}");
}

TEST(TwoCenterPairComponentTest, ShiftTargetsSelectedCenter)
{
    const auto pair = make_pair(TensorComponent(1, 0, 0), TensorComponent(0, 1, 0));

    const auto first = pair.shift('x', 1, 0);
    ASSERT_TRUE(first.has_value());
    EXPECT_EQ((*first)[0], TensorComponent(2, 0, 0));
    EXPECT_EQ((*first)[1], TensorComponent(0, 1, 0));  // untouched

    const auto second = pair.shift('y', 1, 1);
    ASSERT_TRUE(second.has_value());
    EXPECT_EQ((*second)[1], TensorComponent(0, 2, 0));

    EXPECT_FALSE(pair.shift('z', -1, 0).has_value());  // would go negative
}
