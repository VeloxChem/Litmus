// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "one_center_component.hpp"

TEST(OneCenterComponentTest, DefaultConstruction)
{
    const OneCenterComponent center;

    EXPECT_EQ(center.name(), "");
    EXPECT_EQ(center.shape(), TensorComponent(0, 0, 0));
    EXPECT_EQ(center.centers(), 1);
}

TEST(OneCenterComponentTest, NameShapeAndIndexing)
{
    const OneCenterComponent center("a", TensorComponent(1, 0, 0));

    EXPECT_EQ(center.name(), "a");
    EXPECT_EQ(center.shape(), TensorComponent(1, 0, 0));
    EXPECT_EQ(center[0], TensorComponent(1, 0, 0));  // shape is returned for any index
}

TEST(OneCenterComponentTest, Equality)
{
    const OneCenterComponent a("a", TensorComponent(1, 0, 0));

    EXPECT_EQ(a, OneCenterComponent("a", TensorComponent(1, 0, 0)));
    EXPECT_NE(a, OneCenterComponent("b", TensorComponent(1, 0, 0)));
    EXPECT_NE(a, OneCenterComponent("a", TensorComponent(0, 1, 0)));
}

TEST(OneCenterComponentTest, OrderingByNameThenShape)
{
    EXPECT_TRUE(OneCenterComponent("a", TensorComponent(2, 0, 0)) <
                OneCenterComponent("b", TensorComponent(0, 0, 0)));  // name dominates
    EXPECT_TRUE(OneCenterComponent("a", TensorComponent(0, 1, 0)) <
                OneCenterComponent("a", TensorComponent(1, 0, 0)));  // then shape
}

TEST(OneCenterComponentTest, SimilarityComparesOrderOnly)
{
    const OneCenterComponent a("a", TensorComponent(1, 0, 0));

    EXPECT_TRUE(a.similar(OneCenterComponent("a", TensorComponent(0, 1, 0))));   // same order
    EXPECT_FALSE(a.similar(OneCenterComponent("a", TensorComponent(2, 0, 0))));  // higher order
    EXPECT_FALSE(a.similar(OneCenterComponent("b", TensorComponent(1, 0, 0))));  // different name
}

TEST(OneCenterComponentTest, Labels)
{
    EXPECT_EQ(OneCenterComponent("a", TensorComponent(0, 0, 0)).label(), "0");
    EXPECT_EQ(OneCenterComponent("a", TensorComponent(1, 1, 0)).label(), "xy");
}

TEST(OneCenterComponentTest, ToString)
{
    EXPECT_EQ(OneCenterComponent("a", TensorComponent(1, 0, 0)).to_string(), "{a:(1,0,0)}");
}

TEST(OneCenterComponentTest, ShiftRaisesAndLowers)
{
    const OneCenterComponent a("a", TensorComponent(1, 0, 0));

    const auto raised = a.shift('x', 1, 0);
    ASSERT_TRUE(raised.has_value());
    EXPECT_EQ(raised->shape(), TensorComponent(2, 0, 0));

    EXPECT_FALSE(a.shift('y', -1, 0).has_value());  // would go negative
}
