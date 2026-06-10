// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.

#include <gtest/gtest.h>

#include "components.hpp"
#include "tensor.hpp"
#include "tensor_component.hpp"

// make_components forms the direct product of the component lists of a vector of
// tensor-like objects, returning one sub-vector per element of the product.

TEST(ComponentsTest, EmptyInputYieldsEmptyProduct)
{
    const auto product = make_components<TensorComponent>(std::vector<Tensor>({}));

    EXPECT_TRUE(product.empty());
}

TEST(ComponentsTest, SingleTensorMirrorsItsComponents)
{
    const auto product = make_components<TensorComponent>(std::vector<Tensor>({Tensor(1)}));

    // A P shell has three Cartesian components, each wrapped in a 1-element vector.
    ASSERT_EQ(product.size(), 3u);

    for (const auto& vcomp : product)
    {
        EXPECT_EQ(vcomp.size(), 1u);
    }
}

TEST(ComponentsTest, ProductSizeIsCartesianProduct)
{
    // P (3 components) x P (3 components) -> 9 pairs.
    const auto product =
        make_components<TensorComponent>(std::vector<Tensor>({Tensor(1), Tensor(1)}));

    ASSERT_EQ(product.size(), 9u);

    for (const auto& vcomp : product)
    {
        EXPECT_EQ(vcomp.size(), 2u);
    }
}

TEST(ComponentsTest, ScalarFactorDoesNotChangeCount)
{
    // S (1 component) x P (3 components) -> 3 pairs of length two.
    const auto product =
        make_components<TensorComponent>(std::vector<Tensor>({Tensor(0), Tensor(1)}));

    ASSERT_EQ(product.size(), 3u);
    EXPECT_EQ(product[0].size(), 2u);
}
