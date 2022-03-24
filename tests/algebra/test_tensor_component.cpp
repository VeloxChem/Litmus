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

#include "test_tensor_component.hpp"

#include "tensor_component.hpp"

#include "setters.hpp"

TEST_F(TensorComponentTest, Constructor)
{
    EXPECT_EQ(TensorComponent(), TensorComponent(0, 0, 0));
}

TEST_F(TensorComponentTest, OperatorBrackets)
{
    const TensorComponent tcomp(1, 2, 3);
    
    EXPECT_EQ(tcomp['x'], 1);
    
    EXPECT_EQ(tcomp['y'], 2);
    
    EXPECT_EQ(tcomp['z'], 3);
    
    EXPECT_EQ(tcomp['g'], -1);
}

TEST_F(TensorComponentTest, OperatorEqual)
{
    const auto pcomps = gset::tensor_components(1);
    
    EXPECT_TRUE(pcomps[0] == TensorComponent(1, 0, 0));
    
    EXPECT_TRUE(pcomps[1] == TensorComponent(0, 1, 0));
    
    EXPECT_TRUE(pcomps[2] == TensorComponent(0, 0, 1));
}

TEST_F(TensorComponentTest, OperatorNotEqual)
{
    const auto scomps = gset::tensor_components(0);
    
    EXPECT_TRUE(scomps[0] != TensorComponent(1, 0, 0));
    
    EXPECT_TRUE(scomps[0] != TensorComponent(0, 1, 0));
    
    EXPECT_TRUE(scomps[0] != TensorComponent(0, 0, 1));
}

TEST_F(TensorComponentTest, OperatorLess)
{
    const auto dcomps = gset::tensor_components(2);
    
    // compare to d_xx component
    
    EXPECT_FALSE(dcomps[0] < dcomps[0]);
    
    EXPECT_TRUE(dcomps[1] < dcomps[0]);
    
    EXPECT_TRUE(dcomps[2] < dcomps[0]);
    
    EXPECT_TRUE(dcomps[3] < dcomps[0]);
    
    EXPECT_TRUE(dcomps[4] < dcomps[0]);
    
    EXPECT_TRUE(dcomps[5] < dcomps[0]);
    
    // compare to d_xy component
    
    EXPECT_FALSE(dcomps[0] < dcomps[1]);
    
    EXPECT_FALSE(dcomps[1] < dcomps[1]);
    
    EXPECT_TRUE(dcomps[2] < dcomps[1]);
    
    EXPECT_TRUE(dcomps[3] < dcomps[1]);
    
    EXPECT_TRUE(dcomps[4] < dcomps[1]);
    
    EXPECT_TRUE(dcomps[5] < dcomps[1]);
    
    // compare to d_xz component
    
    EXPECT_FALSE(dcomps[0] < dcomps[2]);
    
    EXPECT_FALSE(dcomps[1] < dcomps[2]);
    
    EXPECT_FALSE(dcomps[2] < dcomps[2]);
    
    EXPECT_TRUE(dcomps[3] < dcomps[2]);
    
    EXPECT_TRUE(dcomps[4] < dcomps[2]);
    
    EXPECT_TRUE(dcomps[5] < dcomps[2]);
    
    // compare to d_yy component
    
    EXPECT_FALSE(dcomps[0] < dcomps[3]);
    
    EXPECT_FALSE(dcomps[1] < dcomps[3]);
    
    EXPECT_FALSE(dcomps[2] < dcomps[3]);
    
    EXPECT_FALSE(dcomps[3] < dcomps[3]);
    
    EXPECT_TRUE(dcomps[4] < dcomps[3]);
    
    EXPECT_TRUE(dcomps[5] < dcomps[3]);
    
    // compare to d_yz component
    
    EXPECT_FALSE(dcomps[0] < dcomps[4]);
    
    EXPECT_FALSE(dcomps[1] < dcomps[4]);
    
    EXPECT_FALSE(dcomps[2] < dcomps[4]);
    
    EXPECT_FALSE(dcomps[3] < dcomps[4]);
    
    EXPECT_FALSE(dcomps[4] < dcomps[4]);
    
    EXPECT_TRUE(dcomps[5] < dcomps[4]);
    
    // compare to d_zz component
    
    EXPECT_FALSE(dcomps[0] < dcomps[5]);
    
    EXPECT_FALSE(dcomps[1] < dcomps[5]);
    
    EXPECT_FALSE(dcomps[2] < dcomps[5]);
    
    EXPECT_FALSE(dcomps[3] < dcomps[5]);
    
    EXPECT_FALSE(dcomps[4] < dcomps[5]);
    
    EXPECT_FALSE(dcomps[5] < dcomps[5]);
}

TEST_F(TensorComponentTest, Similar)
{
    const auto dcomps = gset::tensor_components(2);
    
    const auto fcomps = gset::tensor_components(3);
    
    // D tensor components
    
    for (const auto& lhs_comp : dcomps)
    {
        for (const auto rhs_comp : dcomps)
        {
            EXPECT_TRUE(lhs_comp.similar(rhs_comp));
        }
    }
    
    // F tensor components
    
    for (const auto& lhs_comp : fcomps)
    {
        for (const auto rhs_comp : fcomps)
        {
            EXPECT_TRUE(lhs_comp.similar(rhs_comp));
        }
    }
    
    // D, F tensor components
    
    for (const auto& lhs_comp : fcomps)
    {
        for (const auto rhs_comp : dcomps)
        {
            EXPECT_FALSE(lhs_comp.similar(rhs_comp));
        }
    }
}

TEST_F(TensorComponentTest, ToString)
{
    const auto fcomps = gset::tensor_components(3);
    
    EXPECT_EQ(fcomps[0].to_string(), "(3,0,0)");
    
    EXPECT_EQ(fcomps[1].to_string(), "(2,1,0)");
    
    EXPECT_EQ(fcomps[2].to_string(), "(2,0,1)");
    
    EXPECT_EQ(fcomps[3].to_string(), "(1,2,0)");
    
    EXPECT_EQ(fcomps[4].to_string(), "(1,1,1)");
    
    EXPECT_EQ(fcomps[5].to_string(), "(1,0,2)");
    
    EXPECT_EQ(fcomps[6].to_string(), "(0,3,0)");
    
    EXPECT_EQ(fcomps[7].to_string(), "(0,2,1)");
    
    EXPECT_EQ(fcomps[8].to_string(), "(0,1,2)");
    
    EXPECT_EQ(fcomps[9].to_string(), "(0,0,3)");
}

TEST_F(TensorComponentTest, Label)
{
    const auto fcomps = gset::tensor_components(3);
    
    EXPECT_EQ(fcomps[0].label(), "xxx");
    
    EXPECT_EQ(fcomps[1].label(), "xxy");
    
    EXPECT_EQ(fcomps[2].label(), "xxz");
    
    EXPECT_EQ(fcomps[3].label(), "xyy");
    
    EXPECT_EQ(fcomps[4].label(), "xyz");
    
    EXPECT_EQ(fcomps[5].label(), "xzz");
    
    EXPECT_EQ(fcomps[6].label(), "yyy");
    
    EXPECT_EQ(fcomps[7].label(), "yyz");
    
    EXPECT_EQ(fcomps[8].label(), "yzz");
    
    EXPECT_EQ(fcomps[9].label(), "zzz");
}

TEST_F(TensorComponentTest, Order)
{
    for (int32_t i = 0; i < 4; i++)
    {
        for (const auto& tcomp : gset::tensor_components(i))
        {
            EXPECT_EQ(tcomp.order(), i);
        }
    }
}

TEST_F(TensorComponentTest, Maximum)
{
    const auto fcomps = gset::tensor_components(3);
    
    EXPECT_EQ(fcomps[0].maximum(), 3);
    
    EXPECT_EQ(fcomps[1].maximum(), 2);
    
    EXPECT_EQ(fcomps[2].maximum(), 2);
    
    EXPECT_EQ(fcomps[3].maximum(), 2);
    
    EXPECT_EQ(fcomps[4].maximum(), 1);
    
    EXPECT_EQ(fcomps[5].maximum(), 2);
    
    EXPECT_EQ(fcomps[6].maximum(), 3);
    
    EXPECT_EQ(fcomps[7].maximum(), 2);
    
    EXPECT_EQ(fcomps[8].maximum(), 2);
    
    EXPECT_EQ(fcomps[9].maximum(), 3);
}

TEST_F(TensorComponentTest, Primary)
{
    const auto fcomps = gset::tensor_components(3);
    
    EXPECT_EQ(fcomps[0].primary(), 'x');
    
    EXPECT_EQ(fcomps[1].primary(), 'x');
    
    EXPECT_EQ(fcomps[2].primary(), 'x');
    
    EXPECT_EQ(fcomps[3].primary(), 'x');
    
    EXPECT_EQ(fcomps[4].primary(), 'x');
    
    EXPECT_EQ(fcomps[5].primary(), 'x');
    
    EXPECT_EQ(fcomps[6].primary(), 'y');
    
    EXPECT_EQ(fcomps[7].primary(), 'y');
    
    EXPECT_EQ(fcomps[8].primary(), 'y');
    
    EXPECT_EQ(fcomps[9].primary(), 'z');
}

TEST_F(TensorComponentTest, Shift)
{
    const auto dcomps = gset::tensor_components(2);
    
    const auto fcomps = gset::tensor_components(3);
    
    EXPECT_EQ(fcomps[0].shift('x', -1), dcomps[0]);
    
    EXPECT_EQ(fcomps[1].shift('x', -1), dcomps[1]);
    
    EXPECT_EQ(fcomps[1].shift('y', -1), dcomps[0]);
    
    EXPECT_EQ(fcomps[4].shift('x', -1), dcomps[4]);
    
    EXPECT_EQ(fcomps[4].shift('y', -1), dcomps[2]);
    
    EXPECT_EQ(fcomps[4].shift('z', -1), dcomps[1]);
    
    EXPECT_FALSE(fcomps[4].shift('x', -2));
    
    EXPECT_FALSE(fcomps[4].shift('y', -2));
    
    EXPECT_FALSE(fcomps[4].shift('z', -2));
}

