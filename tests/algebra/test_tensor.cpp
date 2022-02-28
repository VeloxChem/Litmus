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

#include "test_tensor.hpp"

#include "tensor.hpp"

#include "setters.hpp"

TEST_F(TensorTest, Constructor)
{
    EXPECT_EQ(Tensor(), Tensor(0));
    
    for (const auto& tcomp : gset::tensor_components(3))
    {
        EXPECT_EQ(Tensor(3), Tensor(tcomp));
    }
}

TEST_F(TensorTest, OperatorEqual)
{
    for (const auto& tcomp : gset::tensor_components(3))
    {
        EXPECT_TRUE(Tensor(3) == Tensor(tcomp));
    }
}

TEST_F(TensorTest, OperatorNotEqual)
{
    for (const auto& tcomp : gset::tensor_components(3))
    {
        EXPECT_TRUE(Tensor(2) != Tensor(tcomp));
    }
}

TEST_F(TensorTest, OperatorLess)
{
    // compare to S tensor
    
    EXPECT_FALSE(Tensor(0) < Tensor(0));
    
    EXPECT_FALSE(Tensor(1) < Tensor(0));
    
    EXPECT_FALSE(Tensor(2) < Tensor(0));
    
    EXPECT_FALSE(Tensor(3) < Tensor(0));
    
    // compare to P tensor
    
    EXPECT_TRUE(Tensor(0) < Tensor(1));
    
    EXPECT_FALSE(Tensor(1) < Tensor(1));
    
    EXPECT_FALSE(Tensor(2) < Tensor(1));
    
    EXPECT_FALSE(Tensor(3) < Tensor(1));
    
    // compare to D tensor
    
    EXPECT_TRUE(Tensor(0) < Tensor(2));
    
    EXPECT_TRUE(Tensor(1) < Tensor(2));
    
    EXPECT_FALSE(Tensor(2) < Tensor(2));
    
    EXPECT_FALSE(Tensor(3) < Tensor(2));
    
    // compare to F tensor
    
    EXPECT_TRUE(Tensor(0) < Tensor(3));
    
    EXPECT_TRUE(Tensor(1) < Tensor(3));
    
    EXPECT_TRUE(Tensor(2) < Tensor(3));
    
    EXPECT_FALSE(Tensor(3) < Tensor(3));
}

TEST_F(TensorTest, Label)
{
    const std::string names("SPDFGHIKLMNOQRTUV");
    
    for (int i = 0; i < 17; i++)
    {
        const Tensor tval(i);
        
        EXPECT_EQ(tval.label(), std::string(1, names[i]));
    }
        
    const Tensor tval(17);
       
    EXPECT_EQ(tval.label(), "l17");
}

TEST_F(TensorTest, Components)
{
    for (int i = 0; i < 4; i++)
    {
        const Tensor tval(i);
        
        EXPECT_EQ(tval.components(), gset::tensor_components(i));
    }
}
