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

#include "test_axes.hpp"

#include "axes.hpp"

TEST_F(AxesTest, ToIndex)
{
    EXPECT_EQ(0, axes::to_index('x'));
    
    EXPECT_EQ(1, axes::to_index('y'));
    
    EXPECT_EQ(2, axes::to_index('z'));
    
    EXPECT_EQ(-1, axes::to_index('q'));
}

TEST_F(AxesTest, ToAxis)
{
    EXPECT_EQ('x', axes::to_axis(0));
    
    EXPECT_EQ('y', axes::to_axis(1));
    
    EXPECT_EQ('z', axes::to_axis(2));
    
    EXPECT_EQ('0', axes::to_axis(3));
}
