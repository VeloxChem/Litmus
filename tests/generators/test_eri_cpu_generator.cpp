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

#include "test_eri_cpu_generator.hpp"

#include "eri_cpu_generator.hpp"
#include "eri_driver.hpp"

TEST_F(EriCPUGeneratorTest, IsHrrRecGroup)
{
    EriDriver eri_drv;
    
    const auto rgraph = eri_drv.create_graph(1, 1, 1, 1, true);
    
    EriCPUGenerator gen_drv;
    
    EXPECT_EQ(rgraph.vertices(), 22);
    
    EXPECT_TRUE(gen_drv.is_hrr_rec_group(rgraph[0]));
    
    EXPECT_TRUE(gen_drv.is_hrr_rec_group(rgraph[1]));
    
    EXPECT_FALSE(gen_drv.is_hrr_rec_group(rgraph[2]));
    
    EXPECT_FALSE(gen_drv.is_hrr_rec_group(rgraph[3]));
    
    EXPECT_TRUE(gen_drv.is_hrr_rec_group(rgraph[4]));
    
    for (int i = 5; i < 22; i++)
    {
        EXPECT_FALSE(gen_drv.is_hrr_rec_group(rgraph[i]));
    }
}

TEST_F(EriCPUGeneratorTest, IsVrrRecGroup)
{
    EriDriver eri_drv;
    
    const auto rgraph = eri_drv.create_graph(1, 1, 1, 1, true);
    
    EriCPUGenerator gen_drv;
    
    EXPECT_EQ(rgraph.vertices(), 22);
    
    EXPECT_FALSE(gen_drv.is_vrr_rec_group(rgraph[0]));
    
    EXPECT_FALSE(gen_drv.is_vrr_rec_group(rgraph[1]));
    
    EXPECT_TRUE(gen_drv.is_vrr_rec_group(rgraph[2]));
    
    EXPECT_TRUE(gen_drv.is_vrr_rec_group(rgraph[3]));
    
    EXPECT_FALSE(gen_drv.is_vrr_rec_group(rgraph[4]));
    
    for (int i = 5; i < 17; i++)
    {
        EXPECT_TRUE(gen_drv.is_vrr_rec_group(rgraph[i]));
    }
    
    for (int i = 17; i < 22; i++)
    {
        EXPECT_FALSE(gen_drv.is_vrr_rec_group(rgraph[i]));
    }
}

TEST_F(EriCPUGeneratorTest, IsAuxRecGroup)
{
    EriDriver eri_drv;
    
    const auto rgraph = eri_drv.create_graph(1, 1, 1, 1, true);
    
    EriCPUGenerator gen_drv;
    
    EXPECT_EQ(rgraph.vertices(), 22);
    
    for (int i = 0; i < 17; i++)
    {
        EXPECT_FALSE(gen_drv.is_aux_rec_group(rgraph[i]));
    }
    
    for (int i = 17; i < 22; i++)
    {
        EXPECT_TRUE(gen_drv.is_aux_rec_group(rgraph[i]));
    }
}
