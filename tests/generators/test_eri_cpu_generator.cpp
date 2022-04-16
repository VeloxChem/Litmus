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

TEST_F(EriCPUGeneratorTest, IsHrrRec)
{
    EriCPUGenerator gen_drv;
    
    const auto operi = Operator("1/|r-r'|");
    
    auto bpair = I2CPair("GA", 0, "GB", 2);
    
    auto kpair = I2CPair("GC", 0, "GD", 4);
        
    auto t4cint = I4CIntegral(bpair, kpair, operi);
    
    EXPECT_FALSE(gen_drv.is_hrr_rec(t4cint));
    
    kpair = I2CPair("GC", 1, "GD", 4);
        
    t4cint = I4CIntegral(bpair, kpair, operi);
    
    EXPECT_TRUE(gen_drv.is_hrr_rec(t4cint));
    
    bpair = I2CPair("GA", 3, "GB", 2);
        
    t4cint = I4CIntegral(bpair, kpair, operi);
    
    EXPECT_TRUE(gen_drv.is_hrr_rec(t4cint));
    
    t4cint = I4CIntegral();
    
    EXPECT_FALSE(gen_drv.is_hrr_rec(t4cint));
}

TEST_F(EriCPUGeneratorTest, IsVrrRec)
{
    EriCPUGenerator gen_drv;
    
    const auto operi = Operator("1/|r-r'|");
    
    auto bpair = I2CPair("GA", 0, "GB", 2);
    
    auto kpair = I2CPair("GC", 0, "GD", 4);
        
    auto t4cint = I4CIntegral(bpair, kpair, operi);
    
    EXPECT_TRUE(gen_drv.is_vrr_rec(t4cint));
    
    kpair = I2CPair("GC", 1, "GD", 4);
        
    t4cint = I4CIntegral(bpair, kpair, operi);
    
    EXPECT_FALSE(gen_drv.is_vrr_rec(t4cint));
    
    bpair = I2CPair("GA", 3, "GB", 2);
        
    t4cint = I4CIntegral(bpair, kpair, operi);
    
    EXPECT_FALSE(gen_drv.is_vrr_rec(t4cint));
    
    t4cint = I4CIntegral();
    
    EXPECT_FALSE(gen_drv.is_vrr_rec(t4cint));
}

TEST_F(EriCPUGeneratorTest, IsAuxRec)
{
    EriCPUGenerator gen_drv;
    
    const auto operi = Operator("1/|r-r'|");
    
    auto bpair = I2CPair("GA", 0, "GB", 2);
    
    auto kpair = I2CPair("GC", 0, "GD", 4);
        
    auto t4cint = I4CIntegral(bpair, kpair, operi);
    
    EXPECT_FALSE(gen_drv.is_aux_rec(t4cint));
    
    kpair = I2CPair("GC", 1, "GD", 4);
        
    t4cint = I4CIntegral(bpair, kpair, operi);
    
    EXPECT_FALSE(gen_drv.is_aux_rec(t4cint));
    
    bpair = I2CPair("GA", 3, "GB", 2);
        
    t4cint = I4CIntegral(bpair, kpair, operi);
    
    EXPECT_FALSE(gen_drv.is_aux_rec(t4cint));
    
    t4cint = I4CIntegral();
    
    EXPECT_TRUE(gen_drv.is_aux_rec(t4cint));
}
