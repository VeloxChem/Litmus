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

#include "test_components.hpp"

#include "components.hpp"
#include "tensor.hpp"
#include "setters.hpp"

TEST_F(ComponentsTest, Components)
{
    const VTensors vtens({Tensor(1), Tensor(2), Tensor(1)});
    
    auto vtcomps = make_components<TensorComponent>(vtens);
    
    EXPECT_EQ(vtcomps.size(), 54);
    
    int idx = 0;
    
    for (const auto& p0comp : gset::tensor_components(1))
    {
        for (const auto& dcomp : gset::tensor_components(2))
        {
            for (const auto& p1comp : gset::tensor_components(1))
            {
                EXPECT_EQ(vtcomps[idx], VTensorComponents({p0comp, dcomp, p1comp}));
            
                idx++;
            }
        }
    }
}
