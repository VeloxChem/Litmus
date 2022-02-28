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

#include "setters.hpp"

namespace gset
{

VTensorComponents
tensor_components(const int order)
{
    VTensorComponents tcomps;
    
    // S angular momentum components
    
    if (order == 0)
    {
        tcomps.push_back(TensorComponent(0, 0, 0));
        
        return tcomps;
    }
    
    // P angular momentum components
    
    if (order == 1)
    {
        tcomps.push_back(TensorComponent(1, 0, 0));
        
        tcomps.push_back(TensorComponent(0, 1, 0));
        
        tcomps.push_back(TensorComponent(0, 0, 1));
        
        return tcomps;
    }
    
    // D angular momentum components
    
    if (order == 2)
    {
        tcomps.push_back(TensorComponent(2, 0, 0));
        
        tcomps.push_back(TensorComponent(1, 1, 0));
        
        tcomps.push_back(TensorComponent(1, 0, 1));
        
        tcomps.push_back(TensorComponent(0, 2, 0));
        
        tcomps.push_back(TensorComponent(0, 1, 1));
        
        tcomps.push_back(TensorComponent(0, 0, 2));
        
        return tcomps;
    }
    
    // F angular momentum components
    
    if (order == 3)
    {
        tcomps.push_back(TensorComponent(3, 0, 0));
        
        tcomps.push_back(TensorComponent(2, 1, 0));
        
        tcomps.push_back(TensorComponent(2, 0, 1));
        
        tcomps.push_back(TensorComponent(1, 2, 0));
        
        tcomps.push_back(TensorComponent(1, 1, 1));
        
        tcomps.push_back(TensorComponent(1, 0, 2));
        
        tcomps.push_back(TensorComponent(0, 3, 0));
        
        tcomps.push_back(TensorComponent(0, 2, 1));
        
        tcomps.push_back(TensorComponent(0, 1, 2));
        
        tcomps.push_back(TensorComponent(0, 0, 3));
        
        return tcomps;
    }
    
    return tcomps;
}

}
