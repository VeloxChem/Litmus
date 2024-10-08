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

#ifndef axes_hpp
#define axes_hpp

namespace axes
{
/// Converts axis label to index.
/// @param axis The label of axis.
/// @return The index of given axis.
inline int to_index(const char axis)
{
    if (axis == 'x') return 0;
    
    if (axis == 'y') return 1;
    
    if (axis == 'z') return 2;
    
    return -1;
}

/// Converts index to axis label.
/// @param index The index of axis.
/// @return The label of axis.
inline char to_axis(const int index)
{
    if (index == 0) return 'x';
    
    if (index == 1) return 'y';
    
    if (index == 2) return 'z';
    
    return '0';
}

}

#endif /* axes_hpp */
