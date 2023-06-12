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

#include "spherical_momentum.hpp"

SphericalMomentum::SphericalMomentum(const int angmom)
{
    // s-type angular momentum

    if (angmom == 0)
    {
        _factors = {"1.0", };

        _indexes = {0, };
        
        _dimensions = {1, };
    }
    
    // p-type angular momentum

    if (angmom == 1)
    {
        // order: p_-1, p_0, p_1 i.e. p_y, p_z, p_x

        _factors = {"1.0", "1.0", "1.0"};

        _indexes = {1, 2, 0};
        
        _dimensions = {1, 1, 1};
    }
    
    // d-type angular momentum

    if (angmom == 2)
    {
        // order: d_-2, d_-1, d_0, d_1, d_2

        _factors = {"f2_3", "f2_3", "-1.0", "-1.0", "2.0", "f2_3", "0.5 * f2_3", "-0.5 * f2_3"};
    
        _indexes = {1, 4, 0, 3, 5, 2, 0, 3};
        
        _dimensions = {1, 1, 3, 1, 2};
    }
    
    if (angmom == 3)
    {
        // order: f_-3, f_-2, f_-1, f_0, f_1, f_2, f_3

        _factors = {"3.0 * f3_5", "-f3_5", "f3_15", "4.0 * f3_3", "-f3_3", "-f3_3", "2.0", "-3.0",
                    "-3.0", "4.0 * f3_3", "-f3_3", "-f3_3", "0.5 * f3_15", "-0.5 * f3_15", "f3_5",
                    "-3.0 * f3_5"};

        _indexes = {1, 6, 4, 8, 1, 6, 9, 2, 7, 5, 0, 3, 2, 7, 0, 3};
        
        _dimensions = {2, 1, 3, 3, 3, 2, 2};
    }
    
    // g-type angular momentum

    if (angmom == 4)
    {
        // order: g_-4, g_-3, g_-2, g_-1, g_0, g_1, g_2, g_3, g_4

        _factors = {"f4_35", "-f4_35", "3.0 * f4_17", "-f4_17", "6.0 * f4_5", "-f4_5", "-f4_5",
                    "4.0 * f4_2", "-3.0 * f4_2", "-3.0 * f4_2", "8.0", "3.0", "3.0", "6.0",
                    "-24.0", "-24.0", "4.0 * f4_2", "-3.0 * f4_2", "-3.0 * f4_2", "3.0 * f4_5",
                    "-3.0 * f4_5", "-0.5 * f4_5", "0.5 * f4_5", "f4_17", "-3.0 * f4_17",
                    "0.25 * f4_35", "0.25 * f4_35", "-1.50 * f4_35"};
                                    
        _indexes = {1, 6, 4, 11, 8, 1, 6, 13, 4, 11, 14, 0, 10, 3, 5, 12, 9, 2, 7, 5, 12, 0,
                    10, 2, 7, 0, 10, 3};
                                        
        _dimensions = {2, 2, 3, 3, 6, 3, 4, 2, 3};
    }
}

std::vector<std::pair<int, std::string>>
SphericalMomentum::select_pairs(const int cartcomp) const
{
    std::vector<std::pair<int, std::string>>  pairs;
    
    int idx = 0;
    
    for (int i = 0; i < _dimensions.size(); i++)
    {
        for (int j = 0; j < _dimensions[i]; j++)
        {
            if (cartcomp == _indexes[idx])
            {
                pairs.push_back({i, _factors[idx]});
            }
            
            idx++;
        }
    }
    
    return pairs;
}

std::vector<std::string>
SphericalMomentum::get_factors(const int angmom) const
{
    if (angmom == 2)
    {
        return std::vector<std::string>({"f2_3 = 2.0 * std::sqrt(3.0)", });
    }
    
    if (angmom == 3)
    {
        return std::vector<std::string>({"f3_5 = std::sqrt(2.5)",
                                         "f3_15 = 2.0 * std::sqrt(15.0)",
                                         "f3 = std::sqrt(1.5)"});
    }
    
    if (angmom == 4)
    {
        return std::vector<std::string>({"f4_35 = 4.0 * std::sqrt(35)",
                                         "f4_17 = 4.0 * std::sqrt(17.5)",
                                         "f4_5 = 4.0 * std::sqrt(5.0)",
                                         "f4_2 = 4.0 * std::sqrt(2.5)"});
    }
    
    return std::vector<std::string>();
}
