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

#include <string>
#include <iostream>
#include <chrono>
#include <utility>
#include <array>

#include "t2c_cpu_generators.hpp"
#include "t4c_cpu_generators.hpp"

int main(int argc, char **argv)
{
    // run configuration
    
    const auto run_type = std::pair<std::string, std::string>({"t4c_cpu", "electron repulsion"});
    
    const int max_ang_mom = 1;
    
    const auto rec_form = std::pair<bool, bool>({false, false});
    
    // set up start timer
    
    auto stime = std::chrono::high_resolution_clock::now();
    
    // case: two-center integrals
    
    if (run_type.first == "t2c_cpu")
    {
        const std::array<int, 3> geom_drvs = {0, 0, 0};
        
        const auto t2c_drv = T2CCPUGenerator();
        
        t2c_drv.generate(run_type.second, max_ang_mom, geom_drvs, rec_form);
    }
    
    // case: four-center integrals
    
    if (run_type.first == "t4c_cpu")
    {
        const std::array<int, 5> geom_drvs = {0, 0, 0, 0, 0};
        
        const auto t4c_drv = T4CCPUGenerator();
        
        t4c_drv.generate(run_type.second, max_ang_mom, geom_drvs);
    }
   
    // set up end timer & compute elapsed time
    
    auto etime = std::chrono::high_resolution_clock::now();
    
    auto dtime = std::chrono::duration_cast<std::chrono::seconds>(etime - stime);
    
    std::cout << "Elapsed time: " << dtime.count() << " seconds." << std::endl;
    
    return 0;
}
