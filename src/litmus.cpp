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

#include "t2c_ovl_driver.hpp"
#include "eri_driver.hpp"
#include "repository.hpp"
#include "t2c_cpu_generator.hpp"

int main(int argc, char **argv)
{
    // select run type
    
    const auto run_type = std::pair<std::string, std::string>({"t2c", "Multipole"});
    
    const int max_angmom = 4;
    
    const int op_gdrv = 3;
    
    // set up start timer
    
    auto stime = std::chrono::high_resolution_clock::now();
    
    if (run_type.first == "t2c")
    {
        const auto t2c_drv = T2CCPUGenerator();
        
        t2c_drv.generate(run_type.second, max_angmom, op_gdrv);
    }
 
    // set up end timer & compute elapsed time
    
    auto etime = std::chrono::high_resolution_clock::now();
    
    auto dtime = std::chrono::duration_cast<std::chrono::seconds>(etime - stime);
    
    std::cout << "Elapsed time: " << dtime.count() << " seconds." << std::endl;
    
    return 0;
}
