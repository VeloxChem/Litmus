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

#include <iostream>
#include <chrono>

#include "eri_driver.hpp"
#include "repository.hpp"
#include "eri_cpu_generator.hpp"

int main(int argc, char **argv)
{
    // set up start timer
    
    auto stime = std::chrono::high_resolution_clock::now();
    
    // four-center integrals repository
    
    Repository<R4Group, T4CIntegral> t4c_repo;
    
    // set up integrals generator parameters
    
    const int mang = 3;
    
    const bool diag_form = false;
    
    // electron repulsion integral recursions
    
    if (true)
    {
        EriDriver eri_drv;
        
        const auto graphs = eri_drv.create_graphs(mang, diag_form);
        
        // t4c_repo.add(eri_drv.create_graphs(mang, diag_form));
        
        EriCPUGenerator gen_drv;
        
        if (diag_form) gen_drv.set_diag_form();
    
        //gen_drv.generate(t4c_repo);
    }
    
    // print summary of integrals repository

    //t4c_repo.summary();
    
    //t4c_repo.details<I4CIntegral>();
 
    // set up end timer & compute elapsed time
    
    auto etime = std::chrono::high_resolution_clock::now();
    
    auto dtime = std::chrono::duration_cast<std::chrono::seconds>(etime - stime);
    
    std::cout << "Elapsed time: " << dtime.count() << " seconds." << std::endl;
    
    return 0;
}
