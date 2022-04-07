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

int main(int argc, char **argv)
{
    // set up start timer
    
    auto stime = std::chrono::high_resolution_clock::now();
    
    // four-center integrals repository
    
    Repository<R4Group, T4CIntegral> t4c_repo;
    
    // set up angular momentum values
    
    const int mang = 2;
    
    // electron repulsion integral recursions
    
    if (true)
    {
        EriDriver eri_drv;
        
        t4c_repo.add(eri_drv.create_graphs(mang, false));
        
        t4c_repo.add(eri_drv.create_graphs(mang, true));
        
        //t4c_repo.add(eri_drv.create_graphs(mang, mang, mang, mang, true));
        
        // testing
        
//        auto vgraphs = eri_drv.create_graph(4, 4, 4, 4, false);
        
//        std::cout << vgraphs.vertices() << std::endl;
    }
    
    // print summary of integrals repository

    t4c_repo.summary();
    
    t4c_repo.details<I4CIntegral>();
 
    // set up end timer & compute elapsed time
    
    auto etime = std::chrono::high_resolution_clock::now();
    
    auto dtime = std::chrono::duration_cast<std::chrono::seconds>(etime - stime);
    
    std::cout << "Elapsed time: " << dtime.count() << " seconds." << std::endl;
    
    return 0;
}
