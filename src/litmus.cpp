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
#include "t2c_geom_cpu_generators.hpp"
#include "t2c_geom_deriv_cpu_generators.hpp"
#include "t4c_cpu_generators.hpp"
#include "t4c_diag_cpu_generators.hpp"
#include "t4c_geom_cpu_generators.hpp"
#include "t4c_geom_hrr_cpu_generators.hpp"
#include "t4c_geom_deriv_cpu_generators.hpp"
#include "t4c_eri_tree_generators.hpp"
#include "t3c_cpu_generators.hpp"
#include "t3c_geom_cpu_generators.hpp"
#include "t3c_geom_hrr_cpu_generators.hpp"
#include "g2c_cpu_generators.hpp"
#include "t2c_ecp_cpu_generators.hpp"

int main(int argc, char **argv)
{
    // run configuration

    const auto run_type = std::pair<std::string, std::string>({"t2c_ecp_cpu", "local"});

    const int max_ang_mom = 4;

    // set up start timer
    
    auto stime = std::chrono::high_resolution_clock::now();
    
    // case: two-center integrals

    if (run_type.first == "t2c_cpu")
    {
        std::array<int, 3> geom_drvs = {0, 0, 0};
        
        const auto rec_form = std::pair<bool, bool>({true, false});
        
        const auto use_rs = false;
        
        if ((geom_drvs[0] + geom_drvs[2]) == 0)
        {
            const auto t2c_drv = T2CCPUGenerator();
            
            t2c_drv.generate(run_type.second, max_ang_mom, geom_drvs, rec_form, use_rs);
        }
        else
        {
            const auto t2c_drv = T2CGeomCPUGenerator();
            
            t2c_drv.generate(run_type.second, max_ang_mom, geom_drvs, rec_form, use_rs);
        }
    }
    
    // case: four-center integrals
    
    if (run_type.first == "t4c_cpu")
    {
        // a, b, operator, c, d
        std::array<int, 5> geom_drvs = {1, 1, 0, 0, 0};
        
        if (geom_drvs == std::array<int, 5>({0, 0, 0, 0, 0}))
        {
            const auto t4c_drv = T4CCPUGenerator();
            
            t4c_drv.generate(run_type.second, max_ang_mom);
        }
        else
        {
            const auto t4c_drv = T4CGeomCPUGenerator();
            
            t4c_drv.generate(run_type.second, max_ang_mom, geom_drvs);
        }
    }
    
    if (run_type.first == "t4c_geom_cpu")
    {
        std::array<int, 4> geom_drvs = {1, 0, 1, 0};
        
        const auto t4c_geom_drv = T4CGeomDerivCPUGenerator();
            
        t4c_geom_drv.generate(max_ang_mom, geom_drvs);
    }
    
    if (run_type.first == "t4c_geom_hrr_cpu")
    {
        std::array<int, 4> geom_drvs = {1, 0, 1, 0};
        
        const auto t4c_geom_drv = T4CGeomHrrCPUGenerator();
            
        t4c_geom_drv.generate(run_type.second, max_ang_mom, geom_drvs);
    }
    
    if (run_type.first == "t2c_geom_cpu")
    {
        std::array<int, 3> geom_drvs = {0, 0, 1};
        
        const auto t2c_geom_drv = T2CGeomDerivCPUGenerator();
            
        t2c_geom_drv.generate(max_ang_mom, geom_drvs);
    }
        
    if (run_type.first == "t4c_diag_cpu")
    {
        const auto t4c_diag_drv = T4CDiagCPUGenerator();
            
        t4c_diag_drv.generate(run_type.second, max_ang_mom);
    }
   
    if (run_type.first == "t4c_call_tree")
    {
        const auto t4c_call_drv = T4CCallTreeGenerator();
            
        t4c_call_drv.generate(run_type.second, max_ang_mom);
    }
    
    // case: three-center integrals
    
    if (run_type.first == "t3c_cpu")
    {
        std::array<int, 3> geom_drvs = {1, 0, 0};
        
        if (geom_drvs == std::array<int, 3>({0, 0, 0}))
        {
            const auto t3c_drv = T3CCPUGenerator();
            
            t3c_drv.generate(run_type.second, max_ang_mom, max_ang_mom + 2);
        }
        else
        {
            const auto t3c_drv = T3CGeomCPUGenerator();
            
            t3c_drv.generate(run_type.second, max_ang_mom, max_ang_mom + 2, geom_drvs);
        }
    }
    
    if (run_type.first == "t3c_geom_hrr_cpu")
    {
        std::array<int, 3> geom_drvs = {1, 0, 0};
        
        const auto t3c_geom_drv = T3CGeomHrrCPUGenerator();
            
        t3c_geom_drv.generate(run_type.second, max_ang_mom + 2, geom_drvs);
    }
    
    // case: two-center integrals on grid

    if (run_type.first == "g2c_cpu")
    {
        std::array<int, 3> geom_drvs = {0, 0, 0};
        
        const auto use_rs = false;
        
        if ((geom_drvs[0] + geom_drvs[2]) == 0)
        {
            const auto g2c_drv = G2CCPUGenerator();
            
            g2c_drv.generate(run_type.second, max_ang_mom, geom_drvs, use_rs);
        }
        else
        {
//            const auto t2c_drv = T2CGeomCPUGenerator();
//            
//            t2c_drv.generate(run_type.second, max_ang_mom, geom_drvs, rec_form, use_rs);
        }
    }
    
    // case: two-center ECP integrals

    if (run_type.first == "t2c_ecp_cpu")
    {
        std::array<int, 3> geom_drvs = {0, 0, 0};
        
        if ((geom_drvs[0] + geom_drvs[2]) == 0)
        {
            const auto t2c_ecp_drv = T2CECPCPUGenerator();
            
            t2c_ecp_drv.generate(run_type.second, max_ang_mom);
        }
        else
        {
            //const auto t2c_drv = T2CGeomCPUGenerator();
            
            //t2c_drv.generate(run_type.second, max_ang_mom, geom_drvs, rec_form, use_rs);
        }
    }
    
    
    // set up end timer & compute elapsed time
    
    auto etime = std::chrono::high_resolution_clock::now();
    
    auto dtime = std::chrono::duration_cast<std::chrono::seconds>(etime - stime);
    
    std::cout << "Elapsed time: " << dtime.count() << " seconds." << std::endl;
    
    return 0;
}
