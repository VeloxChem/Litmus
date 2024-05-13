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
#include "t4c_geom_cpu_generators.hpp"

int main(int argc, char **argv)
{
    // run configuration

    // "Which kind of integral? # of centers", "which type of integral" (which operator is associated with it?)
    // 2c: "overlap" "kinetic energy" "nuclear potential" "dipole moment" "linear momentum"
    // 4c: "electron repulsion"

    const auto run_type = std::pair<std::string, std::string>({"t2c_cpu", "electric field"});
    //const auto run_type = std::pair<std::string, std::string>({"t2c_cpu", "nuclear potential"});

    const int max_ang_mom = 2;

    // To add new integral
    // (Be careful about scalar vs non-scalar integrals (see dipole for example of non-scalar)
    // 1. Add to is_available in t2(4)c_cpu_generators
    // 2. In get_integral (t2(4)c_cpu_generators), add appropriate entry (need operator name, rank of tensor representing operator rank if not scalar)
    // 3. In /recursions, make a driver file to house the new recursion that is needed and fill it with the appropriate recursion code
    // 4. Fill that file according to the example with comments for dipole integrals
        // In some cases may need to get fixed-axis operation (see linmom_driver vs dip_driver)
        // NB: In header file, may need to make changes like initalize to default (see difference like linmom_driver vs dip_driver)
        // In some cases (e.g. linmom, there is no choice of axes and one does not need a "trial recursion" setup from which to choose
    // 5. Make a v2(4)i driver file and fill it according to the example for dipole
    // 6. In t2(4)c_utils, register any new operator "name cases" (various places in this file)
    // 7. In t2(4)c_body, register any new operator "name cases" (up to three (two?) places in this file, marked) including variable declarations
    // 8. In t2(4)c_cpu_generators, register any new operator "name cases", (marked)
    // 9. In t2(4)c_prim_body, make any (marked) changes (remember includes)
    // 10. (see dipole for example) Add the primitive case to the "bottom" prim_ss autogen file (there will be blanks to fill in and you need to match this and the autogen func call that will call it)
    // 11. Modify the function declaration for the non-primitive calls in t2(4)c_decl (and (cosmetic) in t2(4)c_docs)

    // from t2c_cpu_generators.hpp:
    /// MR: First index: Summation or not (explicitly incorporate any multi-term nature in the operator associated with the integrals (if any))
    /// Second index: False: Return as matrix(ces); true: Return as scalars ("Contracted"/"distributed")
    const auto rec_form = std::pair<bool, bool>({true, false});
    
    // set up start timer
    
    auto stime = std::chrono::high_resolution_clock::now();
    
    // case: two-center integrals

    // Three-center also needed for RI generation; otherwise we get what we need for now with 2, 4

    if (run_type.first == "t2c_cpu")
    {
    // a, operator, b
        const std::array<int, 3> geom_drvs = {0, 2, 0};

        const auto t2c_drv = T2CCPUGenerator();
        
        t2c_drv.generate(run_type.second, max_ang_mom, geom_drvs, rec_form);
    }
    
    // case: four-center integrals
    
    if (run_type.first == "t4c_cpu")
    {
    // a, b, operator, c, d
        std::array<int, 5> geom_drvs = {1, 2, 0, 0, 2};
        
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
   
    // set up end timer & compute elapsed time
    
    auto etime = std::chrono::high_resolution_clock::now();
    
    auto dtime = std::chrono::duration_cast<std::chrono::seconds>(etime - stime);
    
    std::cout << "Elapsed time: " << dtime.count() << " seconds." << std::endl;
    
    return 0;
}
