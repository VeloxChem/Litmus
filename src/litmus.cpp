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

    // "Which kind of integral? # of centers", "which type of integral" (which operator is associated with it?)
    // Keywords are listed in
    // 2c: "overlap" "kinetic energy" "nuclear potential" "dipole moment"
    // 4c: "electron repulsion"

    const auto run_type = std::pair<std::string, std::string>({"t4c_cpu", "electron repulsion"});

    const int max_ang_mom = 4;

    // To add new integral
    // (Be careful about scalar vs non-scalar integrals (see dipole for example of non-scalar)
    // 1. Add to is_available
    // 2. In get_integral, add appropriate entry (need operator name, rank of tensor representing operator rank if not scalar)
    // 3. In /recursions, make a driver file to house the new recursion that is needed and fill it with the appropriate recursion code
    // 4. Fill that file according to the example with comments for dipole integrals
    // 5. Make a v2i driver file or eqv for 4-center and fill it according to the according example for dipole
    // 6. In t2(4)c_utils, register any new operator "name cases" (various places in this file)
    // 7. In t2(4)c_body, register any new operator "name cases" (up to three (two?) places in this file, marked)
    // 8. In t2(4)c_cpu_generators, register any new operator "name cases", (marked)
    // 9. In t2(4)c_prim_body, make any (marked) changes (remember includes)

    // from t2c_cpu_generators.hpp:
    /// MR: First index: Summation or not (explicitly incorporate any multi-term nature in the operator associated with the integrals (if any))
    /// Second index: False: Return as matrix(ces); true: Return as scalars ("Contracted"/"distributed")
    const auto rec_form = std::pair<bool, bool>({false, false});
    
    // set up start timer
    
    auto stime = std::chrono::high_resolution_clock::now();
    
    // case: two-center integrals

    // Three-center also needed for RI generation; otherwise we get what we need for now with 2, 4

    if (run_type.first == "t2c_cpu")
    {
    // a, operator, b
        const std::array<int, 3> geom_drvs = {0, 0, 0};


        const auto t2c_drv = T2CCPUGenerator();
        
        t2c_drv.generate(run_type.second, max_ang_mom, geom_drvs, rec_form);
    }
    
    // case: four-center integrals
    
    if (run_type.first == "t4c_cpu")
    {
    // a, b, operator, c, d
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
