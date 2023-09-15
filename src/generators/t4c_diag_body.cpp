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

#include "t4c_diag_body.hpp"

#include "string_formater.hpp"
#include "spherical_momentum.hpp"
#include "t2c_utils.hpp"

void
T4CDiagFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                       const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    //for (const auto& label : _get_angmom_def(integral))
    //{
    //    lines.push_back({1, 0, 2, label});
    //}
    
    //for (const auto& label : _get_gtos_def(diagonal))
    //{
    //    lines.push_back({1, 0, 2, label});
    //}
    
    //for (const auto& label : _get_ket_variables_def())
    //{
    //    lines.push_back({1, 0, 2, label});
    //}
    
    //for (const auto& label : _get_buffers_def(integral))
    //{
    //    lines.push_back({1, 0, 2, label});
    //}
    
    //for (const auto& label : _get_batches_def(diagonal))
    //{
    //    lines.push_back({1, 0, 2, label});
    //}
    
    //_add_batches_loop_start(lines);
    
    //_add_batches_loop_body(lines, diagonal);
    
    //_add_bra_loop_start(lines, diagonal);
   
    //_add_bra_loop_body(lines, integral, diagonal);
    
    //_add_bra_loop_end(lines);
    
    //_add_batches_loop_end(lines);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}
