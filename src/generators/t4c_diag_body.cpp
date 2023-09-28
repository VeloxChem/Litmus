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
#include "t4c_utils.hpp"

void
T4CDiagFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                       const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_gtos_def())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_vars_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_batches_def())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_batches_loop_start(lines);
    
    _add_batches_loop_body(lines, integral);
    
    _add_batches_loop_end(lines);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CDiagFuncBodyDriver::_get_gtos_def() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// intialize GTO pairs data");

    vstr.push_back("const auto gpair_coords = gto_pair_block.getCoordinates();");

    vstr.push_back("const auto gpair_exps = gto_pair_block.getExponents();");

    vstr.push_back("const auto gpair_norms = gto_pair_block.getNormalizationFactors();");

    vstr.push_back("const auto nppairs = gto_pair_block.getNumberOfPrimitivePairs();");
    
    vstr.push_back("const auto ncpairs = gto_pair_block.getNumberOfContractedPairs();");
    
    return vstr;
}

std::vector<std::string>
T4CDiagFuncBodyDriver::_get_vars_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// set up maximum Cartesian integrals vector");

    vstr.push_back("std::vector<double> max_tints(ncpairs, 0.0);");

    vstr.push_back("// initialize aligned arrays for A and B centers");

    vstr.push_back("alignas(64) TDoubleArray coords_a_x;");

    vstr.push_back("alignas(64) TDoubleArray coords_a_y;");

    vstr.push_back("alignas(64) TDoubleArray coords_a_z;");

    vstr.push_back("alignas(64) TDoubleArray coords_b_x;");

    vstr.push_back("alignas(64) TDoubleArray coords_b_y;");

    vstr.push_back("alignas(64) TDoubleArray coords_b_z;");

    vstr.push_back("// initialize aligned arrays for bra side");

    vstr.push_back("alignas(64) TDoubleArray bra_exps_a;");

    vstr.push_back("alignas(64) TDoubleArray bra_exps_b;");

    vstr.push_back("alignas(64) TDoubleArray bra_norms;");

    vstr.push_back("// initialize aligned arrays for ket side");

    vstr.push_back("alignas(64) TDoubleArray ket_exps_c;");

    vstr.push_back("alignas(64) TDoubleArray ket_exps_d;");

    vstr.push_back("alignas(64) TDoubleArray ket_norms;");

    vstr.push_back("// initialize contracted integrals buffer");

    vstr.push_back("alignas(64) TDoubleArray buffer;");
    
    if ((integral[0] + integral[1] + integral[2] + integral[3]) > 0)
    {
        vstr.push_back("alignas(64) TDoubleArray max_buffer;");
    }
    
    return vstr;
}

std::vector<std::string>
T4CDiagFuncBodyDriver::_get_batches_def() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// loop over integral batches");
    
    vstr.push_back("const auto nbatches = batch::getNumberOfBatches(ncpairs, simd_width);");
       
    return vstr;
}

void
T4CDiagFuncBodyDriver::_add_batches_loop_start(VCodeLines& lines) const
{
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < nbatches; i++)"});
        
    lines.push_back({1, 0, 1, "{"});
}

void
T4CDiagFuncBodyDriver::_add_batches_loop_body(      VCodeLines&  lines,
                                              const I4CIntegral& integral) const
{
    lines.push_back({2, 0, 2, "const auto [first, last] = batch::getBatchRange(i, ncpairs, simd_width);"});
    
    lines.push_back({2, 0, 2, "const auto ndim = last - first;"});

    lines.push_back({2, 0, 2, "// load coordinates data"});
    
    lines.push_back({2, 0, 2, "simd::loadCoordinates(coords_a_x, coords_a_y, coords_a_z, coords_b_x, coords_b_y, coords_b_z, gpair_coords, first, last);"});
    
    for (const auto& tcomp : integral.diag_components<T2CPair, T2CPair>())
    {
        _add_component_body(lines, integral, tcomp);
    }
}

void
T4CDiagFuncBodyDriver::_add_batches_loop_end(VCodeLines& lines) const
{
    lines.push_back({2, 0, 1, "t4cfunc::distribute(max_tints, max_buffer, first, last);"});
    
    lines.push_back({1, 0, 2, "}"});
    
    lines.push_back({1, 0, 1, "return max_tints;"});
}

void
T4CDiagFuncBodyDriver::_add_component_body(      VCodeLines&  lines,
                                           const I4CIntegral& integral, 
                                           const T4CIntegral& component) const
{
    auto [nsize, name] = t4c::prim_diag_compute_func_name(component, integral);
    
    name = t4c::namespace_label(integral) + "::" + name;
    
    lines.push_back({2, 0, 2, "// compute primitive integrals block (" + fstr::upcase(component.label()) + ")"});
    
    lines.push_back({2, 0, 2, "simd::zero(buffer);"});
    
    lines.push_back({2, 0, 1, "for (int64_t j = 0; j < nppairs;  j++)"});
    
    lines.push_back({2, 0, 1, "{"});
        
    lines.push_back({3, 0, 2, "simd::loadPrimitiveGTOsData(bra_norms, gpair_norms, j, ncpairs, first, last);"});
        
    lines.push_back({3, 0, 2, "simd::loadPrimitiveGTOsPairsData(bra_exps_a, bra_exps_b, gpair_exps, j, ncpairs, first, last);"});
        
    lines.push_back({3, 0, 2, name + "(buffer, coords_a_x, coords_a_y, coords_a_z, coords_b_x, coords_b_y, coords_b_z, bra_exps_a, bra_exps_b, bra_norms, ndim);"});
    
    lines.push_back({3, 0, 1, "for (int64_t k = j + 1; k < nppairs; k++)"});
    
    lines.push_back({3, 0, 1, "{"});
    
    lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsData(ket_norms, gpair_norms, k, ncpairs, first, last);"});
        
    lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsPairsData(ket_exps_c, ket_exps_d, gpair_exps, k, ncpairs, first, last);"});
    
    lines.push_back({4, 0, 1, name + "(buffer, coords_a_x, coords_a_y, coords_a_z, coords_b_x, coords_b_y, coords_b_z, bra_exps_a, bra_exps_b, bra_norms, ket_exps_c, ket_exps_d, ket_norms, ndim);"});
    
    lines.push_back({3, 0, 1, "}"});
    
    lines.push_back({2, 0, 2, "}"});
    
    lines.push_back({2, 0, 2, "simd::max_update(max_buffer, buffer);"});
}


