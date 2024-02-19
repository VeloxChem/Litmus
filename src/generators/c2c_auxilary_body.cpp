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

#include "c2c_auxilary_body.hpp"

void
C2CAuxilaryBodyDriver::write_aux_body(      std::ofstream& fstream,
                                      const R2Group&       rgroup,
                                      const I2CIntegral&   integral,
                                      const bool           diagonal) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_math_def())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_gtos_def(diagonal))
    {
        lines.push_back({1, 0, 2, label});
    }

    for (const auto& label : _get_ket_variables_def())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_auxilaries_def(rgroup))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_bra_coords(diagonal))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_ket_coords(diagonal))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_ket_pointers_def())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_prim_loop_start(lines, diagonal);
    
    _add_aux_loop_body(lines, rgroup, integral); 
    
    _add_prim_loop_end(lines);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
C2CAuxilaryBodyDriver::_get_gtos_def(const bool diagonal) const
{
    std::vector<std::string> vstr;
    
    if (diagonal)
    {
        vstr.push_back("// intialize GTOs data");
        
        vstr.push_back("const auto gto_coords = gto_block.getCoordinates();");
       
        vstr.push_back("const auto gto_exps = gto_block.getExponents();");
       
        vstr.push_back("const auto gto_norms = gto_block.getNormalizationFactors();");
        
        vstr.push_back("const auto gto_indexes = gto_block.getOrbitalIndexes();");
       
        vstr.push_back("const auto ncgtos = gto_block.getNumberOfBasisFunctions();");

        vstr.push_back("const auto npgtos = gto_block.getNumberOfPrimitives();");
    }
    else
    {
        vstr.push_back("// intialize GTOs data on bra side");
        
        vstr.push_back("const auto bra_gto_coords = bra_gto_block.getCoordinates();");
       
        vstr.push_back("const auto bra_gto_exps = bra_gto_block.getExponents();");
       
        vstr.push_back("const auto bra_gto_norms = bra_gto_block.getNormalizationFactors();");
        
        vstr.push_back("const auto bra_gto_indexes = bra_gto_block.getOrbitalIndexes();");
       
        vstr.push_back("const auto bra_ncgtos = bra_gto_block.getNumberOfBasisFunctions();");

        vstr.push_back("const auto bra_npgtos = bra_gto_block.getNumberOfPrimitives();");
        
        vstr.push_back("// intialize GTOs data on ket side");
        
        vstr.push_back("const auto ket_gto_coords = ket_gto_block.getCoordinates();");
       
        vstr.push_back("const auto ket_gto_exps = ket_gto_block.getExponents();");
       
        vstr.push_back("const auto ket_gto_norms = ket_gto_block.getNormalizationFactors();");
        
        vstr.push_back("const auto ket_gto_indexes = ket_gto_block.getOrbitalIndexes();");
       
        vstr.push_back("const auto ket_ncgtos = ket_gto_block.getNumberOfBasisFunctions();");

        vstr.push_back("const auto ket_npgtos = ket_gto_block.getNumberOfPrimitives();");
    }
    
    return vstr;
}

std::vector<std::string>
C2CAuxilaryBodyDriver::_get_ket_variables_def() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// initialize aligned arrays for ket side");
        
    vstr.push_back("alignas(64) TArray<double> ket_coords_x;");
        
    vstr.push_back("alignas(64) TArray<double> ket_coords_y;");
        
    vstr.push_back("alignas(64) TArray<double> ket_coords_z;");
        
    vstr.push_back("alignas(64) TArray<double> ket_exps;");
        
    vstr.push_back("alignas(64) TArray<double> ket_norms;");
        
    return vstr;
}

std::vector<std::string>
C2CAuxilaryBodyDriver::_get_auxilaries_def(const R2Group& rgroup) const
{
    std::vector<std::string> vstr;

    auto ndims = (t2c::get_unique_auxilaries(rgroup)).size();
    
    if (ndims == 0) ndims = 1;

    vstr.push_back("// set up pointers to auxilary buffers");

    for (size_t i = 0; i < ndims; i++)
    {
        vstr.push_back("auto avals_" + std::to_string(i) + " = auxilaries[" + std::to_string(i) + "].data();");
    }
    
    vstr.push_back("// zero auxilary buffers");
    
    vstr.push_back("simd::zero(auxilaries);");

    return vstr;
}

std::vector<std::string>
C2CAuxilaryBodyDriver::_get_bra_coords(const bool diagonal) const
{
    std::vector<std::string> vstr;

    vstr.push_back("// set up coordinates on bra side");
    
    if (diagonal)
    {
        vstr.push_back("const auto a_x = gto_coords[bra_index][0];");
        
        vstr.push_back("const auto a_y = gto_coords[bra_index][1];");
        
        vstr.push_back("const auto a_z = gto_coords[bra_index][2];");
    }
    else
    {
        vstr.push_back("const auto a_x = bra_gto_coords[bra_index][0];");
        
        vstr.push_back("const auto a_y = bra_gto_coords[bra_index][1];");
        
        vstr.push_back("const auto a_z = bra_gto_coords[bra_index][2];");
    }
    
    return vstr;
}

std::vector<std::string>
C2CAuxilaryBodyDriver::_get_ket_coords(const bool diagonal) const
{
    std::vector<std::string> vstr;

    vstr.push_back("// set up coordinates on ket side");
    
    std::string label = "simd::loadCoordinates(ket_coords_x, ket_coords_y, ket_coords_z, ";
    
    if (diagonal)
    {
        label +=  "gto_coords,";
    }
    else
    {
        label += "ket_gto_coords,";
    }
    
    label += " ket_first, ket_last);";
    
    vstr.push_back(label);
        
    vstr.push_back("auto ket_rx = ket_coords_x.data();");
        
    vstr.push_back("auto ket_ry = ket_coords_y.data();");
        
    vstr.push_back("auto ket_rz = ket_coords_z.data();");
        
    return vstr;
}

std::vector<std::string>
C2CAuxilaryBodyDriver::_get_ket_pointers_def() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// set up pointers to GTOs data for ket side");
        
    vstr.push_back("auto ket_fe = ket_exps.data();");
        
    vstr.push_back("auto ket_fn = ket_norms.data();");
    
    vstr.push_back("// set up ket dimensions");
    
    vstr.push_back("const auto ket_dim = ket_last - ket_first;");
    
    return vstr;
}

void
C2CAuxilaryBodyDriver::_add_prim_loop_start(      VCodeLines& lines,
                                            const bool        diagonal) const
{
    lines.push_back({1, 0, 2, "// compute auxilary integrals"});
    
    if (diagonal)
    {
        lines.push_back({1, 0, 1, "for (int i = 0; i < npgtos; i++)"});
            
        lines.push_back({1, 0, 1, "{"});
            
        lines.push_back({2, 0, 2, "simd::loadPrimitiveGTOsData(ket_exps, gto_exps, i, ncgtos, ket_first, ket_last);"});
                
        lines.push_back({2, 0, 2, "simd::loadPrimitiveGTOsData(ket_norms, gto_norms, i, ncgtos, ket_first, ket_last);"});
                
        lines.push_back({2, 0, 1, "for (int j = 0; j < npgtos; j++)"});
            
        lines.push_back({2, 0, 1, "{"});
            
        lines.push_back({3, 0, 2, "const auto bra_idx = j * ncgtos + bra_index;"});
                    
        lines.push_back({3, 0, 2, "const auto bra_exp = gto_exps[bra_idx];"});
                    
        lines.push_back({3, 0, 2, "const auto bra_norm = gto_norms[bra_idx];"});
    }
    else
    {
        lines.push_back({1, 0, 1, "for (int i = 0; i < ket_npgtos; i++)"});
            
        lines.push_back({1, 0, 1, "{"});
            
        lines.push_back({2, 0, 2, "simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, i, ket_ncgtos, ket_first, ket_last);"});
                
        lines.push_back({2, 0, 2, "simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, i, ket_ncgtos, ket_first, ket_last);"});
                
        lines.push_back({2, 0, 1, "for (int j = 0; j < bra_npgtos; j++)"});
            
        lines.push_back({2, 0, 1, "{"});
            
        lines.push_back({3, 0, 2, "const auto bra_idx = j * bra_ncgtos + bra_index;"});
                    
        lines.push_back({3, 0, 2, "const auto bra_exp = bra_gto_exps[bra_idx];"});
                    
        lines.push_back({3, 0, 2, "const auto bra_norm = bra_gto_norms[bra_idx];"});
    }
}

void
C2CAuxilaryBodyDriver::_add_aux_loop_body(      VCodeLines&  lines,
                                          const R2Group&     rgroup,
                                          const I2CIntegral& integral) const
{
    lines.push_back({3, 0, 1, "#pragma omp simd aligned(ket_rx, ket_ry, ket_rz, ket_fe, ket_fn : 64)"});
    
    lines.push_back({3, 0, 1, "for (int k = 0; k < ket_dim; k++)"});
    
    lines.push_back({3, 0, 1, "{"});
    
    lines.push_back({4, 0, 2, "const auto ab_x = bra_rx - ket_rx[k];"});
    
    lines.push_back({4, 0, 2, "const auto ab_y = bra_ry - ket_ry[k];"});
    
    lines.push_back({4, 0, 2, "const auto ab_z = bra_rz - ket_rz[k];"});
    
    // determine auxilaries
    
    const auto auxilaries = t2c::get_unique_auxilaries(rgroup);
    
    _add_aux_overlap_factor(lines, auxilaries);
    
    _add_aux_polynomial_factors(lines, auxilaries);
    
    _add_aux_values(lines, auxilaries);
    
    lines.push_back({3, 0, 1, "}"});
}

void
C2CAuxilaryBodyDriver::_add_prim_loop_end(VCodeLines& lines) const
{
    lines.push_back({2, 0, 1, "}"});
    
    lines.push_back({1, 0, 2, "}"});
}

void
C2CAuxilaryBodyDriver::_add_aux_overlap_factor(      VCodeLines&   lines,
                                               const V4Auxilaries& auxilaries) const
{
    lines.push_back({4, 0, 2, "const auto ket_exp = ket_fe[k];"});
    
    lines.push_back({4, 0, 2, "const auto fe_0 = 1.0 / (bra_exp + ket_exp);"});
    
    lines.push_back({4, 0, 2, "const auto fz_0 = bra_exp * ket_exp * fe_0 * (ab_x * ab_x + ab_y * ab_y + ab_z * ab_z);"});
    
    lines.push_back({4, 0, 2, "const auto fmpi = fpi * fe_0;"});
    
    if (const auto ndims = auxilaries.size(); ndims == 0)
    {
        lines.push_back({4, 0, 2, "avals_0[k] += bra_norm * ket_fn[k] * fmpi * std::sqrt(fmpi) * std::exp(-fz_0);"});
    }
    else
    {
        lines.push_back({4, 0, 2, "const auto fss = bra_norm * ket_fn[k] * fmpi * std::sqrt(fmpi) * std::exp(-fz_0);"});
    }
}

void
C2CAuxilaryBodyDriver::_add_aux_polynomial_factors(      VCodeLines&   lines,
                                                   const V4Auxilaries& auxilaries) const
{
    const auto mvals = t2c::get_maximum_decomposition(auxilaries);
    
    if (mvals[0] > 0)
    {
        lines.push_back({4, 0, 2, "const auto ft_0 = bra_exp * ket_exp * fe_0;"});
    }
    
    if (mvals[1] > 0)
    {
        lines.push_back({4, 0, 2, "const auto fm_0 = bra_exp * fe_0;"});
    }
    
    if (mvals[2] > 0)
    {
        lines.push_back({4, 0, 2, "const auto fn_0 = ket_exp * fe_0;"});
    }
}

void
C2CAuxilaryBodyDriver::_add_aux_values(      VCodeLines&   lines,
                                       const V4Auxilaries& auxilaries) const
{
    int index = 0;
    
    for (const auto& taux : auxilaries)
    {
        auto label = "avals_" + std::to_string(index) + "[k] += fss";
        
        const auto mvals = t2c::get_factor_decomposition(taux);
        
        for (int i = 0; i < mvals[0]; i++) label += " * ft_0";
        
        for (int i = 0; i < mvals[1]; i++) label += " * fm_0";
        
        for (int i = 0; i < mvals[2]; i++) label += " * fn_0";
        
        for (int i = 0; i < (taux[2] - mvals[0] - mvals[1] - mvals[2]); i++) label += " * fe_0";
        
        for (int i = 0; i < (taux[0] - mvals[0] - mvals[1]); i++) label += " * bra_exp";
        
        for (int i = 0; i < (taux[1] - mvals[0] - mvals[2]); i++) label += " * ket_exp";
        
        label += ";";
        
        lines.push_back({4, 0, 2, label});
        
        index++;
    }
}

std::vector<std::string>
C2CAuxilaryBodyDriver::_get_math_def() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// set up math constants");
        
    vstr.push_back("const auto fpi = mathconst::getPiValue();");
 
    return vstr;
}

