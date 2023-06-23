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

#include "t2c_body.hpp"

#include "string_formater.hpp"
#include "spherical_momentum.hpp"
#include "t2c_utils.hpp"

void
T2CFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                   const I2CIntegral&   integral,
                                   const bool           diagonal) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_angmom_def(integral))
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
    
    for (const auto& label : _get_buffers_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_batches_def(diagonal))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_batches_loop_start(lines);
    
    _add_batches_loop_body(lines, diagonal); 
    
    _add_bra_loop_start(lines, diagonal);
   
    _add_bra_loop_body(lines, integral, diagonal);
    
    _add_bra_loop_end(lines);
    
    _add_batches_loop_end(lines);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CFuncBodyDriver::_get_angmom_def(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if ((integral[0] > 1) || (integral[1] > 1))
    {
        const auto angmom = SphericalMomentum(0);
            
        vstr.push_back("// spherical transformation factors");
        
        if (integral[0] > 1)
        {
            for (const auto& label : angmom.get_factors(integral[0]))
            {
                 vstr.push_back("const double " + label + ";");
            }
        }
        
        if ((integral[1] > 1) && (integral[0] != integral[1]))
        {
            for (const auto& label : angmom.get_factors(integral[1]))
            {
                vstr.push_back("const double " + label + ";");
            }
        }
    }
    
    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_gtos_def(const bool diagonal) const
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
T2CFuncBodyDriver::_get_ket_variables_def() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// initialize aligned arrays for ket side");
        
    vstr.push_back("alignas(64) TDoubleArray ket_coords_x;");
        
    vstr.push_back("alignas(64) TDoubleArray ket_coords_y;");
        
    vstr.push_back("alignas(64) TDoubleArray ket_coords_z;");
        
    vstr.push_back("alignas(64) TDoubleArray ket_exps;");
        
    vstr.push_back("alignas(64) TDoubleArray ket_norms;");
        
    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_buffers_def(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// initialize contracted integrals buffer");
    
    std::vector<std::string> labels;
    
    if (integral.is_simple() && integral.is_simple_integrand())
    {
        const auto bra = Tensor(integral[0]);
        
        const auto ket = Tensor(integral[1]);
        
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            labels = {"buffer", }; 
            
            if (integral[0] > 0) labels = t2c::tensor_components(bra, "buffer");
            
            if (integral[1] > 0) labels = t2c::tensor_components(ket, "buffer");
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                labels = t2c::tensor_components(ket, "buffer");
            }
            else
            {
                labels = t2c::tensor_components(bra, "buffer");
            }
        }
    }
    else
    {
        labels = t2c::integrand_components(integral.integrand(), "buffer");
    }
    
    for (const auto& label : labels)
    {
        vstr.push_back("alignas(64) TDoubleArray " + label + ";");
    }
    
    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_batches_def(const bool diagonal) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// loop over integral batches");
    
    if (diagonal)
    {
        vstr.push_back("const auto nbatches = batch::getNumberOfBatches(ncgtos, simd_width);");
    }
    else
    {
        vstr.push_back("const auto nbatches = batch::getNumberOfBatches(ket_ncgtos, simd_width);");
    }
        
    return vstr;
}

void
T2CFuncBodyDriver::_add_batches_loop_start(VCodeLines& lines) const
{
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < nbatches; i++)"});
        
    lines.push_back({1, 0, 1, "{"});
}

void 
T2CFuncBodyDriver::_add_batches_loop_body(      VCodeLines& lines,
                                          const bool        diagonal) const
{
    if (diagonal)
    {
        lines.push_back({2, 0, 2, "const auto [ket_first, ket_last] = batch::getBatchRange(i, ncgtos, simd_width);"});
    }
    else
    {
        lines.push_back({2, 0, 2, "const auto [ket_first, ket_last] = batch::getBatchRange(i, ket_ncgtos, simd_width);"});
    }
            
    lines.push_back({2, 0, 2, "const auto ket_dim = ket_last - ket_first;"});
            
    lines.push_back({2, 0, 1, "simd::loadCoordinates(ket_coords_x,"});
        
    lines.push_back({2, 22, 1, "ket_coords_y,"});
        
    lines.push_back({2, 22, 1, "ket_coords_z,"});

    if (diagonal)
    {
        lines.push_back({2, 22, 1, "gto_coords,"});
    }
    else
    {
        lines.push_back({2, 22, 1, "ket_gto_coords,"});
    }
        
    lines.push_back({2, 22, 1, "ket_first,"});
        
    lines.push_back({2, 22, 2, "ket_last);"});
}

void
T2CFuncBodyDriver::_add_batches_loop_end(VCodeLines& lines) const
{
    lines.push_back({1, 0, 1, "}"});
}

void
T2CFuncBodyDriver::_add_bra_loop_start(      VCodeLines& lines,
                                       const bool        diagonal) const
{
    lines.push_back({2, 0, 1, "for (int64_t j = bra_first; j < bra_last; j++) "});
        
    lines.push_back({2, 0, 1, "{"});
        
    if (diagonal)
    {
        lines.push_back({3, 0, 2, "const auto bra_coord = gto_coords[j];"});
    }
    else
    {
        lines.push_back({3, 0, 2, "const auto bra_coord = bra_gto_coords[j];"});
    }
}

void
T2CFuncBodyDriver::_add_bra_loop_end(VCodeLines& lines) const
{
    lines.push_back({2, 0, 1, "}"});
}

void
T2CFuncBodyDriver::_add_bra_loop_body(      VCodeLines&  lines,
                                      const I2CIntegral& integral,
                                      const bool         diagonal) const
{
    if (integral.is_simple() && integral.is_simple_integrand())
    {
        if ((integral[0] == 0) || (integral[1] == 0))
        {
            _add_loop_call_tree(lines, integral, diagonal);
        }
        else
        {
            if (integral[0] >= integral[1])
            {
                for (const auto& bcomp: Tensor(integral[0]).components())
                {
                    _add_loop_call_tree(lines, bcomp, integral, true, diagonal);
                }
            }
            else
            {
                for (const auto& kcomp: Tensor(integral[1]).components())
                {
                    _add_loop_call_tree(lines, kcomp, integral, false, diagonal);
                }
            }
        }
    }
    else
    {
        for (const auto& bcomp: Tensor(integral[0]).components())
        {
            for (const auto& kcomp: Tensor(integral[1]).components())
            {
                _add_loop_call_tree(lines, bcomp, kcomp, integral, diagonal);
            }
        }
    }
}

void
T2CFuncBodyDriver::_add_loop_call_tree(      VCodeLines&  lines,
                                       const I2CIntegral& integral,
                                       const bool         diagonal) const
{
    std::vector<std::string> labels({"buffer", });
    
    if (integral[0] > 0) labels = t2c::tensor_components(Tensor(integral[0]), "buffer");
    
    if (integral[1] > 0) labels = t2c::tensor_components(Tensor(integral[1]), "buffer");
    
    lines.push_back({3, 0, 2, "// compute primitive integrals block"});
    
    for (const auto& label : labels)
    {
        lines.push_back({3, 0, 2, "simd::zero(" + label + ");"});
    }
    
    _add_prim_loop_start(lines, diagonal);
    
    auto [nsize, name] = t2c::prim_compute_func_name(integral);
    
    name = t2c::namespace_label(integral) + "::" + name;
    
    nsize = name.size() + 1 ;
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        if (i == 0)
        {
            lines.push_back({5, 0, 1, name + "(" + labels[i] + ","});
        }
        else
        {
            lines.push_back({5, nsize, 1, labels[i] + ","});
        }
    }
    
    _add_prim_call_data(lines, nsize); 
    
    _add_prim_loop_end(lines);
    
    _write_block_distributor(lines, integral, diagonal); 
}

void
T2CFuncBodyDriver::_add_loop_call_tree(      VCodeLines&      lines,
                                       const TensorComponent& component,
                                       const I2CIntegral&     integral,
                                       const bool             bra_first,
                                       const bool             diagonal) const
{
    const auto labels = (bra_first) ? t2c::tensor_components(Tensor(integral[1]), "buffer")
                                    : t2c::tensor_components(Tensor(integral[0]), "buffer");
    
    lines.push_back({3, 0, 2, "// compute primitive integrals block (" + fstr::upcase(component.label()) + ")"});
    
    for (const auto& label : labels)
    {
        lines.push_back({3, 0, 2, "simd::zero(" + label + ");"});
    }
    
    _add_prim_loop_start(lines, diagonal);
    
    auto [nsize, name] = t2c::prim_compute_func_name(component, integral, bra_first);
    
    name = t2c::namespace_label(integral) + "::" + name;
    
    nsize = name.size() + 1 ;
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        if (i == 0)
        {
            lines.push_back({5, 0, 1, name + "(" + labels[i] + ","});
        }
        else
        {
            lines.push_back({5, nsize, 1, labels[i] + ","});
        }
    }
    
    _add_prim_call_data(lines, nsize);
    
    _add_prim_loop_end(lines);
    
    _write_block_distributor(lines, component, integral, bra_first, diagonal);
}

void
T2CFuncBodyDriver::_add_loop_call_tree(      VCodeLines&      lines,
                                       const TensorComponent& bra_component,
                                       const TensorComponent& ket_component,
                                       const I2CIntegral&     integral,
                                       const bool             diagonal) const
{
    const auto labels = t2c::integrand_components(integral.integrand(), "buffer");
    
    lines.push_back({3, 0, 2, "// compute primitive integrals block (" +
                   
                    fstr::upcase(bra_component.label()) + "_" +
                   
                    fstr::upcase(ket_component.label()) + ")"});
    
    for (const auto& label : labels)
    {
        lines.push_back({3, 0, 2, "simd::zero(" + label + ");"});
    }
    
    _add_prim_loop_start(lines, diagonal);
    
    auto [nsize, name] = t2c::prim_compute_func_name(bra_component, ket_component, integral);
    
    name = t2c::namespace_label(integral) + "::" + name;
    
    nsize = name.size() + 1 ;
    
    for (size_t i = 0; i < labels.size(); i++)
    {
        if (i == 0)
        {
            lines.push_back({5, 0, 1, name + "(" + labels[i] + ","});
        }
        else
        {
            lines.push_back({5, nsize, 1, labels[i] + ","});
        }
    }
    
    _add_prim_call_data(lines, nsize);
    
    _add_prim_loop_end(lines);
    
    _write_block_distributor(lines, bra_component, ket_component, integral, diagonal);
}

void
T2CFuncBodyDriver::_add_prim_loop_start(      VCodeLines& lines,
                                        const bool        diagonal) const
{
    if (diagonal)
    {
        lines.push_back({3, 0, 1, "for (int64_t k = 0; k < npgtos; k++)"});
            
        lines.push_back({3, 0, 1, "{"});
            
        lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsData(ket_exps, gto_exps, k, ncgtos, ket_first, ket_last);"});
                
        lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsData(ket_norms, gto_norms, k, ncgtos, ket_first, ket_last);"});
                
        lines.push_back({4, 0, 1, "for (int64_t l = 0; l < npgtos; l++)"});
            
        lines.push_back({4, 0, 1, "{"});
            
        lines.push_back({5, 0, 2, "const auto bra_index = l * ncgtos + j;"});
                    
        lines.push_back({5, 0, 2, "const auto bra_exp = gto_exps[bra_index];"});
                    
        lines.push_back({5, 0, 2, "const auto bra_norm = gto_norms[bra_index];"});
    }
    else
    {
        lines.push_back({3, 0, 1, "for (int64_t k = 0; k < ket_npgtos; k++)"});
            
        lines.push_back({3, 0, 1, "{"});
            
        lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsData(ket_exps, ket_gto_exps, k, ket_ncgtos, ket_first, ket_last);"});
                
        lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsData(ket_norms, ket_gto_norms, k, ket_ncgtos, ket_first, ket_last);"});
                
        lines.push_back({4, 0, 1, "for (int64_t l = 0; l < bra_npgtos; l++)"});
            
        lines.push_back({4, 0, 1, "{"});
            
        lines.push_back({5, 0, 2, "const auto bra_index = l * bra_ncgtos + j;"});
                    
        lines.push_back({5, 0, 2, "const auto bra_exp = bra_gto_exps[bra_index];"});
                    
        lines.push_back({5, 0, 2, "const auto bra_norm = bra_gto_norms[bra_index];"});
    }
}

void
T2CFuncBodyDriver::_add_prim_loop_end(VCodeLines& lines) const
{
    lines.push_back({4, 0, 1, "}"});
    
    lines.push_back({3, 0, 2, "}"});
}

void
T2CFuncBodyDriver::_add_prim_call_data(      VCodeLines& lines,
                                       const size_t      spacer) const
{
    lines.push_back({5, spacer, 1, "bra_exp,"});
    
    lines.push_back({5, spacer, 1, "bra_norm,"});
    
    lines.push_back({5, spacer, 1, "bra_coord,"});
    
    lines.push_back({5, spacer, 1, "ket_exps,"});
    
    lines.push_back({5, spacer, 1, "ket_norms,"});
    
    lines.push_back({5, spacer, 1, "ket_coords_x,"});
    
    lines.push_back({5, spacer, 1, "ket_coords_y,"});
    
    lines.push_back({5, spacer, 1, "ket_coords_z,"});
    
    lines.push_back({5, spacer, 1, "ket_dim);"});
}

void
T2CFuncBodyDriver::_write_block_distributor(      VCodeLines&  lines,
                                            const I2CIntegral& integral,
                                            const bool         diagonal) const
{
    if ((integral[0] + integral[1]) == 0)
    {
        if (diagonal)
        {
            lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, buffer, gto_indexes,"});
                            
            lines.push_back({3, 20, 1, "0, 0, j, ket_first, ket_last);"});
        }
        else
        {
            lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, buffer, bra_gto_indexes, ket_gto_indexes,"});
            
            lines.push_back({3, 20, 1, "0, 0, j, ket_first, ket_last, mat_type);"});
        }
    }
    
    if (integral[0] > 0)
    {
        const auto bra_mom = SphericalMomentum(integral[0]);
        
        const auto labels = t2c::tensor_components(Tensor(integral[0]), "buffer");
        
        for (int i = 0; i < labels.size(); i++)
        {
            for (const auto& pair : bra_mom.select_pairs(i))
            {
                auto flabel = labels[i];
                
                flabel += (pair.second == "1.0")  ? "" : ", "  + pair.second;
                
                const auto findex = std::to_string(pair.first);
                  
                lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, " + flabel +
                                          ", bra_gto_indexes, ket_gto_indexes,"});
                    
                lines.push_back({3, 20, 2, findex + ", 0, j, ket_first, ket_last, ang_order);"});
            }
        }
    }
    
    if (integral[1] > 0)
    {
        const auto ket_mom = SphericalMomentum(integral[1]);
        
        const auto labels = t2c::tensor_components(Tensor(integral[1]), "buffer");
        
        for (int i = 0; i < labels.size(); i++)
        {
            for (const auto& pair : ket_mom.select_pairs(i))
            {
                auto flabel = labels[i];
                
                flabel += (pair.second == "1.0")  ? "" : ", "  + pair.second;
                
                const auto findex = std::to_string(pair.first);
              
                lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, " + flabel +
                                          ", bra_gto_indexes, ket_gto_indexes,"});
                    
                lines.push_back({3, 20, 2, "0, " + findex + ", j, ket_first, ket_last, ang_order);"});
            }
        }
    }
}

void
T2CFuncBodyDriver::_write_block_distributor(      VCodeLines&      lines,
                                            const TensorComponent& component,
                                            const I2CIntegral&     integral,
                                            const bool             bra_first,
                                            const bool             diagonal) const
{
    const auto labels = (bra_first) ? t2c::tensor_components(Tensor(integral[1]), "buffer")
                                    : t2c::tensor_components(Tensor(integral[0]), "buffer");
    
    const auto bra_mom = (bra_first) ? SphericalMomentum(integral[0])
                                     : SphericalMomentum(integral[1]);
    
    const auto ket_mom = (bra_first) ? SphericalMomentum(integral[1])
                                     : SphericalMomentum(integral[0]);
    
    const auto index = t2c::tensor_component_index(component);
    
    const auto bra_pairs = bra_mom.select_pairs(index);
    
    for (int i = 0; i < labels.size(); i++)
    {
        for (const auto& ket_pair : ket_mom.select_pairs(i))
        {
            for (const auto& bra_pair : bra_pairs)
            {
                const auto lfactor = t2c::combine_factors(bra_pair.second, ket_pair.second);
                
                auto flabel = labels[i];
                
                flabel += (lfactor == "1.0")  ? "" : ", "  + lfactor;
                
                auto ijlabel = (bra_first) ? std::to_string(bra_pair.first) + ", " + std::to_string(ket_pair.first)
                                           : std::to_string(ket_pair.first) + ", " + std::to_string(bra_pair.first);
                
                if  (integral[0] ==  integral[1])
                {
                    if (diagonal)
                    {
                        lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, " + flabel + ", gto_indexes,"});
                                        
                        lines.push_back({3, 20, 2, ijlabel + ", j, ket_first, ket_last);"});
                    }
                    else
                    {
                        lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, " + flabel + ", bra_gto_indexes, ket_gto_indexes,"});
                        
                        lines.push_back({3, 20, 2, ijlabel + ", j, ket_first, ket_last, mat_type);"});
                    }
                }
                else
                {
                    lines.push_back({3, 0, 1, "t2cfunc::distribute(matrix, " + flabel +
                                               ", bra_gto_indexes, ket_gto_indexes,"});
                        
                    lines.push_back({3, 20, 2, ijlabel + ", j, ket_first, ket_last, ang_order);"});
                    
                }
            }
        }
    }
}

void
T2CFuncBodyDriver::_write_block_distributor(      VCodeLines&      lines,
                                            const TensorComponent& bra_component,
                                            const TensorComponent& ket_component,
                                            const I2CIntegral&     integral,
                                            const bool             diagonal) const
{
    const auto labels = t2c::integrand_components(integral.integrand(), "buffer");
    
    const auto matrices = t2c::integrand_components(integral.integrand(), "matrix");
    
    const auto bra_mom = SphericalMomentum(integral[0]);
    
    const auto bra_index = t2c::tensor_component_index(bra_component);
    
    const auto bra_pairs = bra_mom.select_pairs(bra_index);
    
    const auto ket_mom = SphericalMomentum(integral[1]);
    
    const auto ket_index = t2c::tensor_component_index(ket_component);
    
    const auto ket_pairs = ket_mom.select_pairs(ket_index);
        
    for (int i = 0; i < labels.size(); i++)
    {
        for (const auto& bra_pair : bra_pairs)
        {
            for (const auto& ket_pair : ket_pairs)
            {
                const auto lfactor = t2c::combine_factors(bra_pair.second, ket_pair.second);
                
                auto flabel = labels[i];
                
                flabel += (lfactor == "1.0")  ? "" : ", "  + lfactor;
                
                auto ijlabel = std::to_string(bra_pair.first) + ", " + std::to_string(ket_pair.first);
                
                if  (integral[0] ==  integral[1])
                {
                    if (diagonal)
                    {
                        lines.push_back({3, 0, 1, "t2cfunc::distribute(" + matrices[i] + ", " + flabel +
                                                  ", gto_indexes,"});
                        
                        lines.push_back({3, 20, 2, ijlabel + ", j, ket_first, ket_last);"});
                    }
                    else
                    {
                        lines.push_back({3, 0, 1, "t2cfunc::distribute(" + matrices[i] + ", " + flabel +
                                                  ", bra_gto_indexes, ket_gto_indexes,"});
                        
                        lines.push_back({3, 20, 2, ijlabel + ", j, ket_first, ket_last, mat_type);"});
                    }
                }
                else
                {
                    lines.push_back({3, 0, 1, "t2cfunc::distribute(" + matrices[i] + ", " + flabel +
                                              ", bra_gto_indexes, ket_gto_indexes,"});
                    
                    lines.push_back({3, 20, 2, ijlabel + ", j, ket_first, ket_last, ang_order);"});
                }
            }
        }
    }
}
