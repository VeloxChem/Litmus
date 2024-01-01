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

#include "c2c_body.hpp"

#include "string_formater.hpp"
#include "spherical_momentum.hpp"
#include "t2c_utils.hpp"

void
C2CFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                   const R2Group&       rgroup,
                                   const I2CIntegral&   integral,
                                   const bool           sum_form,
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

    for (const auto& label : _get_buffers_def(rgroup))
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

    _add_bra_loop_body(lines, rgroup, integral, sum_form, diagonal);

    _add_bra_loop_end(lines);

    _add_batches_loop_end(lines);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
C2CFuncBodyDriver::_get_angmom_def(const I2CIntegral& integral) const
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
C2CFuncBodyDriver::_get_gtos_def(const bool diagonal) const
{
    std::vector<std::string> vstr;
    
    if (diagonal)
    {
        vstr.push_back("// intialize GTOs data");
        
        vstr.push_back("const auto gto_coords = gto_block.getCoordinates();");
        
        vstr.push_back("const auto gto_indexes = gto_block.getOrbitalIndexes();");
       
        vstr.push_back("const auto ncgtos = gto_block.getNumberOfBasisFunctions();");
    }
    else
    {
        vstr.push_back("// intialize GTOs data on bra side");
        
        vstr.push_back("const auto bra_gto_coords = bra_gto_block.getCoordinates();");
        
        vstr.push_back("const auto bra_gto_indexes = bra_gto_block.getOrbitalIndexes();");
       
        vstr.push_back("const auto bra_ncgtos = bra_gto_block.getNumberOfBasisFunctions();");
        
        vstr.push_back("// intialize GTOs data on ket side");
        
        vstr.push_back("const auto ket_gto_coords = ket_gto_block.getCoordinates();");
        
        vstr.push_back("const auto ket_gto_indexes = ket_gto_block.getOrbitalIndexes();");
       
        vstr.push_back("const auto ket_ncgtos = ket_gto_block.getNumberOfBasisFunctions();");
    }
    
    return vstr;
}

std::vector<std::string>
C2CFuncBodyDriver::_get_ket_variables_def() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// initialize aligned arrays for ket side");
        
    vstr.push_back("alignas(64) TDoubleArray ket_coords_x;");
        
    vstr.push_back("alignas(64) TDoubleArray ket_coords_y;");
        
    vstr.push_back("alignas(64) TDoubleArray ket_coords_z;");
    
    return vstr;
}

std::vector<std::string>
C2CFuncBodyDriver::_get_buffers_def(const R2Group& rgroup) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// initialize contracted integrals buffer");
    
    size_t ndims = _get_block_size();
    
    if (const auto nterms = rgroup.expansions(); ndims > nterms)
    {
        ndims = nterms;
    }
 
    vstr.push_back("TDoubleArray2D<" + std::to_string(ndims) + "> buffers;");
    
    vstr.push_back("// set up pointers to contracted integrals buffer");
   
    for (size_t i = 0; i < ndims; i++)
    {
        vstr.push_back("auto buffer_" + std::to_string(i) + " = buffers[" + std::to_string(i) + "].data();");
    }

    return vstr;
}

std::vector<std::string>
C2CFuncBodyDriver::_get_batches_def(const bool diagonal) const
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
C2CFuncBodyDriver::_add_batches_loop_start(VCodeLines& lines) const
{
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < nbatches; i++)"});
        
    lines.push_back({1, 0, 1, "{"});
}

void
C2CFuncBodyDriver::_add_batches_loop_body(      VCodeLines& lines,
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
            
    lines.push_back({2, 0, 2, label});
}

void
C2CFuncBodyDriver::_add_batches_loop_end(VCodeLines& lines) const
{
    lines.push_back({1, 0, 1, "}"});
}

void
C2CFuncBodyDriver::_add_bra_loop_start(      VCodeLines& lines,
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
C2CFuncBodyDriver::_add_bra_loop_end(VCodeLines& lines) const
{
    lines.push_back({2, 0, 1, "}"});
}

void
C2CFuncBodyDriver::_add_bra_loop_body(      VCodeLines&  lines,
                                      const R2Group&     rgroup, 
                                      const I2CIntegral& integral,
                                      const bool         sum_form,
                                      const bool         diagonal) const
{
    const auto ndims = _get_block_size();
    
    const auto rterms = rgroup.expansions();
    
    const auto nblocks = rterms / ndims;
    
    const auto rblocks = rterms % ndims;
    
    if (nblocks > 0)
    {
        for (size_t i = 0; i < nblocks; i++)
        {
            const auto first = i * ndims;
            
            const auto last = first + ndims;
            
            _add_bra_loop_block(lines, rgroup, integral, sum_form, diagonal, first, last);
        }
    }
    
    if (rterms > 0)
    {
        _add_bra_loop_block(lines, rgroup, integral, sum_form, diagonal, nblocks * ndims, rterms);
    }
}

void
C2CFuncBodyDriver::_add_bra_loop_block(      VCodeLines&  lines,
                                       const R2Group&     rgroup,
                                       const I2CIntegral& integral,
                                       const bool         sum_form,
                                       const bool         diagonal,
                                       const size_t       first,
                                       const size_t       last) const
{
    
}

size_t
C2CFuncBodyDriver::_get_block_size() const
{
    return 15;
}
