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
#include "fraction.hpp"
#include "t2c_utils.hpp"

#include <iostream>

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
    
    if (_need_auxilaries(integral))
    {
        for (const auto& label : _get_auxilaries_def(rgroup))
        {
            lines.push_back({1, 0, 2, label});
        }
    }

    for (const auto& label : _get_batches_def(diagonal))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_batches_loop_start(lines);

    _add_batches_loop_body(lines, diagonal);

    _add_bra_loop_start(lines, integral, diagonal);

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
    
    vstr.push_back("// initialize contracted integral buffers");
    
    size_t ndims = _get_block_size();
    
    if (const auto nterms = rgroup.expansions(); ndims > nterms)
    {
        ndims = nterms;
    }
 
    vstr.push_back("TDoubleArray2D<" + std::to_string(ndims) + "> buffers;");
    
    vstr.push_back("// set up pointers to contracted integral buffers");
   
    for (size_t i = 0; i < ndims; i++)
    {
        vstr.push_back("auto bvals_" + std::to_string(i) + " = buffers[" + std::to_string(i) + "].data();");
    }

    return vstr;
}

std::vector<std::string>
C2CFuncBodyDriver::_get_auxilaries_def(const R2Group& rgroup) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// initialize auxilary buffers");

    const auto ndims = (_get_unique_auxilaries(rgroup)).size();

    vstr.push_back("TDoubleArray2D<" + std::to_string(ndims) + "> auxilaries;");

    vstr.push_back("// set up pointers to auxilary buffers");

    for (size_t i = 0; i < ndims; i++)
    {
        vstr.push_back("auto avals_" + std::to_string(i) + " = auxilaries[" + std::to_string(i) + "].data();");
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
C2CFuncBodyDriver::_add_bra_loop_start(      VCodeLines&  lines,
                                       const I2CIntegral& integral,
                                       const bool         diagonal) const
{
    lines.push_back({2, 0, 1, "for (int64_t j = bra_first; j < bra_last; j++) "});
        
    lines.push_back({2, 0, 1, "{"});
    
    if (_need_auxilaries(integral))
    {
        if (diagonal)
        {
            lines.push_back({3, 0, 2, "const auto bra_coord = gto_coords[j];"});
        }
        else
        {
            lines.push_back({3, 0, 2, "const auto bra_coord = bra_gto_coords[j];"});
        }
        
        lines.push_back({3, 0, 2, "const auto a_x = bra_coord[0];"});
        
        lines.push_back({3, 0, 2, "const auto a_y = bra_coord[1];"});
        
        lines.push_back({3, 0, 2, "const auto a_z = bra_coord[2];"});
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
    
    if (nblocks > 0)
    {
        for (size_t i = 0; i < nblocks; i++)
        {
            const auto first = i * ndims;
            
            const auto last = first + ndims;
            
            _add_bra_loop_block(lines, rgroup, integral, sum_form, diagonal, first, last);
            
            _write_block_distributor(lines, rgroup, integral, sum_form, diagonal, first, last);
        }
    }
    
    if ((rterms % ndims) > 0)
    {
        if (_need_auxilaries(integral))
        {
            _add_bra_loop_block(lines, rgroup, integral, sum_form, diagonal, nblocks * ndims, rterms);
        }
        
        _write_block_distributor(lines, rgroup, integral, sum_form, diagonal, nblocks * ndims, rterms);
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
    std::string blabel = (first != last) ? "(" + std::to_string(first) + "-" + std::to_string(last)  + ")" : "";
    
    lines.push_back({3, 0, 2, "// compute integrals batch " + blabel});
    
    lines.push_back({3, 0, 1, "#pragma omp simd aligned(ket_coords_x, ket_coords_y, ket_coords_z : 64)"});
    
    lines.push_back({3, 0, 1, "for (int64_t k = 0; k < ket_dim; k++)"});
    
    lines.push_back({3, 0, 1, "{"});
    
    _add_loop_prefactors(lines, rgroup, sum_form, first, last); 
    
    const auto auxilaries = _get_unique_auxilaries(rgroup);
    
    for (auto i = first; i < last; i++)
    {
        _add_bra_loop_line(lines, rgroup[i], integral, auxilaries,  i - first, sum_form);
    }
    
    lines.push_back({3, 0, 2, "}"});
}

void
C2CFuncBodyDriver::_write_block_distributor(      VCodeLines&  lines,
                                            const R2Group&     rgroup,
                                            const I2CIntegral& integral,
                                            const bool         sum_form,
                                            const bool         diagonal,
                                            const size_t       first,
                                            const size_t       last) const
{
    const auto bra_mom = SphericalMomentum(integral[0]);
    
    const auto ket_mom = SphericalMomentum(integral[1]);
    
    for (auto i = first; i < last; i++)
    {
        const auto tint = rgroup[i].root().integral();
        
        const auto bra_index = t2c::tensor_component_index(tint[0]);
        
        const auto ket_index = t2c::tensor_component_index(tint[1]);
        
        const auto mlabel = _get_matrix_label(tint);
        
        for (const auto& bra_pair : bra_mom.select_pairs(bra_index))
        {
            for (const auto& ket_pair : ket_mom.select_pairs(ket_index))
            {
                const auto lfactor = t2c::combine_factors(bra_pair.second, ket_pair.second);
                
                auto flabel = "buffer[" + std::to_string(i - first) + "]";
                
                flabel += (lfactor == "1.0")  ? "" : ", "  + lfactor;
                
                auto ijlabel = std::to_string(bra_pair.first) + ", " + std::to_string(ket_pair.first);
                
                if  (integral[0] ==  integral[1])
                {
                    if (diagonal)
                    {
                        lines.push_back({3, 0, 2, "t2cfunc::distribute(" + mlabel + ", " + flabel +
                                                  ", gto_indexes," + ijlabel + ", j, ket_first, ket_last);"});
                    }
                    else
                    {
                        lines.push_back({3, 0, 2, "t2cfunc::distribute(" + mlabel + ", " + flabel +
                                                  ", bra_gto_indexes, ket_gto_indexes," + ijlabel + ", j, ket_first, ket_last, mat_type);"});
                    }
                }
                else
                {
                    lines.push_back({3, 0, 2, "t2cfunc::distribute(" + mlabel + ", " + flabel +
                                              ", bra_gto_indexes, ket_gto_indexes," + ijlabel + ", j, ket_first, ket_last, ang_order);"});
                }
            }
        }
    }
}

void
C2CFuncBodyDriver::_add_bra_loop_line(      VCodeLines&   lines,
                                      const R2CDist&      rdist,
                                      const I2CIntegral&  integral,
                                      const V4Auxilaries& auxilaries,
                                      const size_t        index,
                                      const bool          sum_form) const
{
    std::string code = "bvals_" + std::to_string(index) + "[k] = ";
    
    for (size_t i = 0; i < rdist.terms(); i++)
    {
        code += _get_rterm_code(rdist[i], auxilaries, i == 0);
    }
        
    lines.push_back({4, 0, 2, code + ";"});
}

std::string
C2CFuncBodyDriver::_get_rterm_code(const R2CTerm&      rterm,
                                   const V4Auxilaries& auxilaries,
                                   const bool          is_first) const
{
    const auto pre_fact = rterm.prefactor();
        
    auto plabel = pre_fact.label();
        
    if (plabel == "1.0")  plabel = "";
        
    if (plabel == "-1.0") plabel = "-";
        
    if (pre_fact.denominator() != 1)
    {
        if (pre_fact.numerator() < 0) plabel.erase(0, 1);
            
        plabel = "(" + plabel + ")";
            
        if (pre_fact.numerator() < 0) plabel = "-" + plabel;
    }
    
    if (!is_first)
    {
        if (plabel[0] == '-')
        {
            plabel.insert(1, " ");
        }
        else
        {
            plabel = "+ " + plabel;
        }
        
        plabel = " " + plabel;
    }
        
    const auto facts = rterm.factors();
        
    std::string flabel;
        
    for (const auto& fact : facts)
    {
        if (fact == Factor("N", "n")) continue;
        
        if (fact == Factor("M", "m")) continue;
        
        if (fact == Factor("T", "t")) continue;
        
        const auto norder = rterm.factor_order(fact);
            
        for (size_t n = 0; n < norder; n++)
        {
            flabel += " * " + fact.label();
        }
    }
        
    // remove multiplication for special cases
        
    if ((pre_fact == Fraction(1)) || (pre_fact == Fraction(-1)))
    {
        flabel.erase(0, 3);
    }
        
    // merge labels

    flabel = plabel + flabel;
    
    // add auxilary label
    
    const auto index = _get_auxilary_index(auxilaries, _get_auxilary(rterm));
    
    flabel += " * f_" + std::to_string(index);
        
    return flabel;
}

void
C2CFuncBodyDriver::_add_loop_prefactors(      VCodeLines& lines,
                                        const R2Group&    rgroup,
                                        const bool        sum_form, 
                                        const size_t      first,
                                        const size_t      last) const
{
    int spacer = (sum_form) ? 5 : 4;
    
    if (t2c::find_factor(rgroup, "rab_x", first, last))
    {
        lines.push_back({spacer, 0, 2, "const auto rab_x = a_x - ket_coords_x[k];"});
    }
    
    if (t2c::find_factor(rgroup, "rab_y", first, last))
    {
        lines.push_back({spacer, 0, 2, "const auto rab_y = a_y - ket_coords_y[k];"});
    }
    
    if (t2c::find_factor(rgroup, "rab_z", first, last))
    {
        lines.push_back({spacer, 0, 2, "const auto rab_z = a_z - ket_coords_z[k];"});
    }
    
    const auto auxilaries = _get_unique_auxilaries(rgroup);
    
    for (const auto& taux : _get_unique_auxilaries(rgroup, first, last))
    {
        const auto ilabel = std::to_string(_get_auxilary_index(auxilaries, taux));
        
        lines.push_back({spacer, 0, 2, "const auto f_" + ilabel + " = avals_" + ilabel + "[k];"});
    }
}

std::string
C2CFuncBodyDriver::_get_matrix_label(const T2CIntegral& integral) const
{
    const auto prefixes = integral.prefixes();
    
    std::string label("matrix");
    
    std::string olabel;
    
    if (const auto integrand = integral.integrand(); integrand.shape().order() > 0)
    {
        olabel = "_" + integrand.shape().label();
    }
    
    if (prefixes.size() == 0)
    {
        return label + olabel;
    }
    
    if (prefixes.size() == 1)
    {
        std::string blabel = "_" + prefixes[0].shape().label();
        
        return label + blabel + olabel;
    }
    
    if (prefixes.size() == 2)
    {
        std::string blabel = "_" + prefixes[0].shape().label();
        
        std::string klabel = "_" + prefixes[1].shape().label();
        
        return label + blabel + klabel + olabel;
    }
    
    return label;
}

size_t
C2CFuncBodyDriver::_get_block_size() const
{
    return 15;
}

V4Auxilaries
C2CFuncBodyDriver::_get_unique_auxilaries(const R2Group& rgroup) const
{
    V4Auxilaries auxs;
    
    for (size_t i = 0; i < rgroup.expansions(); i++)
    {
        for (size_t j = 0; j < rgroup[i].terms(); j++)
        {
            auxs.insert(_get_auxilary(rgroup[i][j]));
        }
    }
   
    return auxs;
}

V4Auxilaries
C2CFuncBodyDriver::_get_unique_auxilaries(const R2Group& rgroup,
                                          const size_t   first,
                                          const size_t   last) const
{
    V4Auxilaries auxs;
    
    for (size_t i = first; i < last; i++)
    {
        for (size_t j = 0; j < rgroup[i].terms(); j++)
        {
            auxs.insert(_get_auxilary(rgroup[i][j]));
        }
    }
   
    return auxs;
}

size_t
C2CFuncBodyDriver::_get_auxilary_index(const V4Auxilaries& auxilaries,
                                       const T4Index&      target) const
{
    size_t index = 0;
    
    for (const auto& taux : auxilaries)
    {
        if (taux == target) return index;
        
        index++;
    }
   
    return -1;
}

T4Index
C2CFuncBodyDriver::_get_auxilary(const R2CTerm& rterm) const
{
    const auto n = rterm.factor_order(Factor("N", "n"));
    
    const auto m = rterm.factor_order(Factor("M", "m"));
    
    const auto t = rterm.factor_order(Factor("T", "t"));
    
    const auto p = rterm.order();
    
    return {n, m, t, p};
}

bool
C2CFuncBodyDriver::_need_auxilaries(const I2CIntegral& integral) const
{
    if ((integral[0] + integral[1]) == 0)
    {
        const auto integrand = integral.integrand();
        
        if (integrand.name() == "1") return false;
    }
    
    return true;
}
