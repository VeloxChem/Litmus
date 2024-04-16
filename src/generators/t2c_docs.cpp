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

#include "t2c_docs.hpp"

#include "file_stream.hpp"
#include "t2c_utils.hpp"
#include "string_formater.hpp"

void
T2CDocuDriver::write_doc_str(      std::ofstream& fstream,
                             const I2CIntegral&   integral,
                             const bool           sum_form,
                             const bool           diagonal) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_compute_str(integral, diagonal)});
    
    for (const auto& label : _get_matrix_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_special_vars_str(integral, sum_form))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_gto_blocks_str(integral, false, diagonal))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_indexes_str())
    {
        lines.push_back({0, 0, 1, label});
    }
        
    ost::write_code_lines(fstream, lines);
}

void
T2CDocuDriver::write_prim_doc_str(      std::ofstream& fstream,
                                  const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, _get_prim_compute_str(integral)});
    
    for (const auto& label : _get_prim_buffer_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    for (const auto& label : _get_prim_variables_str(integral))
    {
        lines.push_back({0, 0, 1, label});
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T2CDocuDriver::write_auxilary_doc_str(      std::ofstream& fstream,
                                      const I2CIntegral&   integral,
                                      const bool           diagonal) const
{
    auto lines = VCodeLines();
        
    lines.push_back({0, 0, 1, "/**"});
        
    lines.push_back({0, 0, 2, _get_auxilary_compute_str(integral, diagonal)});
    
    lines.push_back({0, 1, 1, "@param auxilaries the buffer for auxilary integrals."});
    
    for (const auto& label : _get_gto_blocks_str(integral, true, diagonal))
    {
        lines.push_back({0, 1, 1, label});
    }
    
    for (const auto& label : _get_auxilary_indexes_str())
    {
        lines.push_back({0, 1, 1, label});
    }
            
    lines.push_back({0, 0, 1, "*/"});
        
    ost::write_code_lines(fstream, lines);
}
void
T2CDocuDriver::write_prim_doc_str(      std::ofstream& fstream,
                                  const I2CIntegral&   integral,
                                  const bool           sum_form) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "/**"});
    
    lines.push_back({0, 1, 2, _get_prim_compute_str(integral)});
    
    for (const auto& label : _get_prim_buffer_str(integral))
    {
        lines.push_back({0, 1, 1, label});
    }
    
    for (const auto& label : _get_special_vars_str(integral, sum_form))
    {
        lines.push_back({0, 1, 1, label});
    }
    
    for (const auto& label : _get_prim_variables_str())
    {
        lines.push_back({0, 1, 1, label});
    }
    
    lines.push_back({0, 0, 1, "*/"});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CDocuDriver::write_prim_doc_str(      std::ofstream&   fstream,
                                  const TensorComponent& component,
                                  const I2CIntegral&     integral,
                                  const bool             sum_form,
                                  const bool             bra_first) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "/**"});
    
    lines.push_back({0, 1, 2, _get_prim_compute_str(component, integral, bra_first)});
    
    for (const auto& label : _get_prim_buffer_str(integral, bra_first))
    {
        lines.push_back({0, 1, 1, label});
    }
    
    for (const auto& label : _get_special_vars_str(integral, sum_form))
    {
        lines.push_back({0, 1, 1, label});
    }
    
    for (const auto& label : _get_prim_variables_str())
    {
        lines.push_back({0, 1, 1, label});
    }
    
    lines.push_back({0, 0, 1, "*/"});
    
    ost::write_code_lines(fstream, lines);
}

void
T2CDocuDriver::write_prim_doc_str(      std::ofstream&   fstream,
                                  const TensorComponent& bra_component,
                                  const TensorComponent& ket_component,
                                  const I2CIntegral&     integral,
                                  const bool             sum_form) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "/**"});
    
    lines.push_back({0, 1, 2, _get_prim_compute_str(bra_component, ket_component, integral)});
    
    for (const auto& label : _get_prim_buffer_str(integral))
    {
        lines.push_back({0, 1, 1, label});
    }
    
    for (const auto& label : _get_special_vars_str(integral, sum_form))
    {
        lines.push_back({0, 1, 1, label});
    }
    
    for (const auto& label : _get_prim_variables_str())
    {
        lines.push_back({0, 1, 1, label});
    }
    
    lines.push_back({0, 0, 1, "*/"});
    
    ost::write_code_lines(fstream, lines);
}

std::string
T2CDocuDriver::_get_compute_str(const I2CIntegral& integral,
                                const bool         diagonal) const
{
    const auto bra = Tensor(integral[0]);
        
    const auto ket = Tensor(integral[1]);
        
    const auto integrand = integral.integrand();
    
    const auto prefixes = integral.prefixes();
    
    auto bra_geom = std::string("");
    
    auto ket_geom = std::string("");
    
    if (const auto nterms = prefixes.size(); nterms > 0)
    {
        if (nterms >= 1)
        {
            const auto border = std::to_string(prefixes[0].shape().order());
            
            bra_geom = "d^(" + border + ")/dA^(" + border + ")";
        }
        
        if (nterms >= 2)
        {
            const auto korder = std::to_string(prefixes[1].shape().order());
            
            ket_geom = "d^(" + korder + ")/dB^(" + korder + ")";
        }
    }
        
    auto label = "/// Evaluates <" + bra_geom + bra.label() + "|";
    
    if (integral.integrand().name() != "1")
    {
        label += t2c::integrand_label(integral.integrand()) + "|";
    }
    
    label += ket_geom + ket.label() + ">  integrals for given ";
        
    label += (diagonal) ? "GTOs block." : "pair of GTOs blocks.";
    
    return label;
}

std::string
T2CDocuDriver::_get_auxilary_compute_str(const I2CIntegral& integral,
                                         const bool         diagonal) const
{
    std::string label = " Evaluates (m|";
        
    label += t2c::integrand_label(integral.integrand());
        
    label += "|n)_t,p  auxilary integrals for given ";
        
    label += (diagonal) ? "GTOs block." : "pair of GTOs blocks.";
    
    return label;
}

std::string
T2CDocuDriver::_get_prim_compute_str(const I2CIntegral& integral) const
{
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    auto label = "/// Evaluates block of primitive <"  + bra.label() + "|" ;
    
    if (integral.integrand().name() != "1")
    {
        label += t2c::integrand_label(integral.integrand()) + "|";
    }
    
    label += ket.label() + "> integrals.";
    
    return label;
}

std::string
T2CDocuDriver::_get_prim_compute_str(const TensorComponent& component,
                                     const I2CIntegral&     integral,
                                     const bool             bra_first) const
{
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    auto label = "Evaluates block of primitive <" + bra.label();
    
    label += (bra_first) ? "_" + fstr::upcase(component.label()) : "";
    
    label += "|" + t2c::integrand_label(integral.integrand()) + "|";
    
    label += ket.label();
    
    label += (bra_first) ? "" : "_" + fstr::upcase(component.label());
    
    label += ">  integrals.";
    
    return label;
}

std::string
T2CDocuDriver::_get_prim_compute_str(const TensorComponent& bra_component,
                                     const TensorComponent& ket_component,
                                     const I2CIntegral&     integral) const
{
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    const auto prefixes = integral.prefixes();
    
    auto bra_geom = std::string("");
    
    auto ket_geom = std::string("");
    
    if (const auto nterms = prefixes.size(); nterms > 0)
    {
        if (nterms >= 1)
        {
            const auto border = std::to_string(prefixes[0].shape().order());
            
            bra_geom = "d^(" + border + ")/dA^(" + border + ")";
        }
        
        if (nterms >= 2)
        {
            const auto korder = std::to_string(prefixes[1].shape().order());
            
            ket_geom = "d^(" + korder + ")/dB^(" + korder + ")";
        }
    }
    
    auto label = "Evaluates block of primitive <" + bra_geom + bra.label();
    
    label += "_" + fstr::upcase(bra_component.label());
    
    label += "|" + t2c::integrand_label(integral.integrand()) + "|";
    
    label += ket_geom +  ket.label() + "_" + fstr::upcase(ket_component.label());
    
    label += "> integrals.";
    
    return label;
}

std::vector<std::string>
T2CDocuDriver::_get_matrix_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// - Parameter matrix: the pointer to matrix for storage of integrals.");
    
    return vstr;
}

std::vector<std::string>
T2CDocuDriver::_get_special_vars_str(const I2CIntegral& integral,
                                     const bool         sum_form) const
{
    std::vector<std::string> vstr;
    
    // nuclear potential integrals
    
    if (integral.integrand() == Operator("A"))
    {
        if (sum_form)
        {
            vstr.push_back("@param charges the vector of charges.");
            
            vstr.push_back("@param points the vector of coordinates of external points.");
        }
        else
        {
            vstr.push_back("@param charge the charge of external point.");
            
            vstr.push_back("@param point the coordinates of external point.");
        }
    }
    
    // nuclear potential first derivative integrals
    
    if (integral.integrand() == Operator("AG", Tensor(1)))
    {
        if (sum_form)
        {
            vstr.push_back("@param dipoles the vector of dipoles.");
            
            vstr.push_back("@param points the vector of coordinates of external points.");
        }
        else
        {
            vstr.push_back("@param dipole the dipole of external point.");
            
            vstr.push_back("@param point the coordinates of external point.");
        }
    }
    
    // nuclear potential second derivative integrals
    
    if (integral.integrand() == Operator("AG", Tensor(2)))
    {
        if (sum_form)
        {
            vstr.push_back("@param quadrupoles the vector of quadrupoles.");
            
            vstr.push_back("@param points the vector of coordinates of external points.");
        }
        else
        {
            vstr.push_back("@param quadrupole the quadrupole of external point.");
            
            vstr.push_back("@param point the coordinates of external point.");
        }
    }
    
    // multipole integrals
    
    if (integral.integrand().name() == "r")
    {
        vstr.push_back("@param point the coordinates of external point.");
    }
    
    // three center overlap integrals
    
    if (integral.integrand().name() == "G(r)")
    {
        vstr.push_back("@param gau_exp the exponent of external Gaussian.");
        
        vstr.push_back("@param gau_center the coordinates of external Gaussian center.");
    }
    
    return vstr;
}

std::vector<std::string>
T2CDocuDriver::_get_gto_blocks_str(const I2CIntegral& integral,
                                   const bool         is_auxilary,
                                   const bool         diagonal) const
{
    std::vector<std::string> vstr;
    
    if (diagonal)
    {
       vstr.push_back("/// - Parameter gto_block: the GTOs block.");
    }
    else
    {
        vstr.push_back("/// - Parameter bra_gto_block: the GTOs block on bra side.");
        
        vstr.push_back("/// - Parameter ket_gto_block: the GTOs block on ket side.");
        
        if (integral[0] != integral[1])
        {
            vstr.push_back("/// - Parameter ang_order: the flag for matching angular order between matrix and pair of GTOs blocks.");
        }
        else
        {
            vstr.push_back("/// - Parameter mat_type: the matrix type.");
        }
    }
        
    return vstr;
}

std::vector<std::string>
T2CDocuDriver::_get_indexes_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("/// - Parameter bra_indices: the range [bra_first, bra_last) of GTOs on bra side.");
    
    vstr.push_back("/// - Parameter ket_indices: the range [ket_first, ket_last) of GTOs on ket side.");
    
    return vstr;
}

std::vector<std::string>
T2CDocuDriver::_get_auxilary_indexes_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("@param bra_index the index of GTO on bra side.");
    
    vstr.push_back("@param ket_igtos the range [ket_first, ket_last) of GTOs on ket side.");
    
    return vstr;
}

std::string
T2CDocuDriver::_get_matrix_type_str(const I2CIntegral& integral,
                                    const bool         diagonal) const
{
    return std::string();
}

std::vector<std::string>
T2CDocuDriver::_get_prim_buffer_str(const I2CIntegral& integral,
                                    const bool         bra_first) const
{
    if (integral.is_simple_integrand() && integral.is_simple())
    {
        std::vector<std::string> vstr;
        
        if (integral[0] > 0)
        {
            for(const auto& label : t2c::tensor_components(Tensor(integral[0]), "buffer"))
            {
                vstr.push_back("@param " + label + " the partial integrals buffer.");
            }
        }
        
        if (integral[1] > 0)
        {
            for(const auto& label : t2c::tensor_components(Tensor(integral[1]), "buffer"))
            {
                vstr.push_back("@param " + label + " the partial integrals buffer.");
            }
        }
        
        if (vstr.empty())
        {
            vstr.push_back("@param buffer the integrals buffer.");
        }
        
        return vstr;
    }
    else
    {
        std::vector<std::string> vstr;
        
        const auto prefixes = integral.prefixes();
        
        if (prefixes.empty())
        {
            for (const auto& label : t2c::integrand_components(integral.integrand(), "buffer"))
            {
                vstr.push_back("@param " + label + " the partial integrals buffer.");
            }
        }
        
        if (prefixes.size() == 1)
        {
            for (const auto& label : t2c::integrand_components(prefixes[0].shape(), integral.integrand(), "buffer"))
            {
                vstr.push_back("@param " + label + " the partial integrals buffer.");
            }
        }
        
        if (prefixes.size() == 2)
        {
            for (const auto& label : t2c::integrand_components(prefixes[0].shape(), prefixes[1].shape(), integral.integrand(), "buffer"))
            {
                vstr.push_back("@param " + label + " the partial integrals buffer.");
            }
        }
        
        return vstr;
    }
}

std::vector<std::string>
T2CDocuDriver::_get_prim_buffer_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    for (const auto& tint : t2c::get_integrals(integral))
    {
        vstr.push_back("/// - Parameter " + t2c::get_buffer_label(tint, "prim") + ": the primitive integrals buffer.");
    }
    
    vstr.push_back("/// - Parameter " + t2c::get_buffer_label(integral, "prim") + ": the primitive integrals buffer.");
    
    return vstr;
}

std::vector<std::string>
T2CDocuDriver::_get_prim_variables_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("@param bra_exp the primitive exponent on bra side.");
    
    vstr.push_back("@param bra_norm the primitive normalization factor on bra side.");
    
    vstr.push_back("@param bra_coord the 3d coordinate of basis function on bra side.");
    
    vstr.push_back("@param ket_exps the array of primitive exponents on ket side.");
    
    vstr.push_back("@param ket_norms the array of primitive normalization factors on ket side.");
    
    vstr.push_back("@param ket_coords_x the array of Cartesian X coordinates on ket side.");
    
    vstr.push_back("@param ket_coords_y the array of Cartesian Y coordinates on ket side.");
    
    vstr.push_back("@param ket_coords_z the array of Cartesian Z coordinates on ket side.");
    
    vstr.push_back("@param ket_dim the end size of ket arrays.");
    
    return vstr;
}

std::vector<std::string>
T2CDocuDriver::_get_prim_variables_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (integral[0] > 0)
    {
        vstr.push_back("/// - Parameter pa_x: the vector of Cartesian X  distances R(PA) = P - A.");
        
        vstr.push_back("/// - Parameter pa_y: the vector of Cartesian Y  distances R(PA) = P - A.");
        
        vstr.push_back("/// - Parameter pa_z: the vector of Cartesian Z  distances R(PA) = P - A.");
    }
    
    if ((integral[0] == 0) && (integral[1] > 0))
    {
        vstr.push_back("/// - Parameter pb_x: the vector of Cartesian X  distances R(PB) = P - B.");
        
        vstr.push_back("/// - Parameter pb_y: the vector of Cartesian Y  distances R(PB) = P - B.");
        
        vstr.push_back("/// - Parameter pb_z: the vector of Cartesian Z  distances R(PB) = P - B.");
    }
    
    if ((integral[0] + integral[1]) == 0)
    {
        if (integral.integrand().name() == "1")
        {
            vstr.push_back("/// - Parameter ab_x: the vector of Cartesian X  distances R(AB) = A - B.");
            
            vstr.push_back("/// - Parameter ab_y: the vector of Cartesian Y  distances R(AB) = A - B.");
            
            vstr.push_back("/// - Parameter ab_z: the vector of Cartesian Z  distances R(AB) = A - B.");
        }
    }
    
    if ((integral[0] + integral[1]) != 1) 
    {
        vstr.push_back("/// - Parameter a_exp: the GTOs exponent on center A.");
        
        vstr.push_back("/// - Parameter b_exps: the vector of GTOs exponents on center B.");
    }
    
    if ((integral[0] + integral[1]) == 0)
    {
        if (integral.integrand().name() == "1")
        {
            vstr.push_back("/// - Parameter a_norm: the GTOs normalization factor on center A.");
            
            vstr.push_back("/// - Parameter b_norms: the vector of GTOs normalization factors on center B.");
        }
    }
    
    return vstr;
}

