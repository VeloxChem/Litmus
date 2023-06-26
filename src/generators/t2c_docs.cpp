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
                             const bool           diagonal) const
{
    auto lines = VCodeLines();
        
    lines.push_back({0, 0, 1, "/**"});
        
    lines.push_back({0, 0, 2, _get_compute_str(integral, diagonal)});
    
    for (const auto& label : _get_matrix_str(integral))
    {
        lines.push_back({0, 1, 1, label});
    }
    
    for (const auto& label : _get_special_vars_str(integral, true))
    {
        lines.push_back({0, 1, 1, label});
    }
    
    for (const auto& label : _get_gto_blocks_str(integral, diagonal))
    {
        lines.push_back({0, 1, 1, label});
    }
    
    for (const auto& label : _get_indexes_str())
    {
        lines.push_back({0, 1, 1, label});
    }
    
    if (const auto label = _get_matrix_type_str(integral, diagonal);
        !label.empty())
    {
        lines.push_back({0, 1, 1, label});
    }
        
    lines.push_back({0, 0, 1, "*/"});
        
    ost::write_code_lines(fstream, lines);
}

void
T2CDocuDriver::write_prim_doc_str(      std::ofstream& fstream,
                                  const I2CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "/**"});
    
    lines.push_back({0, 1, 2, _get_prim_compute_str(integral)});
    
    for (const auto& label : _get_prim_buffer_str(integral))
    {
        lines.push_back({0, 1, 1, label});
    }
    
    for (const auto& label : _get_special_vars_str(integral, true))
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
                                  const bool             bra_first) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "/**"});
    
    lines.push_back({0, 1, 2, _get_prim_compute_str(component, integral, bra_first)});
    
    for (const auto& label : _get_prim_buffer_str(integral, bra_first))
    {
        lines.push_back({0, 1, 1, label});
    }
    
    for (const auto& label : _get_special_vars_str(integral, true))
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
                                  const I2CIntegral&     integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "/**"});
    
    lines.push_back({0, 1, 2, _get_prim_compute_str(bra_component, ket_component, integral)});
    
    for (const auto& label : _get_prim_buffer_str(integral))
    {
        lines.push_back({0, 1, 1, label});
    }
    
    for (const auto& label : _get_special_vars_str(integral, true))
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
        
    auto label = " Evaluates <" + bra.label() + "|";
        
    label += t2c::integrand_label(integral.integrand());
        
    label += "|" + ket.label() + ">  integrals for given ";
        
    label += (diagonal) ? "GTOs block." : "pair of GTOs blocks.";
    
    return label;
}

std::string
T2CDocuDriver::_get_prim_compute_str(const I2CIntegral& integral) const
{
    const auto bra = Tensor(integral[0]);
    
    const auto ket = Tensor(integral[1]);
    
    auto label = "Evaluates block of primitive <" + bra.label() + "|" ;
    
    label += t2c::integrand_label(integral.integrand());
    
    label += "|" + ket.label() + "> integrals.";
    
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
    
    auto label = "Evaluates block of primitive <" + bra.label();
    
    label += "_" + fstr::upcase(bra_component.label());
    
    label += "|" + t2c::integrand_label(integral.integrand()) + "|";
    
    label += ket.label() + "_" + fstr::upcase(ket_component.label());
    
    label += "> integrals.";
    
    return label;
}

std::vector<std::string>
T2CDocuDriver::_get_matrix_str(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (const auto labels = t2c::integrand_components(integral.integrand(), "matrix");
        labels.size() == 1)
    {
        vstr.push_back("@param matrix the pointer to matrix for storage of integrals.");
    }
    else
    {
        for (const auto& label : labels)
        {
            auto lcomp = fstr::upcase(label);
            
            lcomp.erase(0, lcomp.find('_') + 1);
            
            vstr.push_back("@param " + label + "the pointer to matrix for storage " +
                           "of Cartesian integral component " + lcomp + ".");
        }
    }
    
    return vstr;
}

std::vector<std::string>
T2CDocuDriver::_get_special_vars_str(const I2CIntegral& integral,
                                     const bool         geom_form) const
{
    std::vector<std::string> vstr;
    
    // nuclear potential integrals
    
    if (integral.integrand() == Operator("A"))
    {
        if (geom_form)
        {
            vstr.push_back("@param charge the charge of external point.");
            
            vstr.push_back("@param point the coordinates of external point.");
        }
        else
        {
            vstr.push_back("@param charges the vector of charges.");
            
            vstr.push_back("@param points the vector of coordinates of external points.");
        }
    }
    
    // nuclear potential first derivative integrals
    
    if (integral.integrand() == Operator("AG", Tensor(1)))
    {
        if (geom_form)
        {
            vstr.push_back("@param dipole the charge of external point.");
            
            vstr.push_back("@param point the coordinates of external point.");
        }
        else
        {
            vstr.push_back("@param dipoles the vector of charges.");
            
            vstr.push_back("@param points the vector of coordinates of external points.");
        }
    }
    
    return vstr;
}

std::vector<std::string>
T2CDocuDriver::_get_gto_blocks_str(const I2CIntegral& integral,
                                   const bool         diagonal) const
{
    std::vector<std::string> vstr;
    
    if (diagonal)
    {
       vstr.push_back("@param gto_block the GTOs block.");
    }
    else
    {
        vstr.push_back("@param bra_gto_block the GTOs block on bra side.");
        
        vstr.push_back("@param ket_gto_block the GTOs block on ket side.");
    }
    
    if (integral[0] != integral[1])
    {
        vstr.push_back("@param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.");
    }
    
    return vstr;
}

std::vector<std::string>
T2CDocuDriver::_get_indexes_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("@param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.");
    
    vstr.push_back("@param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.");
    
    return vstr;
}

std::string
T2CDocuDriver::_get_matrix_type_str(const I2CIntegral& integral,
                                    const bool         diagonal) const
{
    if ((!diagonal) && (integral[0] == integral[1]))
    {
        return std::string("@param mat_type the matrix type.");
    }
    else
    {
        return std::string();
    }
}

std::vector<std::string>
T2CDocuDriver::_get_prim_buffer_str(const I2CIntegral& integral) const
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
        
        for (const auto& label : t2c::integrand_components(integral.integrand(), "buffer"))
        {
            vstr.push_back("@param " + label + " the partial integrals buffer.");
        }
        
        return vstr;
    }
}

std::vector<std::string>
T2CDocuDriver::_get_prim_buffer_str(const I2CIntegral& integral,
                                    const bool         bra_first) const
{
    std::vector<std::string> vstr;
    
    const auto tensor = (bra_first) ? Tensor(integral[1]) : Tensor(integral[0]);

    for (const auto& label : t2c::tensor_components(tensor, "buffer"))
    {
        vstr.push_back("@param " + label + " the partial integrals buffer.");
    }
 
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
