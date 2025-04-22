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

#include "g2c_body.hpp"

#include "t2c_utils.hpp"

void
G2CFuncBodyDriver::write_func_body(      std::ofstream&         fstream,
                                   const SI2CIntegrals&         geom_integrals,
                                   const SI2CIntegrals&         vrr_integrals,
                                   const I2CIntegral&           integral,
                                   const std::array<int, 3>& geom_drvs,
                                   const bool                   use_rs) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_gtos_def())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_variables_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }

    _add_loop_start(lines, integral);
    
    _add_call_tree(lines, vrr_integrals, integral);
    
    _add_loop_end(lines, vrr_integrals, integral);
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
G2CFuncBodyDriver::_get_gtos_def() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// intialize GTOs data on bra side");
        
    vstr.push_back("const auto bra_gto_exps = bra_gto_block.exponents();");
        
    vstr.push_back("const auto bra_gto_norms = bra_gto_block.normalization_factors();");
        
    vstr.push_back("const auto bra_ncgtos = bra_gto_block.number_of_basis_functions();");
       
    vstr.push_back("const auto bra_npgtos = bra_gto_block.number_of_primitives();");
        
    vstr.push_back("// intialize GTOs data on ket side");
        
    vstr.push_back("const auto ket_gto_exps = ket_gto_block.exponents();");
        
    vstr.push_back("const auto ket_gto_norms = ket_gto_block.normalization_factors();");
        
    vstr.push_back("const auto ket_ncgtos = ket_gto_block.number_of_basis_functions();");
       
    vstr.push_back("const auto ket_npgtos = ket_gto_block.number_of_primitives();");
    
    return vstr;
}

std::vector<std::string>
G2CFuncBodyDriver::_get_variables_def(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// define pi constant");
        
    vstr.push_back("const double fpi = mathconst::pi_value();");
    
    vstr.push_back("// set A and B centers");
    
    vstr.push_back("const auto r_a = bra_gto_block.coordinates()[bra_igto];");
    
    vstr.push_back("const auto r_b = ket_gto_block.coordinates()[ket_igto];");
    
    vstr.push_back("// set up Cartesian A coordinates");

    vstr.push_back("const auto a_xyz = r_a.coordinates();");

    vstr.push_back("const auto a_x = a_xyz[0];");

    vstr.push_back("const auto a_y = a_xyz[1];");

    vstr.push_back("const auto a_z = a_xyz[2];");
    
    vstr.push_back("// set up Cartesian B coordinates");

    vstr.push_back("const auto b_xyz = r_b.coordinates();");

    vstr.push_back("const auto b_x = b_xyz[0];");

    vstr.push_back("const auto b_y = b_xyz[1];");

    vstr.push_back("const auto b_z = b_xyz[2];");
    
    vstr.push_back("// compute overlap between A and B centers");
    
    vstr.push_back("const auto ab_x = a_x - b_x;");
    
    vstr.push_back("const auto ab_y = a_y - b_y;");
    
    vstr.push_back("const auto ab_z = a_z - b_z;");
    
    vstr.push_back("const double rab2 = ab_x * ab_x + ab_y * ab_y + ab_z * ab_z;");
    
    if (_need_boys_func(integral))
    {
        auto order = integral[0] + integral[1] + integral.integrand().shape().order();
        
        if (auto prefixes = integral.prefixes(); !prefixes.empty())
        {
            for (const auto& prefix : prefixes)
            {
                order += prefix.shape().order();
            }
        }
        
        vstr.push_back("// setup Boys function data");
        
        vstr.push_back("const CBoysFunc<" + std::to_string(order) + "> bf_table;");
    }
    
    return vstr;
}

bool
G2CFuncBodyDriver::_need_boys_func(const I2CIntegral& integral) const
{
    const auto integrand = integral.integrand();
        
    if (integrand.name() == "A") return true;
    
    return false;
}

void
G2CFuncBodyDriver::_add_loop_start(      VCodeLines&  lines,
                                   const I2CIntegral& integral) const
{
    lines.push_back({1, 0, 2, "// loop over primitives"});
    
    lines.push_back({1, 0, 1, "for (size_t i = 0; i < bra_npgtos; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "// set up primitive exponents and norms of center A"});
    
    lines.push_back({2, 0, 2, "const auto a_exp = bra_gto_exps[i * bra_ncgtos + bra_igto];"});
    
    lines.push_back({2, 0, 2, "const auto a_norm = bra_gto_norms[i * bra_ncgtos + bra_igto];"});
    
    lines.push_back({2, 0, 1, "for (size_t j = 0; j < ket_npgtos; j++)"});
    
    lines.push_back({2, 0, 1, "{"});
    
    lines.push_back({3, 0, 2, "// set up primitive exponents and norms of center B"});
    
    lines.push_back({3, 0, 2, " const auto b_exp = ket_gto_exps[j * ket_ncgtos + ket_igto];"});
    
    lines.push_back({3, 0, 2, "const auto b_norm = ket_gto_norms[j * ket_ncgtos + ket_igto];"});
    
    lines.push_back({3, 0, 2, "// compute exponential factors"});
    
    lines.push_back({3, 0, 2, "auto finv = 1.0 / (a_exp + b_exp);"});
    
    lines.push_back({3, 0, 2, "const double fzeta = a_exp * b_exp * finv;"});
    
    lines.push_back({3, 0, 2, "// compute P center coordinates"});
    
    lines.push_back({3, 0, 2, "const auto p_x = finv * (a_exp * a_x + b_exp * b_x);"});
    
    lines.push_back({3, 0, 2, "const auto p_y = finv * (a_exp * a_y + b_exp * b_y);"});
    
    lines.push_back({3, 0, 2, "const auto p_z = finv * (a_exp * a_z + b_exp * b_z);"});
    
    lines.push_back({3, 0, 2, "// compute overlap integral"});
    
    lines.push_back({3, 0, 2, "finv *= fpi;"});
    
    lines.push_back({3, 0, 2, "const auto fovl = a_norm * b_norm * finv * std::sqrt(finv) * std::exp(-fzeta * rab2);"});
    
    if (_need_distances_pa(integral))
    {
        lines.push_back({3, 0, 2, "// compute R(PA) = P - A distances"});
        
        lines.push_back({3, 0, 2, "const auto pa_x = p_x - a_x;"});
        
        lines.push_back({3, 0, 2, "const auto pa_y = p_y - a_y;"});
        
        lines.push_back({3, 0, 2, "const auto pa_z = p_z - a_z;"});
    }
    
    if (_need_distances_pb(integral))
    {
        lines.push_back({3, 0, 2, "// compute R(PB) = P - B distances"});
        
        lines.push_back({3, 0, 2, "const auto pb_x = p_x - b_x;"});
        
        lines.push_back({3, 0, 2, "const auto pb_y = p_y - b_y;"});
        
        lines.push_back({3, 0, 2, "const auto pb_z = p_z - b_z;"});
    }
    
    lines.push_back({3, 0, 2, "// compute R(PC) = P - C distances"});
    
    lines.push_back({3, 0, 2, "t2cfunc::comp_distances_pc(cart_buffer, 0, gcoords_x, gcoords_y, gcoords_z, p_x, p_y, p_z);"});
    
    if (_need_boys_func(integral))
    {
        lines.push_back({3, 0, 2, "// compute Boys function arguments"});
        
        lines.push_back({3, 0, 2, "t2cfunc::comp_boys_args(cart_buffer, 3, 0, a_exp + b_exp);"});
        
        lines.push_back({3, 0, 2, "// compute Boys function values"});
        
        lines.push_back({3, 0, 2, "bf_table.compute(cart_buffer, 4, 3);"});
    }
}

void
G2CFuncBodyDriver::_add_loop_end(      VCodeLines&    lines,
                                 const SI2CIntegrals& integrals,
                                 const I2CIntegral&   integral) const
{
   
    
    lines.push_back({3, 0, 2, "// reduce integrals"});
    
    const auto refpos = _get_position(integral, integrals, integral);
    
    const auto ncomps = integral.components<T1CPair, T1CPair>().size();
    
    std::string label = "t2cfunc::reduce(cart_buffer, ";
    
    label += std::to_string(refpos + ncomps) + ", ";
    
    label += std::to_string(refpos) + ", ";
    
    label += std::to_string(ncomps) + ");";
    
    lines.push_back({3, 0, 2, label});
    
    lines.push_back({2, 0, 1, "}"});
    
    lines.push_back({1, 0, 2, "}"});
    
    if ((integral[0] + integral[1]) > 0)
    {
        lines.push_back({1, 0, 2, "// transform integrals"});
        
        label = "t2cfunc::transform<"  + std::to_string(integral[0]);
            
        label += ", " + std::to_string(integral[1]) + ">(spher_buffer, cart_buffer, ";
        
        label += std::to_string(refpos + ncomps) + ");";
            
        lines.push_back({1, 0, 1, label});
    }
    
    std::cout << " *** (" << std::to_string(integral[0]) << "," << std::to_string(integral[1]) << ") = " << std::to_string(refpos + 2 * ncomps) << std::endl; 
}

void
G2CFuncBodyDriver::_add_call_tree(      VCodeLines&            lines,
                                  const SI2CIntegrals&         integrals,
                                  const I2CIntegral&           integral) const
{
    const int spacer = 3;
    
    lines.push_back({spacer, 0, 2, "// compute primitive integrals"});
    
    for (const auto& tint : integrals)
    {
        if (!tint.is_simple()) continue;
        
        const auto name = t2c::grid_prim_compute_func_name(tint);
            
        auto label = t2c::namespace_label(tint) + "::" + name + "(cart_buffer, ";
        
        label += _get_arguments(tint, integrals, integral);
        
        if ((tint[0] + tint[1]) == 0)
        {
            label += std::to_string(4 + tint.order()) +  ", ";
        }
       
        if (_need_distances_pa(tint))
        {
            label += "pa_x, pa_y, pa_z, ";
        }
        else
        {
            if (_need_distances_pb(tint))
            {
                label += "pb_x, pb_y, pb_z, ";
            }
        }
        
        if (_need_exponents(tint))
        {
            label += "a_exp + b_exp);";
        }
        else
        {
            if ((tint[0] + tint[1]) == 0)
            {
                label += "fovl, a_exp + b_exp);";
            }
            else
            {
                label.pop_back();
                
                label.pop_back();
                
                label += ");";
            }
        }
        
        lines.push_back({spacer, 0, 2, label});
    }
}


bool
G2CFuncBodyDriver::_need_distances_pa(const I2CIntegral& integral) const
{
    if (integral.is_simple())
    {
        return integral[0] > 0;
    }
    else
    {
        return (integral[0] + (integral.prefixes())[0].shape().order()) > 0;
    }
}

bool
G2CFuncBodyDriver::_need_distances_pb(const I2CIntegral& integral) const
{
    if (integral.is_simple())
    {
        return integral[1] > 0;
    }
    else
    {
        if (auto prefixes = integral.prefixes(); prefixes.size() == 2)
        {
            return (integral[1] + prefixes[1].shape().order()) > 0;
        }
        else
        {
            return integral[1] > 0;
        }
    }
}

bool
G2CFuncBodyDriver::_need_exponents(const I2CIntegral& integral) const
{
    int order = 0;
    
    for (const auto& prefix : integral.prefixes())
    {
        order += prefix.shape().order();
    }
    
    return (order + integral[0] + integral[1]) > 1;
}

std::string
G2CFuncBodyDriver::_get_arguments(const I2CIntegral&   integral,
                                  const SI2CIntegrals& integrals,
                                  const I2CIntegral&   ref_integral) const
{
    auto label = std::to_string(_get_position(integral, integrals, ref_integral)) + ", ";
    
    if ((integral[0] + integral[1]) > 0)
    {
        for (const auto& tint : t2c::get_integrals(integral))
        {
            label += std::to_string(_get_position(tint, integrals, ref_integral)) + ", ";
        }
    }
    
    return label;
}

size_t
G2CFuncBodyDriver::_get_position(const I2CIntegral&   integral,
                                 const SI2CIntegrals& integrals,
                                 const I2CIntegral&   ref_integral) const
{
    size_t pos = 4;
    
    auto order = ref_integral[0] + ref_integral[1] + ref_integral.integrand().shape().order() + 1;
    
    if (auto prefixes = ref_integral.prefixes(); !prefixes.empty())
    {
        for (const auto& prefix : prefixes)
        {
            order += prefix.shape().order();
        }
    }
    
    pos += order;
    
    for (const auto& tint : integrals)
    {
        if (tint == integral) return pos;
        
        pos += tint.components<T1CPair, T1CPair>().size();
    }
    
    return 0;
}
