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

#include <iostream>
#include "t2c_body.hpp"

#include "t2c_utils.hpp"

void
T2CFuncBodyDriver::write_func_body(      std::ofstream&         fstream,
                                   const SI2CIntegrals&         integrals,
                                   const I2CIntegral&           integral,
                                   const std::pair<bool, bool>& rec_form,
                                   const bool                   diagonal) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_external_data_def(integral, rec_form))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_gtos_def(diagonal))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_ket_variables_def(diagonal))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_coordinates_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_buffers_def(integrals, integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_boys_function_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_loop_start(lines, integral, diagonal);
    
    _add_ket_loop_start(lines, integral, rec_form, diagonal);
    
    _add_auxilary_integrals(lines, integrals, rec_form, false);
    
    _add_sum_loop_start(lines, integral, rec_form); 
    
    _add_auxilary_integrals(lines, integrals, rec_form, true);
    
    _add_call_tree(lines, integrals, rec_form);
    
    _add_geom_call_tree(lines, integrals, integral, rec_form);
    
    _add_sum_loop_end(lines, integrals, integral, rec_form, diagonal);
    
    _add_ket_loop_end(lines, integrals, integral, rec_form, diagonal);
    
    _add_loop_end(lines, integral, rec_form, diagonal);
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T2CFuncBodyDriver::_get_external_data_def(const I2CIntegral&           integral,
                                          const std::pair<bool, bool>& rec_form) const
{
    std::vector<std::string> vstr;
    
    if ((integral.integrand().name() == "A"))
    {
        if (rec_form.first)
        {
            vstr.push_back("// intialize external charges data");
            
            vstr.push_back("const auto charges = distributor->data();");
            
            vstr.push_back("const auto coords_x = distributor->coordinates_x();");
            
            vstr.push_back("const auto coords_y = distributor->coordinates_y();");
            
            vstr.push_back("const auto coords_z = distributor->coordinates_z();");
        }
        else
        {
            vstr.push_back("// intialize external charge data");
            
            vstr.push_back("const double charge = distributor->data()[0];");
            
            vstr.push_back("const double coord_x = distributor->coordinates_x()[0];");
            
            vstr.push_back("const double coord_y = distributor->coordinates_y()[0];");
            
            vstr.push_back("const double coord_z = distributor->coordinates_z()[0];");
        }
    }
    
    if ((integral.integrand().name() == "r"))
    {
        vstr.push_back("// intialize origin data");
            
        vstr.push_back("const double coord_x = distributor->coordinates_x()[0];");
            
        vstr.push_back("const double coord_y = distributor->coordinates_y()[0];");
            
        vstr.push_back("const double coord_z = distributor->coordinates_z()[0];");
    }
    
    if ((integral.integrand().name() == "AG"))
    {
        const auto iorder = integral.integrand().shape().order();
        
        if (rec_form.first)
        {
            if (iorder == 1)
            {
                vstr.push_back("// intialize external dipoles data");
                
                vstr.push_back("const auto dipoles = distributor->data();");
            }
            
            if (iorder == 2)
            {
                vstr.push_back("// intialize external quadrupoles data");
                
                vstr.push_back("const auto quadrupoles = distributor->data();");
            }
            
            if (iorder == 3)
            {
                vstr.push_back("// intialize external octupoles data");
                
                vstr.push_back("const auto quadrupoles = distributor->data();");
            }
           
            vstr.push_back("const auto coords_x = distributor->coordinates_x();");
            
            vstr.push_back("const auto coords_y = distributor->coordinates_y();");
            
            vstr.push_back("const auto coords_z = distributor->coordinates_z();");
        }
        else
        {
            if (iorder == 1)
            {
                vstr.push_back("// intialize external dipole data");
                
                vstr.push_back("const auto dipole = distributor->data();");
            }
            
            if (iorder == 2)
            {
                vstr.push_back("// intialize external quadrupole data");
                
                vstr.push_back("const auto quadrupole = distributor->data();");
            }
            
            if (iorder == 3)
            {
                vstr.push_back("// intialize external octupole data");
                
                vstr.push_back("const auto quadrupole = distributor->data();");
            }
            
            vstr.push_back("const double coord_x = distributor->coordinates_x()[0];");
            
            vstr.push_back("const double coord_y = distributor->coordinates_y()[0];");
            
            vstr.push_back("const double coord_z = distributor->coordinates_z()[0];");
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

        vstr.push_back("const auto gto_coords_x = gto_block.coordinates_x();");

        vstr.push_back("const auto gto_coords_y = gto_block.coordinates_y();");

        vstr.push_back("const auto gto_coords_z = gto_block.coordinates_z();");

        vstr.push_back("const auto gto_exps = gto_block.exponents();");

        vstr.push_back("const auto gto_norms = gto_block.normalization_factors();");

        vstr.push_back("const auto gto_indices = gto_block.orbital_indices();");

        vstr.push_back("const auto ncgtos = gto_block.number_of_basis_functions();");

        vstr.push_back("const auto npgtos = gto_block.number_of_primitives();");
    }
    else
    {
        vstr.push_back("// intialize GTOs data on bra side");
        
        vstr.push_back("const auto bra_gto_coords_x = bra_gto_block.coordinates_x();");
        
        vstr.push_back("const auto bra_gto_coords_y = bra_gto_block.coordinates_y();");
        
        vstr.push_back("const auto bra_gto_coords_z = bra_gto_block.coordinates_z();");
        
        vstr.push_back("const auto bra_gto_exps = bra_gto_block.exponents();");
        
        vstr.push_back("const auto bra_gto_norms = bra_gto_block.normalization_factors();");
       
        vstr.push_back("const auto bra_gto_indices = bra_gto_block.orbital_indices();");
        
        vstr.push_back("const auto bra_ncgtos = bra_gto_block.number_of_basis_functions();");

        vstr.push_back("const auto bra_npgtos = bra_gto_block.number_of_primitives();");
        
        vstr.push_back("// intialize GTOs data on ket side");
        
        vstr.push_back("const auto ket_gto_coords_x = ket_gto_block.coordinates_x();");
        
        vstr.push_back("const auto ket_gto_coords_y = ket_gto_block.coordinates_y();");
        
        vstr.push_back("const auto ket_gto_coords_z = ket_gto_block.coordinates_z();");
        
        vstr.push_back("const auto ket_gto_exps = ket_gto_block.exponents();");
        
        vstr.push_back("const auto ket_gto_norms = ket_gto_block.normalization_factors();");
       
        vstr.push_back("const auto ket_gto_indices = ket_gto_block.orbital_indices();");

        vstr.push_back("const auto ket_npgtos = ket_gto_block.number_of_primitives();");
    }
    
    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_ket_variables_def(const bool diagonal) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned 2D arrays for ket side");
    
    if (diagonal)
    {
        vstr.push_back("const auto ket_dim = gto_range[1] - gto_range[0];");
        
        vstr.push_back("const auto ket_pdim = ket_dim * npgtos;");
    }
    else
    {
        vstr.push_back("const auto ket_dim = ket_range[1] - ket_range[0];");
        
        vstr.push_back("const auto ket_pdim = ket_dim * ket_npgtos;");
    }
    
    vstr.push_back("CSimdArray<double> b_x(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> b_y(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> b_z(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> b_exps(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> b_norms(1, ket_pdim);");

    vstr.push_back(" // load GTOs data for ket side");

    if (diagonal)
    {
        vstr.push_back("b_x.replicate(gto_coords_x, gto_range, npgtos);");

        vstr.push_back("b_y.replicate(gto_coords_y, gto_range, npgtos);");

        vstr.push_back("b_z.replicate(gto_coords_z, gto_range, npgtos);");

        vstr.push_back("b_exps.load(gto_exps, gto_range, npgtos);");

        vstr.push_back("b_norms.load(gto_norms, gto_range, npgtos);");
    }
    else
    {
        vstr.push_back("b_x.replicate(ket_gto_coords_x, ket_range, ket_npgtos);");

        vstr.push_back("b_y.replicate(ket_gto_coords_y, ket_range, ket_npgtos);");

        vstr.push_back("b_z.replicate(ket_gto_coords_z, ket_range, ket_npgtos);");

        vstr.push_back("b_exps.load(ket_gto_exps, ket_range, ket_npgtos);");

        vstr.push_back("b_norms.load(ket_gto_norms, ket_range, ket_npgtos);");
    }
    
    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_coordinates_def(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (_need_center_p(integral))
    {
        vstr.push_back("// allocate aligned coordinates of P center");

        vstr.push_back("CSimdArray<double> p_x(1, ket_pdim);");

        vstr.push_back("CSimdArray<double> p_y(1, ket_pdim);");

        vstr.push_back("CSimdArray<double> p_z(1, ket_pdim);");
    }
    

    vstr.push_back("// allocate aligned distances R(AB) = A - B");

    vstr.push_back("CSimdArray<double> ab_x(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> ab_y(1, ket_pdim);");

    vstr.push_back("CSimdArray<double> ab_z(1, ket_pdim);");
    
    if (integral[1] > 0)
    {
        vstr.push_back("// allocate aligned distances R(PB) = P - B");
        
        vstr.push_back("CSimdArray<double> pb_x(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> pb_y(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> pb_z(1, ket_pdim);");
    }
    
    if (integral[0] > 0)
    {
        vstr.push_back("// allocate aligned distances R(PA) = P - A");
        
        vstr.push_back("CSimdArray<double> pa_x(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> pa_y(1, ket_pdim);");
        
        vstr.push_back("CSimdArray<double> pa_z(1, ket_pdim);");
    }

    if (_need_distances_pc(integral))
    {
        vstr.push_back("// allocate aligned distances R(PC) = P - C");

        vstr.push_back("CSimdArray<double> pc_x(1, ket_pdim);");

        vstr.push_back("CSimdArray<double> pc_y(1, ket_pdim);");

        vstr.push_back("CSimdArray<double> pc_z(1, ket_pdim);");
    }
    
    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_buffers_def(const SI2CIntegrals& integrals,
                                    const I2CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned primitive integrals");
    
    size_t icomps = 0;
    
    for (const auto& tint : integrals)
    {
        icomps += tint.components<T1CPair, T1CPair>().size();
    }
    
    auto label = "CSimdArray<double> pbuffer(" + std::to_string(icomps) + ", ket_pdim);";
    
    vstr.push_back(label);
    
    vstr.push_back("// allocate aligned contracted integrals");
    
    icomps = integral.components<T1CPair, T1CPair>().size();
    
    label = "CSimdArray<double> cbuffer(" + std::to_string(icomps) + ", ket_dim);";
    
    vstr.push_back(label);
    
    if ((integral[0] + integral[1]) > 0)
    {
        const auto angpair = std::array<int, 2>({integral[0], integral[1]});
        
        icomps = t2c::number_of_spherical_components(angpair);
        
        icomps *= integral.integrand().components().size();
        
        if (const auto prefixes = integral.prefixes(); !prefixes.empty())
        {
            icomps *= make_components<OperatorComponent>(prefixes).size();
        }
        
        label = "CSimdArray<double> sbuffer(" + std::to_string(icomps) + ", ket_dim);";
        
        vstr.push_back(label);
    }
    
    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_boys_function_def(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if (_need_boys_func(integral))
    {
        auto order = integral[0] + integral[1] + integral.integrand().shape().order();
        
        vstr.push_back("// setup Boys function data");
        
        vstr.push_back("const CBoysFunc<" + std::to_string(order) + "> bf_table;");
        
        vstr.push_back("CSimdArray<double> bf_args(1, ket_pdim);");

        vstr.push_back("CSimdArray<double> bf_values(" + std::to_string(order + 1) + ", ket_pdim);");
    }
    
    return vstr;
}

void
T2CFuncBodyDriver::_add_loop_start(      VCodeLines&  lines,
                                   const I2CIntegral& integral,
                                   const bool         diagonal) const
{
    lines.push_back({1, 0, 2, "// loop over contracted GTOs on bra side"});
   
    if (diagonal)
    {
        lines.push_back({1, 0, 1, "for (auto i = gto_range[0]; i < gto_range[1]; i++)"});
    }
    else
    {
        lines.push_back({1, 0, 1, "for (auto i = bra_range[0]; i < bra_range[1]; i++)"});
    }
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "cbuffer.zero();"});

    if ((integral[0] + integral[1]) > 0)
    {
        lines.push_back({2, 0, 2, "sbuffer.zero();"});
    }

    if (diagonal && (integral[0] == integral[1]))
    {
        lines.push_back({2, 0, 2, "const auto a_x = gto_coords_x[i];"});

        lines.push_back({2, 0, 2, "const auto a_y = gto_coords_y[i];"});

        lines.push_back({2, 0, 2, "const auto a_z = gto_coords_z[i];"});
    }
    else
    {
        lines.push_back({2, 0, 2, "const auto a_x = bra_gto_coords_x[i];"});

        lines.push_back({2, 0, 2, "const auto a_y = bra_gto_coords_y[i];"});

        lines.push_back({2, 0, 2, "const auto a_z = bra_gto_coords_z[i];"});
    }

    lines.push_back({2, 0, 2, "t2cfunc::comp_distances_ab(ab_x[0], ab_y[0], ab_z[0], a_x, a_y, a_z, b_x[0], b_y[0], b_z[0], ket_pdim);"});
}

void
T2CFuncBodyDriver::_add_loop_end(      VCodeLines&  lines,
                                 const I2CIntegral& integral,
                                 const std::pair<bool, bool>& rec_form,
                                 const bool         diagonal) const
{
    std::string label;
    
    if ((integral[0] + integral[1]) > 0)
    {
        label = "t2cfunc::transform<"  + std::to_string(integral[0]);
            
        label += ", " + std::to_string(integral[1]) + ">(sbuffer, cbuffer);";
            
        lines.push_back({2, 0, 2, label});
    }
    
    if ((integral[0] + integral[1]) > 0)
    {
        label = "distributor->distribute(sbuffer, ";
    }
    else
    {
        label = "distributor->distribute(cbuffer, ";
    }
    
    if (diagonal)
    {
        label += "gto_indices, ";
    }
    else
    {
        label += "bra_gto_indices, ket_gto_indices, ";
    }
            
    label += std::to_string(integral[0]) + ", ";
            
    if (!diagonal)
    {
        label += std::to_string(integral[1]) + ", ";
    }
    
    if (diagonal)
    {
        label += "i, gto_range);";
    }
    else
    {
        label += "i, ket_range);";
    }
    
    lines.push_back({2, 0, 1, label});
   
    lines.push_back({1, 0, 1, "}"});
}

void
T2CFuncBodyDriver::_add_ket_loop_start(      VCodeLines&            lines,
                                       const I2CIntegral&           integral,
                                       const std::pair<bool, bool>& rec_form,
                                       const bool                   diagonal) const
{
    if (diagonal)
    {
        lines.push_back({2, 0, 1, "for (int j = 0; j < npgtos; j++)"});
    }
    else
    {
        lines.push_back({2, 0, 1, "for (int j = 0; j < bra_npgtos; j++)"});
    }
    
    lines.push_back({2, 0, 1, "{"});
    
    if (diagonal && (integral[0] == integral[1]))
    {
        lines.push_back({3, 0, 2, "const auto a_exp = gto_exps[j * ncgtos + i];"});
            
        lines.push_back({3, 0, 2, "const auto a_norm = gto_norms[j * ncgtos + i];"});
    }
    else
    {
        lines.push_back({3, 0, 2, "const auto a_exp = bra_gto_exps[j * bra_ncgtos + i];"});
            
        lines.push_back({3, 0, 2, "const auto a_norm = bra_gto_norms[j * bra_ncgtos + i];"});
    }

    if (_need_center_p(integral))
    {
        lines.push_back({3, 0, 2, "t2cfunc::comp_coordinates_p(p_x[0], p_y[0], p_z[0], a_x, a_y, a_z, b_x[0], b_y[0], b_z[0], a_exp, b_exps[0], ket_pdim);"});
    }

    if (integral[0] > 0)
    {
        if (_need_center_p(integral))
        {
            lines.push_back({3, 0, 2, "t2cfunc::comp_distances_pa(pa_x[0], pa_y[0], pa_z[0], p_x[0], p_y[0], p_z[0], a_x, a_y, a_z, ket_pdim);"});
        }
        else
        {
            lines.push_back({3, 0, 2, "t2cfunc::comp_distances_pa(pa_x[0], pa_y[0], pa_z[0], ab_x[0], ab_y[0], ab_z[0], a_exp, b_exps[0], ket_pdim);"});
        }
    }

    if (integral[1] > 0)
    {
        if (_need_center_p(integral))
        {
            lines.push_back({3, 0, 2, "t2cfunc::comp_distances_pb(pb_x[0], pb_y[0], pb_z[0], p_x[0], p_y[0], p_z[0], b_x[0], b_y[0], b_z[0], ket_pdim);"});
        }
        else
        {
            lines.push_back({3, 0, 2, "t2cfunc::comp_distances_pb(pb_x[0], pb_y[0], pb_z[0], ab_x[0], ab_y[0], ab_z[0], a_exp, b_exps[0], ket_pdim);"});
        }
        
    }

    if (_need_distances_pc(integral) && !rec_form.first)
    {
        lines.push_back({3, 0, 2, "t2cfunc::comp_distances_pc(pc_x[0], pc_y[0], pc_z[0], p_x[0], p_y[0], p_z[0], coord_x, coord_y, coord_z, ket_pdim);"});
    }
}

void
T2CFuncBodyDriver::_add_ket_loop_end(      VCodeLines&            lines,
                                     const SI2CIntegrals&         integrals,
                                     const I2CIntegral&           integral,
                                     const std::pair<bool, bool>& rec_form,
                                     const bool                   diagonal) const
{
    if (!rec_form.first)
    {
        std::string label = "t2cfunc::reduce(cbuffer, pbuffer, ";
        
        label += std::to_string(_get_position(integral, integrals)) + ", ";
        
        if (diagonal)
        {
            label += "ket_dim, npgtos);";
        }
        else
        {
            label += "ket_dim, ket_npgtos);";
        }
        
        lines.push_back({3, 0, 1, label});
    }
    
    lines.push_back({2, 0, 2, "}"});
}

void
T2CFuncBodyDriver::_add_sum_loop_start(      VCodeLines&            lines,
                                       const I2CIntegral&           integral,
                                       const std::pair<bool, bool>& rec_form) const
{
    if (rec_form.first)
    {
        if (integral.integrand().name() == "A")
        {
            lines.push_back({3, 0, 2, "const auto ncenters = static_cast<int>(charges.size());"});
        }
        
        if (integral.integrand().name() == "AG")
        {
            lines.push_back({3, 0, 2, "const auto ncenters = static_cast<int>(coords_x.size());"});
        }
        
        lines.push_back({3, 0, 1, "for (int k = 0; k < ncenters; k++)"});
       
        lines.push_back({3, 0, 1, "{"});
       
        if (_need_distances_pc(integral))
        {
            lines.push_back({4, 0, 2, "t2cfunc::comp_distances_pc(pc_x[0], pc_y[0], pc_z[0], p_x[0], p_y[0], p_z[0], coords_x[k], coords_y[k], coords_z[k], ket_pdim);"});
        }
        
        if (_need_boys_func(integral))
        {
            lines.push_back({4, 0, 2, "t2cfunc::comp_boys_args(bf_args, pc_x[0], pc_y[0], pc_z[0], a_exp, b_exps[0]);"});
            
            lines.push_back({4, 0, 2, "bf_table.compute(bf_values, bf_args);"});
        }
    }
}

void
T2CFuncBodyDriver::_add_sum_loop_end(      VCodeLines&            lines,
                                     const SI2CIntegrals&         integrals,
                                     const I2CIntegral&           integral,
                                     const std::pair<bool, bool>& rec_form,
                                     const bool                   diagonal) const
{
    if (rec_form.first)
    {
        std::string label = "t2cfunc::reduce(cbuffer, pbuffer, ";
        
        label += std::to_string(_get_position(integral, integrals)) + ", ";
        
        if (integral.integrand().name() == "A")
        {
            label += "charges[k], ";
        }
        
        if (integral.integrand().name() == "AG")
        {
            const auto iorder = integral.integrand().shape().order();
            
            if (iorder == 1) label += "dipoles, 3, k, ";
            
            if (iorder == 2) label += "quadrupoles, 3, k, ";
        }
        
        if (diagonal)
        {
            label += "ket_dim, npgtos);";
        }
        else
        {
            label += "ket_dim, ket_npgtos);";
        }
        
        lines.push_back({4, 0, 1, label});
        
        lines.push_back({3, 0, 1, "}"});
    }
}

void
T2CFuncBodyDriver::_add_auxilary_integrals(      VCodeLines&            lines,
                                           const SI2CIntegrals&         integrals,
                                           const std::pair<bool, bool>& rec_form,
                                           const bool                   in_sum_loop) const
{
    const auto spacer = (rec_form.first) ? 4 : 3;
    
    for (const auto& tint : integrals)
    {
        if (!tint.is_simple()) continue;
        
        if ((tint[0] == 0) && (tint[1] == 0))
        {
            if ((tint.integrand().name() == "1") && (!in_sum_loop))
            {
                lines.push_back({3, 0, 2, "ovlrec::comp_prim_overlap_ss(pbuffer, " + std::to_string(_get_position(tint, integrals)) + ", ab_x[0], ab_y[0], ab_z[0], a_exp, b_exps[0], a_norm, b_norms[0]);"});
            }
            
            if ((tint.integrand().name() == "T") && (!in_sum_loop))
            {
                const auto sint = tint.replace(Operator("1"));
                
                const auto label = std::to_string(_get_position(tint, integrals)) + ", " + std::to_string(_get_position(sint, integrals));
                
                lines.push_back({3, 0, 2, "kinrec::comp_prim_kinetic_energy_ss(pbuffer, " + label + ", ab_x[0], ab_y[0], ab_z[0], a_exp, b_exps[0]);"});
            }

            if ((tint.integrand().name() == "r") && (!in_sum_loop))
            {
                const auto sint = tint.replace(Operator("1"));
                
                const auto label = std::to_string(_get_position(tint, integrals)) + ", " + std::to_string(_get_position(sint, integrals));
                
                lines.push_back({3, 0, 2, "diprec::comp_prim_dipole_ss(pbuffer, " + label + ", pc_x[0], pc_y[0], pc_z[0]);"});
            }

            if ((tint.integrand().name() == "p") && (!in_sum_loop))
            {
                lines.push_back({3, 0, 2, "linmomrec::comp_prim_dipole_ss(prim_buffer_dip_ss, prim_buffer_ovl_ss);"});
            }
            
            if ((tint.integrand().name() == "A"))
            {
                const auto sint = tint.replace(Operator("1"));
                
                if (rec_form.first && in_sum_loop)
                {
                    const auto label = std::to_string(_get_position(tint, integrals)) + ", " + std::to_string(_get_position(sint, integrals));
                    
                    lines.push_back({spacer, 0, 2, "npotrec::comp_prim_nuclear_potential_ss(pbuffer, " + label + ", bf_values[" + std::to_string(tint.order()) + "], a_exp, b_exps[0]);"});
                }
            }
            
            if ((tint.integrand().name() == "AG"))
            {
                const auto iorder = tint.integrand().shape().order();
                
                if (rec_form.first && in_sum_loop)
                {
                    if (iorder == 1)
                    {
                        auto xint = tint.replace(Operator("A"));
                        
                        auto label = std::to_string(_get_position(tint, integrals)) + ", ";
                        
                        xint.set_order(tint.order() + 1);
                        
                        label += std::to_string(_get_position(xint, integrals));
                        
                        lines.push_back({spacer, 0, 2, "npotrec::comp_prim_nuclear_potential_geom010_ss(pbuffer, " + label + ", pc_x[0], pc_y[0], pc_z[0], a_exp, b_exps[0]);"});
                    }
                    
                    if (iorder == 2)
                    {
                        auto xint = tint.replace(Operator("A"));
                        
                        auto label = std::to_string(_get_position(tint, integrals)) + ", ";
                        
                        label += std::to_string(_get_position(xint, integrals)) + ", ";
                        
                        xint.set_order(tint.order() + 2);
                        
                        label += std::to_string(_get_position(xint, integrals));
                        
                        lines.push_back({spacer, 0, 2, "npotrec::comp_prim_nuclear_potential_geom020_ss(pbuffer, " + label + ", pc_x[0], pc_y[0], pc_z[0], a_exp, b_exps[0]);"});
                    }
                }
            }
            
            // TODO: other integrals...
        }
    }
}

// May need to change this for new integral situations
void
T2CFuncBodyDriver::_add_call_tree(      VCodeLines&  lines,
                                  const SI2CIntegrals& integrals,
                                  const std::pair<bool, bool>& rec_form) const
{
    const int spacer = (rec_form.first) ? 4 : 3;
    
    for (const auto& tint : integrals)
    {
        if (!tint.is_simple()) continue;
        
        if ((tint[0] != 0) || (tint[1] != 0))
        {
            const auto name = t2c::prim_compute_func_name(tint);
            
            auto label = t2c::namespace_label(tint) + "::" + name + "(pbuffer, ";
            
            label += _get_arguments(tint, integrals);
            
            if (tint[0] > 0)
            {
                if ((tint[0] == 1) && (tint[1] == 0) && (tint.integrand().name() != "T") && (tint.integrand().name() != "A") && (tint.integrand().name() != "AG"))
                {
                    label += "pa_x[0], pa_y[0], pa_z[0]";
                }
                else
                {
                    label += "pa_x[0], pa_y[0], pa_z[0], ";
                }
            }
            
            if ((tint[1] > 0) && (tint[0] == 0))
            {
                if ((tint[1] == 1) && (tint[0] == 0) && (tint.integrand().name() != "T") && (tint.integrand().name() != "A") && (tint.integrand().name() != "AG"))
                {
                    label += "pb_x[0], pb_y[0], pb_z[0]";
                }
                else
                {
                    label += "pb_x[0], pb_y[0], pb_z[0], ";
                }
            }
            
            if (tint.integrand().name() == "A")
            {
                if ((tint[0] + tint[1]) > 1)
                {
                    label += "pc_x[0], pc_y[0], pc_z[0], ";
                }
                
                if ((tint[0] + tint[1]) == 1)
                {
                    label += "pc_x[0], pc_y[0], pc_z[0]";
                }
            }

            if (tint.integrand().name() == "AG")
            {
                if ((tint[0] + tint[1]) > 1)
                {
                    label += "pc_x[0], pc_y[0], pc_z[0], ";
                }

                if ((tint[0] + tint[1]) == 1)
                {
                    label += "pc_x[0], pc_y[0], pc_z[0]";
                }
            }

            if (((tint[0] + tint[1]) > 1) || (tint.integrand().name() == "T"))
            {
                label += "a_exp, b_exps[0]";
            }
           else if (((tint[0] + tint[1]) > 0) && (tint.integrand().name() == "r"))
            {
                label += ", a_exp, b_exps[0]";
            }
            
            label += ");";
            
            lines.push_back({spacer, 0, 2, label});
        }
    }
}

void
T2CFuncBodyDriver::_add_geom_call_tree(      VCodeLines&            lines,
                                       const SI2CIntegrals&         integrals,
                                       const I2CIntegral&           integral,
                                       const std::pair<bool, bool>& rec_form) const
{

    if (integral.prefixes().size() > 0)
    {
        const int spacer = (rec_form.first) ? 4 : 3;

        auto geom_order = 0;

        if (integral.prefixes().size() == 1)
        {

            if (integral.prefixes()[0].shape().order() == 0)
            {
                geom_order = geom_order + integral.prefixes()[1].shape().order();
            }
            else
            {
                geom_order = integral.prefixes()[0].shape().order();
            }
        }

        if (integral.prefixes().size() == 2)
        {
           geom_order = integral.prefixes()[0].shape().order() + integral.prefixes()[1].shape().order();
        }

        for (int i = 1; i <= geom_order; i++)
        {
            for (const auto& tint : integrals)
            {
                if (tint.prefixes().size() > 0)
                {
                    auto tint_geom_order = 0;

                    if (tint.prefixes().size() == 1)
                    {

                        if (tint.prefixes()[0].shape().order() == 0)
                        {
                            tint_geom_order = tint_geom_order + tint.prefixes()[1].shape().order();
                        }
                        else
                        {
                            tint_geom_order = tint.prefixes()[0].shape().order();
                        }
                    }

                    if (tint.prefixes().size() == 2)
                    {
                       tint_geom_order = tint.prefixes()[0].shape().order() + tint.prefixes()[1].shape().order();
                    }

                    if (tint_geom_order == i)
                    {
                        const auto name = t2c::prim_compute_func_name(tint);
                            
                        auto label = t2c::namespace_label(tint) + "::" + name + "(";
                            
                        label += _get_arguments(tint);
                        
                        label += "a_exp, b_exps[0]);";
                    
                        lines.push_back({spacer, 0, 2, label});
                    }
                }
            }
        }
    }
}

std::string
T2CFuncBodyDriver::_get_arguments(const I2CIntegral& integral) const
{
    std::string label = t2c::get_buffer_label(integral, {"prim"}) + ", ";;
   
    for (const auto& tint : t2c::get_integrals(integral))
    {
        label += t2c::get_buffer_label(tint, {"prim"}) + ", ";
    }
    
    return label;
}

std::string
T2CFuncBodyDriver::_get_arguments(const I2CIntegral&   integral,
                                  const SI2CIntegrals& integrals) const
{
    auto label = std::to_string(_get_position(integral, integrals)) + ", ";
    
    for (const auto& tint : t2c::get_integrals(integral))
    {
        label += std::to_string(_get_position(tint, integrals)) + ", ";
    }
    
    return label;
}

size_t
T2CFuncBodyDriver::_get_position(const I2CIntegral&   integral,
                                 const SI2CIntegrals& integrals) const
{
    size_t pos = 0;
    
    for (const auto& tint : integrals)
    {
        if (tint == integral) return pos;
        
        pos += tint.components<T1CPair, T1CPair>().size();
    }
    
    return 0;
}

bool
T2CFuncBodyDriver::_need_center_p(const I2CIntegral& integral) const
{
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "r") return true;
    
    if (integrand.name() == "A") return true;
    
    if (integrand.name() == "AG") return true;
    
    return false;
}

bool
T2CFuncBodyDriver::_need_distances_pc(const I2CIntegral& integral) const
{
    const auto integrand = integral.integrand();
    
    if (integrand.name() == "r") return true;
    
    if (integrand.name() == "A") return true;
    
    if (integrand.name() == "AG") return true;
    
    return false;
}

bool
T2CFuncBodyDriver::_need_boys_func(const I2CIntegral& integral) const
{
    const auto integrand = integral.integrand();
        
    if (integrand.name() == "A") return true;
    
    if (integrand.name() == "AG") return true;
    
    return false;
}
