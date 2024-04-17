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
    
    
    lines.push_back({0, 0, 1, "}"});
    
    ost::write_code_lines(fstream, lines);
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
        
        vstr.push_back("const auto ket_ncgtos = ket_gto_block.number_of_basis_functions();");

        vstr.push_back("const auto ket_npgtos = ket_gto_block.number_of_primitives();");
    }
    
    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_ket_variables_def(const bool diagonal) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned 2D arrays for ket side");
    
    vstr.push_back("const auto ket_dim = ket_indices[1] - ket_indices[0];");
    
    if (diagonal)
    {
        vstr.push_back("const auto ket_pdim = ket_dim * npgtos;");
    }
    else
    {
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
        vstr.push_back("b_x.replicate(gto_coords_x, ket_indices, npgtos);");

        vstr.push_back("b_y.replicate(gto_coords_y, ket_indices, npgtos);");

        vstr.push_back("b_z.replicate(gto_coords_z, ket_indices, npgtos);");

        vstr.push_back("b_exps.load(gto_exps, ket_indices, npgtos);");

        vstr.push_back("b_norms.load(gto_norms, ket_indices, npgtos);");
    }
    else
    {
        vstr.push_back("b_x.replicate(ket_gto_coords_x, ket_indices, ket_npgtos);");

        vstr.push_back("b_y.replicate(ket_gto_coords_y, ket_indices, ket_npgtos);");

        vstr.push_back("b_z.replicate(ket_gto_coords_z, ket_indices, ket_npgtos);");

        vstr.push_back("b_exps.load(ket_gto_exps, ket_indices, ket_npgtos);");

        vstr.push_back("b_norms.load(ket_gto_norms, ket_indices, ket_npgtos);");
    }
    
    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_coordinates_def(const I2CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
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
    
    return vstr;
}

std::vector<std::string>
T2CFuncBodyDriver::_get_buffers_def(const SI2CIntegrals& integrals,
                                    const I2CIntegral&   integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// allocate aligned primitive integrals");
    
    for (const auto& tint : integrals)
    {
        std::string label = "CSimdArray<double> ";
        
        label += t2c::get_buffer_label(tint, "prim");
        
        const auto angpair = std::array<int, 2>({tint[0], tint[1]});
        
        const auto tcomps = t2c::number_of_cartesian_components(angpair);
        
        label += "(" + std::to_string(tcomps) + ", ket_pdim);";
        
        vstr.push_back(label);
    }
    
    vstr.push_back("// allocate aligned contracted integrals");
    
    const auto angpair = std::array<int, 2>({integral[0], integral[1]});
    
    auto icomps = t2c::number_of_cartesian_components(angpair);
    
    std::string label = "CSimdArray<double> ";
    
    label += t2c::get_buffer_label(integral, "cart");
    
    label += "(" + std::to_string(icomps) + ", ket_dim);";
    
    vstr.push_back(label);
    
    if ((integral[0] > 1) || (integral[1] > 1))
    {
        icomps = t2c::number_of_spherical_components(angpair);
        
        std::string label = "CSimdArray<double> ";
        
        label += t2c::get_buffer_label(integral, "spher");
        
        label += "(" + std::to_string(icomps) + ", ket_dim);";
        
        vstr.push_back(label);
    }
    
    return vstr;
}
