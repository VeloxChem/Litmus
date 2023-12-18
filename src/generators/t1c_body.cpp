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

#include "t1c_body.hpp"

#include "t2c_center_driver.hpp"
#include "t2c_utils.hpp"
#include "spherical_momentum.hpp"

#include <iostream>

void
T1CFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                   const int            angmom,
                                   const int            gdrv) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_angmom_def(angmom))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_gtos_def(angmom, gdrv))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_loop_body(lines, angmom, gdrv);
        
    lines.push_back({1, 0, 1, "return gto_values;"});
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

I2CIntegral
T1CFuncBodyDriver::_get_integral(const int angmom,
                                 const int gdrv) const
{
    // bra and ket sides
    
    const auto bra = I1CPair("GA", angmom);
    
    const auto ket = I1CPair("GB", 0);
    
    // prefixes of integral bra, ket order
    
    VOperators prefixes;
    
    prefixes.push_back(Operator("d/dR", Tensor(gdrv)));

    return I2CIntegral(bra, ket, Operator("1"), 0, prefixes);
}

std::vector<std::string>
T1CFuncBodyDriver::_get_angmom_def(const int angmom) const
{
    std::vector<std::string> vstr;
    
    if (angmom > 1)
    {
        const auto ang_mom = SphericalMomentum(angmom);
            
        vstr.push_back("// spherical transformation factors");
       
        for (const auto& label : ang_mom.get_factors(angmom))
        {
            vstr.push_back("const double " + label + ";");
        }
    }
    
    return vstr;
}

std::vector<std::string>
T1CFuncBodyDriver::_get_gtos_def(const int angmom,
                                 const int gdrv) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// set up GTO values storage");

    vstr.push_back("const auto nrows = mathfunc::countSignificantElements(gtos_mask);");

    vstr.push_back("const auto ncols = static_cast<int64_t>(grid_coords_x.size());");

    vstr.push_back("// set up GTOs data");

    vstr.push_back("const auto gto_exps = gto_block.getExponents();");

    vstr.push_back("const auto gto_norms = gto_block.getNormalizationFactors();");

    vstr.push_back("const auto gto_coords = gto_block.getCoordinates();");

    vstr.push_back("// set up grid data");

    vstr.push_back("auto g_x = grid_coords_x.data();");

    vstr.push_back("auto g_y = grid_coords_y.data();");

    vstr.push_back("auto g_z = grid_coords_z.data();");

    vstr.push_back("// set GTOs block dimensions");

    vstr.push_back("const auto ncgtos = gto_block.getNumberOfBasisFunctions();");

    vstr.push_back("const auto npgtos = gto_block.getNumberOfPrimitives();");
    
    vstr.push_back("// set storage matrix");
    
    const auto geom_comps = Tensor(gdrv).components();
    
    auto label = std::to_string(2 * angmom + 1) + " * ";
    
    vstr.push_back("auto gto_values = matfunc::makeMatrix(" + std::to_string(angmom) + ", " + label + "nrows, ncols);");

    vstr.push_back("gto_values.zero();");
    
    vstr.push_back("// set submatrices");

    for (size_t i = 0; i < geom_comps.size(); i++)
    {
        label = "auto submat_" + geom_comps[i].label() + " = ";
        
        label += "gto_values.getSubMatrix({" + std::to_string(gdrv) + ", ";
        
        label += std::to_string(i) + "});";
        
        vstr.push_back(label);
    }
    
    vstr.push_back("// set temporary buffer for contracted GTOs");

    for (size_t i = 0; i < geom_comps.size(); i++)
    {
        vstr.push_back("std::vector<double> buffer_" + geom_comps[i].label() + "(ncols);");
    }
    
    for (size_t i = 0; i < geom_comps.size(); i++)
    {
        vstr.push_back("auto ptr_buffer_" + geom_comps[i].label() +  " = buffer_" + geom_comps[i].label() + ".data();");
    }
    
    return vstr;
}

void
T1CFuncBodyDriver::_add_loop_body(      VCodeLines&  lines,
                                  const int          angmom,
                                  const int          gdrv) const
{
    const auto tint = _get_integral(angmom, gdrv);
    
    const auto gten = Tensor(angmom);
    
    const auto gcomps = gten.components();
    
    for (size_t i = 0; i < gcomps.size(); i++)
    {
        auto label = gten.label();
        
        label += (gten.order() > 0) ? "_" + gcomps[i].label() : "";
        
        lines.push_back({1, 0, 2, "// compute geometrical derivatives for " + label + " type GTOs" });
        
        if (i == 0)
        {
            lines.push_back({1, 0, 2, "int64_t irow = 0;" });
        }
        else
        {
            lines.push_back({1, 0, 2, "irow = 0;" });
        }
        
        lines.push_back({1, 0, 1, "for (int64_t i = 0; i < ncgtos; i++)" });
        
        lines.push_back({1, 0, 1, "{" });
       
        lines.push_back({2, 0, 1, "if (gtos_mask[i] == 1)"});
    
        lines.push_back({2, 0, 1, "{"});
        
        lines.push_back({3, 0, 2, "// set up GTO coordinates"});

        lines.push_back({3, 0, 2, "const auto r_x = gto_coords[i][0];"});

        lines.push_back({3, 0, 2, "const auto r_y = gto_coords[i][1];"});

        lines.push_back({3, 0, 2, "const auto r_z = gto_coords[i][2];"});

        lines.push_back({3, 0, 2, "// compute GTO values on grid"});
        
        for (const auto& gcomp : Tensor(gdrv).components())
        {
            lines.push_back({3, 0, 2, "mathfunc::zero(buffer_" + gcomp.label() + ");"});
        }
        
        lines.push_back({3, 0, 1, "for (int64_t j = 0; j < npgtos; j++)"});
            
        lines.push_back({3, 0, 1, "{"});
                
        lines.push_back({4, 0, 2, "const auto tbe_0 = gto_exps[j * ncgtos + i];"});

        lines.push_back({4, 0, 2, "const auto fnorm = gto_norms[j * ncgtos + i];"});

        lines.push_back({4, 0, 1, "#pragma omp simd"});
                            
        lines.push_back({4, 0, 1, "for (int64_t k = 0; k < ncols; k++)"});
                                
        lines.push_back({4, 0, 1, "{"});
        
        lines.push_back({5, 0, 2, "const auto gr_x = g_x[k] - r_x;"});

        lines.push_back({5, 0, 2, "const auto gr_y = g_y[k] - r_y;"});

        lines.push_back({5, 0, 2, "const auto gr_z = g_z[k] - r_z;"});

        lines.push_back({5, 0, 2, "const auto fss = fnorm * std::exp(-tbe_0 * (gr_x * gr_x + gr_y * gr_y + gr_z * gr_z));"});
        
        _add_simd_code(lines, tint, gcomps[i], gdrv);
        
        lines.push_back({4, 0, 1, "}"});
        
        lines.push_back({3, 0, 2, "}"});
        
        _add_distribution_code(lines, gcomps[i], gdrv); 
        
        lines.push_back({3, 0, 1, "irow++;" });
        
        lines.push_back({2, 0, 1, "}" });

        lines.push_back({1, 0, 2, "}" });
    }
}

R2Group
T1CFuncBodyDriver::_generate_integral_group(const VT2CIntegrals& components) const
{
    T2CCenterDriver t2c_geom_drv;
        
    auto rgroup = t2c_geom_drv.create_recursion(components);
    
    rgroup.simplify();
    
    return rgroup;
}

VT2CIntegrals
T1CFuncBodyDriver::_select_integral_components(const TensorComponent& component,
                                               const I2CIntegral&     integral) const
{
    VT2CIntegrals tcomps;
    
    for (const auto& tcomp : integral.components<T1CPair, T1CPair>())
    {
        if (tcomp.bra().shape() == component) tcomps.push_back(tcomp);
    }
        
    return tcomps;
}

void
T1CFuncBodyDriver::_add_simd_code(      VCodeLines&      lines,
                                  const I2CIntegral&     integral,
                                  const TensorComponent& angcomp,
                                  const int              gdrv) const
{
    const auto tints = _select_integral_components(angcomp, integral);
    
    const auto rgroup = _generate_integral_group(tints);
    
    for (size_t i = 0; i < rgroup.expansions(); i++)
    {
        _add_simd_line(lines, rgroup[i]);
    }
}

void
T1CFuncBodyDriver::_add_simd_line(      VCodeLines&  lines,
                                  const R2CDist&     rdist) const
{
    const auto prefix = rdist.root().integral().prefixes()[0];
    
    auto label = "ptr_buffer_" + prefix.shape().label() + "[k] += ";
    
    const auto nterms = rdist.terms();
    
    if (nterms > 1) label += "(";
    
    for (size_t i = 0; i < nterms; i++)
    {
        const auto tcomp = rdist[i][0];
        
        const auto tlabel = _polynomial_string(tcomp);
        
        const auto flabel = t2c::get_factor_label(rdist[i], i == 0);
        
        if (flabel.empty())
        {
            if (tlabel.empty())
            {
                label += "1.0";
            }
        }
        else
        {
            label += flabel;
        }
        
        if (tlabel.size() > 0)
        {
            if (!flabel.empty())
            {
                label += " * ";
            }
        }
        
        label += tlabel;
    }
    
    if (nterms > 1) label += ")";
    
    label += " * fss;";
    
    lines.push_back({5, 0, 2, label});
}

std::string
T1CFuncBodyDriver::_polynomial_string(const TensorComponent& component) const
{
    std::string label;
    
    for (int i = 0; i < component['x']; i++)
    {
        label += " * gr_x";
    }
    
    for (int i = 0; i < component['y']; i++)
    {
        label += " * gr_y";
    }
    
    for (int i = 0; i < component['z']; i++)
    {
        label += " * gr_z";
    }
    
    return label.erase(0, 3);
}

void
T1CFuncBodyDriver::_add_distribution_code(      VCodeLines&      lines,
                                          const TensorComponent& component,
                                          const int              gdrv) const
{
    lines.push_back({3, 0, 2, "// distribute GTO values into submatrices"});
    
    const auto blabels = t2c::tensor_components(Tensor(gdrv), "buffer");
    
    const auto mlabels = t2c::tensor_components(Tensor(gdrv), "submat");
    
    const auto ang_mom = SphericalMomentum(Tensor(component).order());
                               
    const auto index = t2c::tensor_component_index(component);
    
    const auto ang_pairs = ang_mom.select_pairs(index);

    for (int i = 0; i < blabels.size(); i++)
    {
        for (const auto& ang_pair : ang_pairs)
        {
            std::string flabel = (ang_pair.second == "1.0")  ? "" : ", "  + ang_pair.second;
            
            std::string ilabel = "irow";
            
            if (ang_pair.first > 1)
            {
                ilabel = std::to_string(ang_pair.first) + " * nrows + " + ilabel;
            }
            
            lines.push_back({3, 0, 2, "gtoval::distribute(" + mlabels[i] + ", " + blabels[i] + flabel + ", "  + ilabel + ");"});
        }
    }
}
