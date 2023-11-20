#include "t4c_full_body.hpp"

#include "string_formater.hpp"
#include "spherical_momentum.hpp"
#include "t4c_utils.hpp"
#include "t2c_utils.hpp"
#include "t4c_hrr_eri_driver.hpp"

void
T4CFullFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                       const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_angmom_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_gtos_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_vars_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_batches_def())
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_batches_loop_start(lines);
    
    _add_batches_loop_body(lines, integral);
    
    _add_batches_loop_end(lines, integral);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

void
T4CFullFuncBodyDriver::write_comp_func_body(      std::ofstream& fstream,
                                            const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_angmom_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_comp_gtos_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    for (const auto& label : _get_vars_def(integral))
    {
        lines.push_back({1, 0, 2, label});
    }
    
    _add_comp_batches_loop_start(lines);
    
    _add_comp_batches_loop_body(lines, integral);
    
    //_add_batches_loop_end(lines, integral);
    
    lines.push_back({1, 0, 1, "}"});
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CFullFuncBodyDriver::_get_angmom_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    if ((integral[0] > 1) || (integral[1] > 1) || (integral[2] > 1) || (integral[3] > 1))
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
        
        if ((integral[2] > 1) && (integral[0] != integral[2]) && (integral[1] != integral[2]))
        {
            for (const auto& label : angmom.get_factors(integral[2]))
            {
                vstr.push_back("const double " + label + ";");
            }
        }
        
        if ((integral[3] > 1) && (integral[0] != integral[3]) && (integral[1] != integral[3]) && (integral[2] != integral[3]))
        {
            for (const auto& label : angmom.get_factors(integral[3]))
            {
                vstr.push_back("const double " + label + ";");
            }
        }
    }
    
    return vstr;
}

std::vector<std::string>
T4CFullFuncBodyDriver::_get_gtos_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// intialize GTO pairs data on bra side");

    vstr.push_back("const auto bra_gpair_coords = bra_gto_pair_block.getCoordinates();");

    vstr.push_back("const auto bra_gpair_exps = bra_gto_pair_block.getExponents();");

    vstr.push_back("const auto bra_gpair_norms = bra_gto_pair_block.getNormalizationFactors();");
    
    vstr.push_back("const auto bra_gpair_ovls = bra_gto_pair_block.getOverlapFactors();");

    vstr.push_back("const auto bra_nppairs = bra_gto_pair_block.getNumberOfPrimitivePairs();");
    
    vstr.push_back("const auto bra_ncpairs = bra_gto_pair_block.getNumberOfContractedPairs();");
    
    vstr.push_back("// intialize GTO pairs data on ket side");

    vstr.push_back("const auto ket_gpair_coords = ket_gto_pair_block.getCoordinates();");

    vstr.push_back("const auto ket_gpair_exps = ket_gto_pair_block.getExponents();");

    vstr.push_back("const auto ket_gpair_norms = ket_gto_pair_block.getNormalizationFactors();");
    
    vstr.push_back("const auto ket_gpair_ovls = ket_gto_pair_block.getOverlapFactors();");

    vstr.push_back("const auto ket_nppairs = ket_gto_pair_block.getNumberOfPrimitivePairs();");
    
    vstr.push_back("const auto ket_ncpairs = ket_gto_pair_block.getNumberOfContractedPairs();");
    
    vstr.push_back("// set up orbital indexes on bra and ket sides");
    
    vstr.push_back("const auto bra_orb_indexes = bra_gto_pair_block.getOrbitalIndexes();");
    
    vstr.push_back("const auto ket_orb_indexes = ket_gto_pair_block.getOrbitalIndexes();");
    
    vstr.push_back("// angular momentum on bra and ket sides");
    
    vstr.push_back("const auto bra_angmom = bra_gto_pair_block.getAngularMomentums();");
    
    vstr.push_back("const auto ket_angmom = ket_gto_pair_block.getAngularMomentums();");
    
    vstr.push_back("// local OMP storage for critical region");
    
    const auto label = std::to_string(integral.components<T2CPair, T2CPair>().size());
    
    vstr.push_back("auto tint_buffers = std::vector<std::vector<double>> (" + label + ", std::vector<double>(ket_ncpairs, 0.0));");
    
    return vstr;
}

std::vector<std::string>
T4CFullFuncBodyDriver::_get_comp_gtos_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// intialize GTO pairs data on bra side");

    vstr.push_back("const auto bra_gpair_coords = bra_gto_pair_block.getCoordinates();");

    vstr.push_back("const auto bra_gpair_exps = bra_gto_pair_block.getExponents();");

    vstr.push_back("const auto bra_gpair_norms = bra_gto_pair_block.getNormalizationFactors();");
    
    vstr.push_back("const auto bra_gpair_ovls = bra_gto_pair_block.getOverlapFactors();");

    vstr.push_back("const auto bra_nppairs = bra_gto_pair_block.getNumberOfPrimitivePairs();");
    
    vstr.push_back("const auto bra_ncpairs = bra_gto_pair_block.getNumberOfContractedPairs();");
    
    vstr.push_back("// intialize GTO pairs data on ket side");

    vstr.push_back("const auto ket_gpair_coords = ket_gto_pair_block.getCoordinates();");

    vstr.push_back("const auto ket_gpair_exps = ket_gto_pair_block.getExponents();");

    vstr.push_back("const auto ket_gpair_norms = ket_gto_pair_block.getNormalizationFactors();");
    
    vstr.push_back("const auto ket_gpair_ovls = ket_gto_pair_block.getOverlapFactors();");

    vstr.push_back("const auto ket_nppairs = ket_gto_pair_block.getNumberOfPrimitivePairs();");
    
    vstr.push_back("const auto ket_ncpairs = ket_gto_pair_block.getNumberOfContractedPairs();");
    
    vstr.push_back("// set up orbital indexes on bra and ket sides");
    
    vstr.push_back("const auto bra_orb_indexes = bra_gto_pair_block.getOrbitalIndexes();");
    
    vstr.push_back("const auto ket_orb_indexes = ket_gto_pair_block.getOrbitalIndexes();");
    
    vstr.push_back("// angular momentum on bra and ket sides");
    
    vstr.push_back("const auto bra_angmom = bra_gto_pair_block.getAngularMomentums();");
    
    vstr.push_back("const auto ket_angmom = ket_gto_pair_block.getAngularMomentums();");
    
    vstr.push_back("// local OMP storage for critical region");
    
    const auto label = std::to_string(integral.components<T2CPair, T2CPair>().size());
    
    vstr.push_back("TDoubleArray2D<" + label + "> tints;");
    
    if ((integral[0] + integral[2]) > 0)
    {
        vstr.push_back("// horizontal recursion storage");
        
        for (auto i = integral[1]; i <= (integral[0] + integral[1]); i++)
        {
            for (auto j = integral[2]; j <= (integral[2] + integral[3]); j++)
            {
                vstr.push_back("alignas(64) TDoubleArray& buffer_" + std::to_string(i) + std::to_string(j) + ";");
            }
        }
    }
    
    return vstr;
}

std::vector<std::string>
T4CFullFuncBodyDriver::_get_vars_def(const I4CIntegral& integral) const
{
    std::vector<std::string> vstr;

    vstr.push_back("// initialize aligned arrays for C and D centers");

    vstr.push_back("alignas(64) TDoubleArray coords_c_x;");

    vstr.push_back("alignas(64) TDoubleArray coords_c_y;");

    vstr.push_back("alignas(64) TDoubleArray coords_c_z;");

    vstr.push_back("alignas(64) TDoubleArray coords_d_x;");

    vstr.push_back("alignas(64) TDoubleArray coords_d_y;");

    vstr.push_back("alignas(64) TDoubleArray coords_d_z;");

    vstr.push_back("// initialize aligned arrays for ket side");

    vstr.push_back("alignas(64) TDoubleArray ket_exps_c;");

    vstr.push_back("alignas(64) TDoubleArray ket_exps_d;");

    vstr.push_back("alignas(64) TDoubleArray ket_norms;");
    
    vstr.push_back("alignas(64) TDoubleArray ket_ovls;");

    vstr.push_back("// initialize contracted integrals buffer");

    vstr.push_back("alignas(64) TDoubleArray buffer;");
    
    return vstr;
}

std::vector<std::string>
T4CFullFuncBodyDriver::_get_batches_def() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// loop over integral batches");
    
    vstr.push_back("const auto nbatches = batch::getNumberOfBatches(ket_ncpairs, simd_width);");
       
    return vstr;
}

void
T4CFullFuncBodyDriver::_add_batches_loop_start(VCodeLines& lines) const
{
    lines.push_back({1, 0, 1, "for (int64_t i = bra_first; i < bra_last; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "const auto [bra_coords_a, bra_coords_b]  = bra_gpair_coords[i];"});
    
    lines.push_back({2, 0, 1, "for (int64_t j = 0; j < nbatches; j++)"});
    
    lines.push_back({2, 0, 1, "{"});
}

void
T4CFullFuncBodyDriver::_add_comp_batches_loop_start(VCodeLines& lines) const
{
    lines.push_back({1, 0, 1, "for (int64_t i = bra_first; i < bra_last; i++)"});
    
    lines.push_back({1, 0, 1, "{"});
    
    lines.push_back({2, 0, 2, "// skip repeating integral buffers in diagonal blocks"});
    
    lines.push_back({2, 0, 1, "if (diagonal)"});
    
    lines.push_back({2, 0, 1, "{"});
    
    lines.push_back({3, 0, 1, "if (ket_last < i) continue;"});
    
    lines.push_back({2, 0, 2, "}"});
    
    lines.push_back({2, 0, 2, "const auto ket_dim = ket_last - ket_first;"});
    
    lines.push_back({2, 0, 2, "// load coordinates data on bra and ket sides"});
    
    lines.push_back({2, 0, 2, "const auto [bra_coords_a, bra_coords_b]  = bra_gpair_coords[i];"});
    
    lines.push_back({2, 0, 2, "simd::loadCoordinates(coords_c_x, coords_c_y, coords_c_z, coords_d_x, coords_d_y, coords_d_z, ket_gpair_coords, ket_first, ket_last);"});
}

void
T4CFullFuncBodyDriver::_add_batches_loop_body(      VCodeLines&  lines,
                                              const I4CIntegral& integral) const
{
    lines.push_back({3, 0, 2, "const auto [ket_first, ket_last] = batch::getBatchRange(j, ket_ncpairs, simd_width);"});
    
    lines.push_back({3, 0, 2, "// load coordinates data on ket side"});
    
    lines.push_back({3, 0, 2, "const auto ket_dim = ket_last - ket_first;"});
    
    lines.push_back({3, 0, 2, "// skip repeating integral buffers in diagonal blocks"});
    
    lines.push_back({3, 0, 1, "if (diagonal)"});
    
    lines.push_back({3, 0, 1, "{"});
    
    lines.push_back({4, 0, 1, "if (ket_last < i) continue;"});
    
    lines.push_back({3, 0, 2, "}"});
    
    lines.push_back({3, 0, 2, "simd::loadCoordinates(coords_c_x, coords_c_y, coords_c_z, coords_d_x, coords_d_y, coords_d_z, ket_gpair_coords, ket_first, ket_last);"});
    
    size_t index = 0;
    
    for (const auto& tcomp : integral.components<T2CPair, T2CPair>())
    {
        _add_component_body(lines, integral, tcomp, index);
    
        index++;
    }
  
//    lines.push_back({2, 0, 1, "#pragma omp critical"});
//
//    lines.push_back({2, 0, 1, "{"});
//
//    index = 0;
//
//    for (const auto& tcomp : integral.components<T2CPair, T2CPair>())
//    {
//        _write_block_distributor(lines, integral, tcomp, index);
//
//       //_add_component_body(lines, integral, tcomp, index);
//
//        index++;
//    }
//
//    //_write_block_distributor(lines, integral, component);
    
    //lines.push_back({2, 0, 1, "}"});
    
    //lines.push_back({1, 0, 1, "}"});
}

void
T4CFullFuncBodyDriver::_add_comp_batches_loop_body(      VCodeLines&  lines,
                                                   const I4CIntegral& integral) const
{
    size_t index = 0;
    
    for (const auto& tcomp : integral.components<T2CPair, T2CPair>())
    {
        _add_comp_component_body(lines, integral, tcomp, index);
    
        index++;
    }
    
    lines.push_back({2, 0, 1, "#pragma omp critical"});
    
    lines.push_back({2, 0, 1, "{"});
    
    index = 0;
    
    for (const auto& tcomp : integral.components<T2CPair, T2CPair>())
    {
        _write_comp_block_distributor(lines, integral, tcomp, index);
    
        index++;
    }
}

void
T4CFullFuncBodyDriver::_add_batches_loop_end(      VCodeLines&  lines,
                                             const I4CIntegral& integral) const
{
    lines.push_back({2, 0, 2, "}"});
    
    lines.push_back({2, 0, 2, "// distribute integrals into Fock matrix"});
    
    lines.push_back({2, 0, 1, "#pragma omp critical"});
    
    lines.push_back({2, 0, 1, "{"});
    
    int64_t index = 0;
    
    for (const auto& tcomp : integral.components<T2CPair, T2CPair>())
    {
        _write_block_distributor(lines, integral, tcomp, index);
    
        index++;
    }
    
    lines.push_back({2, 0, 1, "}"});
    
    lines.push_back({1, 0, 1, "}"});
}

void
T4CFullFuncBodyDriver::_add_component_body(      VCodeLines&  lines,
                                           const I4CIntegral& integral,
                                           const T4CIntegral& component,
                                           const size_t       index) const
{
    auto [nsize, name] = t4c::prim_full_compute_func_name(component, integral);
    
    name = t4c::namespace_label(integral) + "::" + name;
    
    lines.push_back({3, 0, 2, "// compute primitive integrals block (" + fstr::upcase(component.label()) + ")"});
    
    lines.push_back({3, 0, 2, "simd::zero(buffer);"});
    
    lines.push_back({3, 0, 1, "for (int64_t k = 0; k < ket_nppairs; k++)"});
    
    lines.push_back({3, 0, 1, "{"});
        
    lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsData(ket_norms, ket_gpair_norms, k, ket_ncpairs, ket_first, ket_last);"});
    
    lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsData(ket_ovls, ket_gpair_ovls, k, ket_ncpairs, ket_first, ket_last);"});
        
    lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsPairsData(ket_exps_c, ket_exps_d, ket_gpair_exps, k, ket_ncpairs, ket_first, ket_last);"});
    
    lines.push_back({4, 0, 1, "for (int64_t l = 0; l < bra_nppairs; l++)"});
    
    lines.push_back({4, 0, 1, "{"});
    
    lines.push_back({5, 0, 2, "const auto bra_index = l * bra_ncpairs + i;"});
    
    lines.push_back({5, 0, 2, "const auto [bra_exp_a, bra_exp_b] = bra_gpair_exps[bra_index];"});
    
    lines.push_back({5, 0, 2, "const auto bra_norm = bra_gpair_norms[bra_index];"});
    
    lines.push_back({5, 0, 2, "const auto bra_ovl = bra_gpair_ovls[bra_index];"});
    
    lines.push_back({5, 0, 1, name + "(buffer, use_rs, omega, bra_coords_a, bra_coords_b, coords_c_x, coords_c_y, coords_c_z, coords_d_x, coords_d_y, coords_d_z, bra_exp_a, bra_exp_b, bra_norm, bra_ovl, ket_exps_c, ket_exps_d, ket_norms, ket_ovls, ket_dim);"});
    
    lines.push_back({4, 0, 1, "}"});
    
    lines.push_back({3, 0, 2, "}"});
    
    lines.push_back({3, 0, 2, "// copy buffer for OMP critical"});
        
    lines.push_back({3, 0, 2, "t4cfunc::distribute(tint_buffers[" + std::to_string(index) + "], buffer, ket_first, ket_last);"});
    
    //lines.push_back({3, 0, 1, "#pragma omp critical"});
    
    //lines.push_back({3, 0, 1, "{"});
    
    //_write_block_distributor(lines, integral, component);
    
    //lines.push_back({3, 0, 2, "}"});
}

void
T4CFullFuncBodyDriver::_add_comp_component_body(      VCodeLines&  lines,
                                                const I4CIntegral& integral,
                                                const T4CIntegral& component,
                                                const size_t       index) const
{
    if ((integral[0] + integral[2]) > 0)
    {
        T4CHrrElectronRepulsionDriver t4c_eri_drv;
        
        const auto rdist = (t4c_eri_drv.create_recursion({component, }))[0];
        
        for (const auto& tcomp : rdist.unique_integrals())
        {
            auto tint = I4CIntegral(tcomp);
                
            auto [nsize, name] = t4c::prim_vrr_compute_func_name(tcomp, tint);
            
            name = t4c::namespace_label(tint) + "::" + name;
            
            lines.push_back({2, 0, 2, "// compute auxilary integrals block (" + fstr::upcase(tcomp.label()) + ")"});
            
            lines.push_back({2, 0, 2, "simd::zero(buffer_" + std::to_string(tint[1]) + std::to_string(tint[3]) + ");"});
            
            lines.push_back({2, 0, 1, "for (int64_t k = 0; k < ket_nppairs; k++)"});
            
            lines.push_back({2, 0, 1, "{"});
                
            lines.push_back({3, 0, 2, "simd::loadPrimitiveGTOsData(ket_norms, ket_gpair_norms, k, ket_ncpairs, ket_first, ket_last);"});
            
            lines.push_back({3, 0, 2, "simd::loadPrimitiveGTOsData(ket_ovls, ket_gpair_ovls, k, ket_ncpairs, ket_first, ket_last);"});
                
            lines.push_back({3, 0, 2, "simd::loadPrimitiveGTOsPairsData(ket_exps_c, ket_exps_d, ket_gpair_exps, k, ket_ncpairs, ket_first, ket_last);"});
            
            lines.push_back({3, 0, 1, "for (int64_t l = 0; l < bra_nppairs; l++)"});
            
            lines.push_back({3, 0, 1, "{"});
            
            lines.push_back({4, 0, 2, "const auto bra_index = l * bra_ncpairs + i;"});
            
            lines.push_back({4, 0, 2, "const auto [bra_exp_a, bra_exp_b] = bra_gpair_exps[bra_index];"});
            
            lines.push_back({4, 0, 2, "const auto bra_norm = bra_gpair_norms[bra_index];"});
            
            lines.push_back({4, 0, 2, "const auto bra_ovl = bra_gpair_ovls[bra_index];"});
            
            lines.push_back({4, 0, 1, name + "(buffer_" + std::to_string(tint[1]) + std::to_string(tint[3]) + ", use_rs, omega, bra_coords_a, bra_coords_b, coords_c_x, coords_c_y, coords_c_z, coords_d_x, coords_d_y, coords_d_z, bra_exp_a, bra_exp_b, bra_norm, bra_ovl, ket_exps_c, ket_exps_d, ket_norms, ket_ovls, ket_dim);"});
            
            lines.push_back({3, 0, 1, "}"});
            
            lines.push_back({2, 0, 2, "}"});
        }
        
        lines.push_back({2, 0, 2, "// compute contracted integrals block (" + fstr::upcase(component.label()) + ")"});
        
        auto [nsize, name] = t4c::contr_hrr_compute_func_name(component, integral);
        
        name = t4c::namespace_label(integral) + "::" + name;
    }
    else
    {
        auto [nsize, name] = t4c::prim_vrr_compute_func_name(component, integral);
        
        name = t4c::namespace_label(integral) + "::" + name;
        
        lines.push_back({2, 0, 2, "// compute primitive integrals block (" + fstr::upcase(component.label()) + ")"});
        
        lines.push_back({2, 0, 2, "simd::zero(buffer);"});
        
        lines.push_back({2, 0, 1, "for (int64_t k = 0; k < ket_nppairs; k++)"});
        
        lines.push_back({2, 0, 1, "{"});
            
        lines.push_back({3, 0, 2, "simd::loadPrimitiveGTOsData(ket_norms, ket_gpair_norms, k, ket_ncpairs, ket_first, ket_last);"});
        
        lines.push_back({3, 0, 2, "simd::loadPrimitiveGTOsData(ket_ovls, ket_gpair_ovls, k, ket_ncpairs, ket_first, ket_last);"});
            
        lines.push_back({3, 0, 2, "simd::loadPrimitiveGTOsPairsData(ket_exps_c, ket_exps_d, ket_gpair_exps, k, ket_ncpairs, ket_first, ket_last);"});
        
        lines.push_back({3, 0, 1, "for (int64_t l = 0; l < bra_nppairs; l++)"});
        
        lines.push_back({3, 0, 1, "{"});
        
        lines.push_back({4, 0, 2, "const auto bra_index = l * bra_ncpairs + i;"});
        
        lines.push_back({4, 0, 2, "const auto [bra_exp_a, bra_exp_b] = bra_gpair_exps[bra_index];"});
        
        lines.push_back({4, 0, 2, "const auto bra_norm = bra_gpair_norms[bra_index];"});
        
        lines.push_back({4, 0, 2, "const auto bra_ovl = bra_gpair_ovls[bra_index];"});
        
        lines.push_back({4, 0, 1, name + "(buffer, use_rs, omega, bra_coords_a, bra_coords_b, coords_c_x, coords_c_y, coords_c_z, coords_d_x, coords_d_y, coords_d_z, bra_exp_a, bra_exp_b, bra_norm, bra_ovl, ket_exps_c, ket_exps_d, ket_norms, ket_ovls, ket_dim);"});
        
        lines.push_back({3, 0, 1, "}"});
        
        lines.push_back({2, 0, 2, "}"});
    }
    
    lines.push_back({2, 0, 2, "tints[" +  std::to_string(index) + "] = buffer;"});
    
    // lines.push_back({3, 0, 2, "// copy buffer for OMP critical"});
        
    //lines.push_back({3, 0, 2, "t4cfunc::distribute(tint_buffers[" + std::to_string(index) + "], buffer, ket_first, ket_last);"});
    
    //lines.push_back({3, 0, 1, "#pragma omp critical"});
    
    //lines.push_back({3, 0, 1, "{"});
    
    //_write_block_distributor(lines, integral, component);
    
    //lines.push_back({3, 0, 2, "}"});
}

void
T4CFullFuncBodyDriver::_write_block_distributor(      VCodeLines&  lines,
                                                const I4CIntegral& integral,
                                                const T4CIntegral& component,
                                                const size_t       index) const
{
    // create angular momentum data
    
    const auto amom = SphericalMomentum(integral[0]);
    
    const auto bmom = SphericalMomentum(integral[1]);
    
    const auto cmom = SphericalMomentum(integral[2]);
    
    const auto dmom = SphericalMomentum(integral[3]);
    
    // set up tensor component indexes
    
    auto idxa = t2c::tensor_component_index(component[0]);
    
    auto idxb = t2c::tensor_component_index(component[1]);
    
    auto idxc = t2c::tensor_component_index(component[2]);
    
    auto idxd = t2c::tensor_component_index(component[3]);
    
    // select angular factor pairs
    
    const auto apairs = amom.select_pairs(idxa);
    
    const auto bpairs = bmom.select_pairs(idxb);
    
    const auto cpairs = cmom.select_pairs(idxc);
    
    const auto dpairs = dmom.select_pairs(idxd);
    
    for (const auto& apair : apairs)
    {
        for (const auto& bpair : bpairs)
        {
            for (const auto& cpair : cpairs)
            {
                for (const auto& dpair : dpairs)
                {
                    auto lfactor = t2c::combine_factors(apair.second, bpair.second);
                    
                    lfactor = t2c::combine_factors(lfactor, cpair.second);
                    
                    lfactor = t2c::combine_factors(lfactor, dpair.second);
                    
                    lfactor = (lfactor == "1.0")  ? "" : lfactor + ", ";
                    
                    auto label = "{"  + std::to_string(apair.first) +
                                 ", " + std::to_string(bpair.first) +
                                 ", " + std::to_string(cpair.first) +
                                 ", " + std::to_string(dpair.first) + "}";
                    
                    lines.push_back({3, 0, 2, "t4cfunc::distribute(fock_matrix, density, tint_buffers[" + std::to_string(index) + "], bra_orb_indexes, ket_orb_indexes, bra_angmom, ket_angmom, " + lfactor + label + ", diagonal, i, 0, ket_ncpairs);"});
                }
            }
        }
    }
}

void
T4CFullFuncBodyDriver::_write_comp_block_distributor(      VCodeLines&  lines,
                                                     const I4CIntegral& integral,
                                                     const T4CIntegral& component,
                                                     const size_t       index) const
{
    // create angular momentum data
    
    const auto amom = SphericalMomentum(integral[0]);
    
    const auto bmom = SphericalMomentum(integral[1]);
    
    const auto cmom = SphericalMomentum(integral[2]);
    
    const auto dmom = SphericalMomentum(integral[3]);
    
    // set up tensor component indexes
    
    auto idxa = t2c::tensor_component_index(component[0]);
    
    auto idxb = t2c::tensor_component_index(component[1]);
    
    auto idxc = t2c::tensor_component_index(component[2]);
    
    auto idxd = t2c::tensor_component_index(component[3]);
    
    // select angular factor pairs
    
    const auto apairs = amom.select_pairs(idxa);
    
    const auto bpairs = bmom.select_pairs(idxb);
    
    const auto cpairs = cmom.select_pairs(idxc);
    
    const auto dpairs = dmom.select_pairs(idxd);
    
    for (const auto& apair : apairs)
    {
        for (const auto& bpair : bpairs)
        {
            for (const auto& cpair : cpairs)
            {
                for (const auto& dpair : dpairs)
                {
                    auto lfactor = t2c::combine_factors(apair.second, bpair.second);
                    
                    lfactor = t2c::combine_factors(lfactor, cpair.second);
                    
                    lfactor = t2c::combine_factors(lfactor, dpair.second);
                    
                    lfactor = (lfactor == "1.0")  ? "" : lfactor + ", ";
                    
                    auto label = "{"  + std::to_string(apair.first) +
                                 ", " + std::to_string(bpair.first) +
                                 ", " + std::to_string(cpair.first) +
                                 ", " + std::to_string(dpair.first) + "}";
                    
                    lines.push_back({3, 0, 2, "t4cfunc::distribute(fock_matrix, density, tint_buffers[" + std::to_string(index) + "], bra_orb_indexes, ket_orb_indexes, bra_angmom, ket_angmom, " + lfactor + label + ", diagonal, i, ket_first, ket_last);"});
                }
            }
        }
    }
}
