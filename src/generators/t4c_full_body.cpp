#include "t4c_full_body.hpp"

#include "string_formater.hpp"
#include "spherical_momentum.hpp"
#include "t4c_utils.hpp"

void
T4CFullFuncBodyDriver::write_func_body(      std::ofstream& fstream,
                                       const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    for (const auto& label : _get_gtos_def())
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
    
    _add_batches_loop_end(lines);
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CFullFuncBodyDriver::_get_gtos_def() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("// intialize GTO pairs data on bra side");

    vstr.push_back("const auto bra_gpair_coords = bra_gto_pair_block.getCoordinates();");

    vstr.push_back("const auto bra_gpair_exps = bra_gto_pair_block.getExponents();");

    vstr.push_back("const auto bra_gpair_norms = bra_gto_pair_block.getNormalizationFactors();");

    vstr.push_back("const auto bra_nppairs = bra_gto_pair_block.getNumberOfPrimitivePairs();");
    
    vstr.push_back("const auto bra_ncpairs = bra_gto_pair_block.getNumberOfContractedPairs();");
    
    vstr.push_back("// intialize GTO pairs data on ket side");

    vstr.push_back("const auto ket_gpair_coords = ket_gto_pair_block.getCoordinates();");

    vstr.push_back("const auto ket_gpair_exps = ket_gto_pair_block.getExponents();");

    vstr.push_back("const auto ket_gpair_norms = ket_gto_pair_block.getNormalizationFactors();");

    vstr.push_back("const auto ket_nppairs = ket_gto_pair_block.getNumberOfPrimitivePairs();");
    
    vstr.push_back("const auto ket_ncpairs = ket_gto_pair_block.getNumberOfContractedPairs();");
    
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
    lines.push_back({1, 0, 1, "for (int64_t i = 0; i < nbatches; i++)"});
        
    lines.push_back({1, 0, 1, "{"});
}

void
T4CFullFuncBodyDriver::_add_batches_loop_body(      VCodeLines&  lines,
                                              const I4CIntegral& integral) const
{
    lines.push_back({2, 0, 2, "const auto [ket_first, ket_last] = batch::getBatchRange(i, ncpairs, simd_width);"});
    
    lines.push_back({2, 0, 2, "const auto ket_dim = last - first;"});

    lines.push_back({2, 0, 2, "// load coordinates data on ket side"});
    
    lines.push_back({2, 0, 2, "simd::loadCoordinates(coords_c_x, coords_c_y, coords_c_z, coords_d_x, coords_d_y, coords_d_z, key_gpair_coords, ket_first, ket_last);"});
    
    lines.push_back({2, 0, 1, "for (int64_t j = bra_first; j < bra_last; j++)"});
    
    lines.push_back({2, 0, 1, "{"});
    
    lines.push_back({3, 0, 2, "const auto [bra_coords_a, bra_coords_b]  = bra_gpair_coords[j];"});
        
    for (const auto& tcomp : integral.components<T2CPair, T2CPair>())
    {
       _add_component_body(lines, integral, tcomp);
    }
    
    lines.push_back({2, 0, 1, "}"});
}

void
T4CFullFuncBodyDriver::_add_batches_loop_end(VCodeLines& lines) const
{
    lines.push_back({1, 0, 1, "}"});
}

void
T4CFullFuncBodyDriver::_add_component_body(      VCodeLines&  lines,
                                           const I4CIntegral& integral,
                                           const T4CIntegral& component) const
{
    auto [nsize, name] = t4c::prim_full_compute_func_name(component, integral);
    
    name = t4c::namespace_label(integral) + "::" + name;
    
    lines.push_back({3, 0, 2, "// compute primitive integrals block (" + fstr::upcase(component.label()) + ")"});
    
    lines.push_back({3, 0, 2, "simd::zero(buffer);"});
    
    lines.push_back({3, 0, 1, "for (int64_t k = 0; k < ket_nppairs; k++)"});
    
    lines.push_back({3, 0, 1, "{"});
        
    lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsData(ket_norms, ket_gpair_norms, k, ket_ncpairs, ket_first, ket_last);"});
        
    lines.push_back({4, 0, 2, "simd::loadPrimitiveGTOsPairsData(ket_exps_c, ket_exps_d, ket_gpair_exps, k, ket_ncpairs, ket_first, ket_last);"});
    
    lines.push_back({4, 0, 1, "for (int64_t l = 0; l < bra_nppairs; l++)"});
    
    lines.push_back({4, 0, 1, "{"});
    
    lines.push_back({5, 0, 2, "const auto bra_index = l * bra_ncpairs + j;"});
    
    lines.push_back({5, 0, 2, "const auto [bra_exp_a, bra_exp_b] = bra_gpair_exps[bra_index];"});
    
    lines.push_back({5, 0, 2, "const auto bra_norm = bra_gpair_norms[bra_index];"});
    
    lines.push_back({5, 0, 1, name + "(buffer, bra_coords_a, bra_coords_b, coords_c_x, coords_c_y, coords_c_z, coords_d_x, coords_d_y, coords_d_z, bra_exps_a, bra_exps_b, bra_norm, ket_exps_c, ket_exps_d, ket_norms, ket_ndim);"});
    
    lines.push_back({4, 0, 1, "}"});
    
    lines.push_back({3, 0, 2, "}"});
    
    lines.push_back({3, 0, 2, "t4c::distribute(fock_matrices, buffer, bra_gpair_block, ket_gpair_block, j, ket_first, ket_last);"});
}
