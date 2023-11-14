#include "t4c_full_decl.hpp"

#include "file_stream.hpp"
#include "t4c_utils.hpp"

void
T4CFullDeclDriver::write_func_decl(      std::ofstream& fstream,
                               const I4CIntegral&   integral,
                               const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_vars_str(integral, terminus))
    {
        if  (label.find(";") == std::string::npos)
        {
            lines.push_back({0, 0, 1, label});
        }
        else
        {
            lines.push_back({0, 0, 2, label});
        }
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T4CFullDeclDriver::write_prim_func_decl(      std::ofstream& fstream,
                                        const T4CIntegral&   component,
                                        const I4CIntegral&   integral,
                                        const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_prim_vars_str(component, integral, terminus))
    {
        if  (label.find(";") == std::string::npos)
        {
            lines.push_back({0, 0, 1, label});
        }
        else
        {
            lines.push_back({0, 0, 2, label});
        }
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T4CFullDeclDriver::write_vrr_func_decl(      std::ofstream& fstream,
                                       const T4CIntegral&   component,
                                       const I4CIntegral&   integral,
                                       const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_vrr_vars_str(component, integral, terminus))
    {
        if  (label.find(";") == std::string::npos)
        {
            lines.push_back({0, 0, 1, label});
        }
        else
        {
            lines.push_back({0, 0, 2, label});
        }
    }
    
    ost::write_code_lines(fstream, lines);
}

void
T4CFullDeclDriver::write_hrr_func_decl(      std::ofstream& fstream,
                                       const T4CIntegral&   component,
                                       const I4CIntegral&   integral,
                                       const bool           terminus) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "auto"});
    
    for (const auto& label : _get_hrr_vars_str(component, integral, terminus))
    {
        if  (label.find(";") == std::string::npos)
        {
            lines.push_back({0, 0, 1, label});
        }
        else
        {
            lines.push_back({0, 0, 2, label});
        }
    }
    
    ost::write_code_lines(fstream, lines);
}

std::vector<std::string>
T4CFullDeclDriver::_get_vars_str(const I4CIntegral& integral,
                                 const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto [nsize, name] = t4c::full_compute_func_name(integral);
    
    vstr.push_back(name + "(CFockMatrix* fock_matrix,");
    
    vstr.push_back(std::string(nsize, ' ') + "const CMatrix* density,");
    
    vstr.push_back(std::string(nsize, ' ') + "const CGtoPairBlock& bra_gto_pair_block,");
    
    vstr.push_back(std::string(nsize, ' ') + "const CGtoPairBlock& ket_gto_pair_block,");
    
    vstr.push_back(std::string(nsize, ' ') + "const bool diagonal,");
    
    vstr.push_back(std::string(nsize, ' ') + "const bool use_rs,");
    
    vstr.push_back(std::string(nsize, ' ') + "const double omega,");
    
    vstr.push_back(std::string(nsize, ' ') + "const int64_t bra_first,");
    
    //vstr.push_back(std::string(nsize, ' ') + "const int64_t bra_last,");
    
    //vstr.push_back(std::string(nsize, ' ') + "const int64_t ket_first,");
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    //vstr.push_back(std::string(nsize, ' ') + "const int64_t ket_last) -> void" +  tsymbol);
    
    vstr.push_back(std::string(nsize, ' ') + "const int64_t bra_last) -> void" +  tsymbol);
    
    return vstr;
}

std::vector<std::string>
T4CFullDeclDriver::_get_prim_vars_str(const T4CIntegral& component,
                                      const I4CIntegral& integral,
                                      const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto [nsize, name] = t4c::prim_full_compute_func_name(component, integral);
    
    vstr.push_back(name + "(TDoubleArray& buffer,");
    
    vstr.push_back(std::string(nsize, ' ') + "const bool use_rs,");
    
    vstr.push_back(std::string(nsize, ' ') + "const double omega,");
   
    vstr.push_back(std::string(nsize, ' ') + "const TPoint3D& coords_a,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TPoint3D& coords_b,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_c_x,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_c_y,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_c_z,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_d_x,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_d_y,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_d_z,");
    
    vstr.push_back(std::string(nsize, ' ') + "const double bra_exp_a,");
    
    vstr.push_back(std::string(nsize, ' ') + "const double bra_exp_b,");
    
    vstr.push_back(std::string(nsize, ' ') + "const double bra_norm,");
    
    vstr.push_back(std::string(nsize, ' ') + "const double bra_ovl,");
 
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& ket_exps_c,");
        
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& ket_exps_d,");
        
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& ket_norms,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& ket_ovls,");
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    vstr.push_back(std::string(nsize, ' ') + "const int64_t       ket_dim) -> void" + tsymbol);
    
    return vstr;
}

std::vector<std::string>
T4CFullDeclDriver::_get_vrr_vars_str(const T4CIntegral& component,
                                     const I4CIntegral& integral,
                                     const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto [nsize, name] = t4c::prim_vrr_compute_func_name(component, integral);
    
    vstr.push_back(name + "(TDoubleArray& buffer,");
    
    vstr.push_back(std::string(nsize, ' ') + "const bool use_rs,");
    
    vstr.push_back(std::string(nsize, ' ') + "const double omega,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TPoint3D& coords_a,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TPoint3D& coords_b,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_c_x,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_c_y,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_c_z,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_d_x,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_d_y,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_d_z,");
    
    vstr.push_back(std::string(nsize, ' ') + "const double bra_exp_a,");
    
    vstr.push_back(std::string(nsize, ' ') + "const double bra_exp_b,");
    
    vstr.push_back(std::string(nsize, ' ') + "const double bra_norm,");
    
    vstr.push_back(std::string(nsize, ' ') + "const double bra_ovl,");
 
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& ket_exps_c,");
        
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& ket_exps_d,");
        
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& ket_norms,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& ket_ovls,");
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    vstr.push_back(std::string(nsize, ' ') + "const int64_t       ket_dim) -> void" + tsymbol);
    
    return vstr;
}

std::vector<std::string>
T4CFullDeclDriver::_get_hrr_vars_str(const T4CIntegral& component,
                                     const I4CIntegral& integral,
                                     const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto [nsize, name] = t4c::contr_hrr_compute_func_name(component, integral);
    
    vstr.push_back(name + "(TDoubleArray& buffer,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TPoint3D& coords_a,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TPoint3D& coords_b,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_c_x,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_c_y,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_c_z,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_d_x,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_d_y,");
    
    vstr.push_back(std::string(nsize, ' ') + "const TDoubleArray& coords_d_z,");
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    vstr.push_back(std::string(nsize, ' ') + "const int64_t       ket_dim) -> void" + tsymbol);
    
    return vstr;
}
