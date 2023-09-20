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

std::vector<std::string>
T4CFullDeclDriver::_get_vars_str(const I4CIntegral& integral,
                                 const bool         terminus) const
{
    std::vector<std::string> vstr;
    
    const auto [nsize, name] = t4c::full_compute_func_name(integral);
    
    vstr.push_back(name + "(CFockMatrices* matrices,");
    
    vstr.push_back(std::string(nsize, ' ') + "const CGtoPairBlock& bra_gto_pair_block,");
    
    vstr.push_back(std::string(nsize, ' ') + "const CGtoPairBlock& ket_gto_pair_block,");
    
    vstr.push_back(std::string(nsize, ' ') + "const int64_t bra_first,");
    
    const auto tsymbol = (terminus) ? ";" : "";
    
    vstr.push_back(std::string(nsize, ' ') + "const int64_t bra_last) -> void" +  tsymbol);
    
    return vstr;
}
