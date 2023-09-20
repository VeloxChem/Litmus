#include "t4c_full_docs.hpp"

#include "file_stream.hpp"
#include "t4c_utils.hpp"
#include "string_formater.hpp"

void
T4CFullDocuDriver::write_doc_str(      std::ofstream& fstream,
                                 const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
        
    lines.push_back({0, 0, 1, "/**"});
        
    lines.push_back({0, 0, 2, _get_compute_str(integral)});
    
    for (const auto& label : _get_vars_str())
    {
        lines.push_back({0, 1, 1, label});
    }
    
    lines.push_back({0, 0, 1, "*/"});
        
    ost::write_code_lines(fstream, lines);
}

std::string
T4CFullDocuDriver::_get_compute_str(const I4CIntegral& integral) const
{
    const auto bra_a = Tensor(integral[0]);
    
    const auto bra_b = Tensor(integral[1]);
        
    const auto ket_a = Tensor(integral[2]);
    
    const auto ket_b = Tensor(integral[3]);
        
    const auto integrand = integral.integrand();
            
    auto label = " Evaluates <"  + bra_a.label() + bra_b.label() + "|";
        
    label += t4c::integrand_label(integral.integrand());
        
    label += "|" + ket_a.label()  + ket_b.label()+ ">  integrals for given ";
        
    label += "GTOs pair blocks.";
    
    return label;
    
}

std::vector<std::string>
T4CFullDocuDriver::_get_vars_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("@param fock_matrices the pointer to Fock matrices.");
    
    vstr.push_back("@param bra_gto_pair_block the GTOs pair block for bra side.");
    
    vstr.push_back("@param ket_gto_pair_block the GTOs pair block for ket side.");
    
    vstr.push_back("@param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.");
    
    vstr.push_back("@param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.");
    
    return vstr;
}
