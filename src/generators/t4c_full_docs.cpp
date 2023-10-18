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

void
T4CFullDocuDriver::write_prim_doc_str(      std::ofstream& fstream,
                                      const T4CIntegral&   component,
                                      const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "/**"});
    
    lines.push_back({0, 1, 2, _get_prim_compute_str(component, integral)});
    
    for (const auto& label : _get_prim_vars_str())
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

std::string
T4CFullDocuDriver::_get_prim_compute_str(const T4CIntegral& component,
                                         const I4CIntegral& integral) const
{
    const auto bra_a = Tensor(integral[0]);
    
    const auto bra_b = Tensor(integral[1]);
        
    const auto ket_a = Tensor(integral[2]);
    
    const auto ket_b = Tensor(integral[3]);

    auto label = "Evaluates block of primitive <" + bra_a.label() + bra_b.label() ;

    label += "|" + t4c::integrand_label(integral.integrand()) + "|";

    label += ket_a.label()  + ket_b.label();

    label += ">  (" +  fstr::upcase(component.label())  +  ") integrals.";

    return label;
}

std::vector<std::string>
T4CFullDocuDriver::_get_vars_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("@param fock_matrix the pointer to Fock matrix.");
    
    vstr.push_back("@param density the AO density matrix.");
    
    vstr.push_back("@param bra_gto_pair_block the GTOs pair block for bra side.");
    
    vstr.push_back("@param ket_gto_pair_block the GTOs pair block for ket side.");
                                            
    vstr.push_back("@param diagonal the flag signaling diagonal contributions.");
                                            
    vstr.push_back("@param use_rs the flag to use range separated form of electron repulsion integrals.");
                                            
    vstr.push_back("@param omega the range separation factor.");
    
    vstr.push_back("@param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.");
    
    vstr.push_back("@param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.");
    
    vstr.push_back("@param ket_first the index of the range [ket_first, ket_last) of GTOs on ket side.");
    
    vstr.push_back("@param ket_last the index of the range [ket_first, ket_last) of GTOs on ket side.");
    
    return vstr;
}

std::vector<std::string>
T4CFullDocuDriver::_get_prim_vars_str() const
{
    std::vector<std::string> vstr;
    
    vstr.push_back("@param buffer the integrals buffer.");
    
    vstr.push_back("@param use_rs the flag to use range separated form of electron repulsion integrals.");
    
    vstr.push_back("@param omega the range separation factor."); 
    
    vstr.push_back("@param coords_a the Cartesian coordinates of center A.");
    
    vstr.push_back("@param coords_b the Cartesian coordinates of center B.");
    
    vstr.push_back("@param coords_c_x the array of Cartesian X coordinates on center C.");
    
    vstr.push_back("@param coords_c_y the array of Cartesian Y coordinates on center C.");
    
    vstr.push_back("@param coords_c_z the array of Cartesian Z coordinates on center C.");
    
    vstr.push_back("@param coords_d_x the array of Cartesian X coordinates on center D.");
    
    vstr.push_back("@param coords_d_y the array of Cartesian Y coordinates on center D.");
    
    vstr.push_back("@param coords_d_z the array of Cartesian Z coordinates on center D.");
    
    vstr.push_back("@param bra_exp_a the exponent on bra center A.");
    
    vstr.push_back("@param bra_exp_b the exponent on bra center B.");
    
    vstr.push_back("@param bra_norm the normalization factor on bra side.");
    
    vstr.push_back("@param bra_ovl the overlap factor on bra side.");
    
    vstr.push_back("@param ket_exps_c the array of exponents on ket center C.");
        
    vstr.push_back("@param ket_exps_d the array of exponents on ket center D.");
        
    vstr.push_back("@param ket_norms the array of normalization factors on ket side.");
    
    vstr.push_back("@param ket_ovls the array of overlap factors on ket side.");
     
    vstr.push_back("@param ket_dim the size of integrals batch on ket side.");
        
    return vstr;
}
