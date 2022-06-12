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

#include "eri_cpu_generator.hpp"

#include "file_stream.hpp"

EriCPUGenerator::EriCPUGenerator()

    : _diag_form(false)
{
    
}

void
EriCPUGenerator::set_diag_form()
{
    _diag_form = true;
}

void
EriCPUGenerator::generate(const Repository<R4Group, T4CIntegral>& repo) const
{
    for (const auto& tint : repo.base<I4CIntegral>())
    {
        // generates VRR recursions
        
        if (_is_vrr_rec(tint))
        {
            _write_vrr_cpp_header(tint, repo);
        }
        
        // generates HRR recursions
        
        if (_is_hrr_rec(tint))
        {
            _write_hrr_cpp_header(tint, repo);
        }
    }
}

void
EriCPUGenerator::_write_vrr_cpp_header(const I4CIntegral&                      integral,
                                       const Repository<R4Group, T4CIntegral>& repo) const
{
    std::string fname = _file_name(integral, "VRR") + ".hpp";
        
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    ost::write_copyright(fstream);
    
    ost::write_vrr_includes(fstream);
    
    ost::write_namespace(fstream, "derirec", true);
    
    const auto tmaps = repo.base_map<I4CIntegral>(integral);
    
    for (const auto& tval : tmaps)
    {
        _write_hvrr_func_decl(fstream, tmaps, tval.first);
        
        _write_vrr_func_body(fstream, tval.first, tval.second);
    }
    
    ost::write_namespace(fstream, "derirec", false);
    
    fstream.close();
}

void
EriCPUGenerator::_write_hrr_cpp_header(const I4CIntegral&                      integral,
                                       const Repository<R4Group, T4CIntegral>& repo) const
{
    std::string fname = _file_name(integral, "HRR") + ".hpp";
        
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    ost::write_copyright(fstream);
    
    ost::write_vrr_includes(fstream);
    
    ost::write_namespace(fstream, "derirec", true);
    
    const auto tmaps = repo.base_map<I4CIntegral>(integral);
    
    for (const auto& tval : tmaps)
    {
        _write_hvrr_func_decl(fstream, tmaps, tval.first);
        
        _write_hrr_func_body(fstream, tval.first, tval.second);
    }
    
    ost::write_namespace(fstream, "derirec", false);
    
    fstream.close();
}


std::string
EriCPUGenerator::_file_name(const I4CIntegral& integral,
                            const std::string& rectype) const
{
    std::string fname;
    
    if (integral.integrand() == Operator("1/|r-r'|")) fname += "Eri";
    
    if (_diag_form) fname += "Diag";
    
    fname += rectype + "For" + integral.label(); 
    
    return fname;
}

std::string
EriCPUGenerator::_buffer_name(const I4CIntegral& integral,
                              const bool         flg_hrr) const
{
    std::string name = "intsBuffer"; 
    
    name += integral.label();
    
    if (!flg_hrr) name += std::to_string(integral.order());
    
    return name;
}

std::string
EriCPUGenerator::_indexes_name(const I4CIntegral& integral,
                               const bool         flg_hrr) const
{
    std::string name = "intsIndexes";
    
    name += integral.label();
    
    if (!flg_hrr) name += std::to_string(integral.order());
    
    return name;
}

std::string
EriCPUGenerator::_factor_name(const std::string& label) const
{
    // distances PB, QD, WP, WQ
    
    if ((label == "PB") || (label == "QD") ||
        (label == "WP") || (label == "WQ") ||
        (label == "AB") || (label == "CD"))
    {
        return "rDistances" + label;
    }
    
    // Obara-Saika factors
    
    if (label == "1/zeta")
    {
        return "osFactorsBraZeta";
    }
    
    if (label == "1/eta")
    {
        return "osFactorsKetZeta";
    }
    
    if (label == "1/(zeta+eta)")
    {
        return "osFactorsZeta";
    }
    
    if (label == "rho/zeta^2")
    {
        return "osFactorsBraRhoZeta";
    }
    
    if (label == "rho/eta^2")
    {
        return "osFactorsKetRhoZeta";
    }
    
    return std::string();
}

std::string
EriCPUGenerator::_fraction_name(const Fraction& fraction) const
{
    std::string label = "fact_";
    
    label += std::to_string(fraction.numerator());
    
    if (fraction.denominator() != 1)
    {
        label += "_" + std::to_string(fraction.denominator());
    }
    
    return label;
}

std::string
EriCPUGenerator::_rec_term_name(const R4CTerm&     recterm,
                                const std::string& index,
                                const bool         first,
                                const bool         flg_hrr) const
{
    std::string label;
    
    // sign of recursion term
    
    auto pfact = recterm.prefactor();
    
    if (pfact.is_negative())
    {
        label = (first) ? "-" : "- ";
            
        pfact = pfact.abs();
    }
    else
    {
        label = (first) ? "" : "+ ";
    }
    
    // prefactor of recursion term
    
    if (pfact != Fraction(1)) label += _fraction_name(pfact) + " * ";
    
    // factors of recursion term
        
    for (const auto& tval : recterm.factors())
    {
        label += tval.label() + index  + " * ";
    }
    
    // integral component of recursion term
    
    label += "t_" + recterm.label(!flg_hrr) + index;
    
    return label;
}


bool
EriCPUGenerator::_is_hrr_rec(const I4CIntegral& integral) const
{
    if ((integral[0] > 0) || (integral[2] > 0))
    {
        return true;
    }
    else
    {
        return false;
    }
}


bool
EriCPUGenerator::_is_vrr_rec(const I4CIntegral& integral) const
{
    if (((integral[0] + integral[2]) == 0) &&
        ((integral[1] + integral[3]) > 0))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool
EriCPUGenerator::_is_aux_rec(const I4CIntegral& integral) const
{
    if ((integral[0] + integral[1] + integral[2] + integral[3]) == 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void
EriCPUGenerator::_write_hvrr_func_decl(      std::ofstream&                             fstream,
                                       const std::map<Signature<T4CIntegral>, R4Group>& signatures,
                                       const Signature<T4CIntegral>&                    signature) const
{
    const std::vector<std::string> vlabels = {"BufferHostXY<T>&      ", "BufferHostX<int32_t>& ",
                                              "BufferHostMY<T, 3>&   ", "T*                    ",
                                              "int32_t               ", "bool                  "};
        
    // write function declaration
    
    fstream << "template <typename T>" << std::endl << "auto" << std::endl;
    
    const auto flabel = _hvrr_func_name(signatures, signature);
    
    fstream << flabel << "(" << std::string(6, ' ');
    
    // accumulation integral
    
    const auto rint = signature.base<I4CIntegral>();
    
    fstream << vlabels[0] << _buffer_name(*rint, true) << "," << std::endl;
    
    const auto space = std::string(flabel.size() + 1, ' ');
        
    fstream << space << "const " << vlabels[1] << _indexes_name(*rint, _is_hrr_rec(*rint)) << "," << std::endl;
    
    // recursion integrals
    
    for (const auto& tint : signature.expansion<I4CIntegral>())
    {
        if ((tint[0] + tint[1] + tint[2] + tint[3]) == 0)
        {
            fstream << space << "const " << vlabels[3] << _buffer_name(tint, _is_hrr_rec(*rint)) << "," << std::endl;
        }
        else
        {
            fstream << space << "const " << vlabels[0] << _buffer_name(tint, _is_hrr_rec(*rint)) << "," << std::endl;
                
            fstream << space << "const " << vlabels[1] << _indexes_name(tint, _is_hrr_rec(*rint)) << "," << std::endl;
        }
    }
    
    // recursion factors
    
    for (const auto& tfact : signature.factor_names())
    {
        if (const auto fname = _factor_name(tfact); fname.find("rDistances") == std::string::npos)
        {
            fstream << space << "const " << vlabels[3] << fname << "," << std::endl;
        }
        else
        {
            fstream << space <<"const " << vlabels[2] << fname << "," << std::endl;
        }
    }
    
    // other input parameters
    
    if (_is_vrr_rec(*rint))
    {
        fstream << space <<"const " << vlabels[5] << "useSummation," << std::endl;
    }
    
    fstream << space <<"const " << vlabels[4] << "nBatchPairs) -> void" << std::endl;
}

void
EriCPUGenerator::_write_vrr_func_body(      std::ofstream&          fstream,
                                      const Signature<T4CIntegral>& signature,
                                      const R4Group&                recgroup) const
{
    fstream << "{" << std::endl;
    
    _write_os_factors(fstream, signature);
    
    _write_distances(fstream, signature);
    
    _write_buffers(fstream, signature, false);
    
    _write_fractions(fstream, recgroup);
    
    _write_vrr_loop(fstream, recgroup);
    
    fstream << "}" << std::endl << std::endl;
}

void
EriCPUGenerator::_write_hrr_func_body(      std::ofstream&          fstream,
                                      const Signature<T4CIntegral>& signature,
                                      const R4Group&                recgroup) const
{
    fstream << "{" << std::endl;
    
    _write_os_factors(fstream, signature);
    
    _write_distances(fstream, signature);
    
    _write_buffers(fstream, signature, true);
    
    _write_fractions(fstream, recgroup);
    
    _write_hrr_loop(fstream, recgroup);
    
    fstream << "}" << std::endl << std::endl;
}


void
EriCPUGenerator::_write_os_factors(      std::ofstream&          fstream,
                                   const Signature<T4CIntegral>& signature) const
{
    const auto space = std::string(4, ' ');
    
    bool header = true;
    
    for (const auto& tfact : signature.factors())
    {
        if (const auto fname = _factor_name(tfact.name()); fname.find("rDistances") == std::string::npos)
        {
            if (header)
            {
                fstream << space << "// set up Obara-Saika factors";
                
                fstream << std::endl << std::endl;
                
                header = false;
            }
            
            fstream << space << "auto " << tfact.label() << " = " << fname << ";";
            
            fstream << std::endl << std::endl;
        }
    }
}

void
EriCPUGenerator::_write_distances(      std::ofstream&          fstream,
                                  const Signature<T4CIntegral>& signature) const
{
    const auto space = std::string(4, ' ');
    
    for (const auto& tname : signature.factor_names())
    {
        if (const auto label = _factor_name(tname); label.find("rDistances") != std::string::npos)
        {
            fstream << space << "// set up R(" << tname << ") distances";
                
            fstream << std::endl << std::endl;
        
            for (const auto& tcomp : signature.factors(tname))
            {
                const auto tlabel = tcomp.label();
            
                fstream << space << "auto " << tlabel << " = ";
            
                fstream << label << ".data(";
            
                if (tlabel.find("_x") != std::string::npos) fstream << "0";
            
                if (tlabel.find("_y") != std::string::npos) fstream << "1";
        
                if (tlabel.find("_z") != std::string::npos) fstream << "2";
            
                fstream << ");" << std::endl << std::endl;
            }
        }
    }
}

void
EriCPUGenerator::_write_buffers(      std::ofstream&          fstream,
                                const Signature<T4CIntegral>& signature,
                                const bool                    flg_hrr) const
{
    // base integral components
    
    _write_intetgrals(fstream, signature.params("out"), flg_hrr);
    
    // recursion integral components
    
    for (const auto& tint : signature.expansion<I4CIntegral>())
    {
        _write_intetgrals(fstream, signature.expansion_components(tint), flg_hrr);
    }
}

void
EriCPUGenerator::_write_intetgrals(      std::ofstream&         fstream,
                                   const std::set<T4CIntegral>& integrals,
                                   const bool                   flg_hrr) const
{
    const auto space = std::string(4, ' ');
    
    bool header = true;
    
    int32_t idx = 0;
    
    for (const auto& tcomp : integrals)
    {
        const auto tint = I4CIntegral(tcomp);
        
        if (header)
        {
            if (flg_hrr)
            {
                fstream << space << "// set up (" << tint.label() << ") ";
                
                fstream << "integral components";
            }
            else
            {
                fstream << space << "// set up [" << tint.label() << "]^(";
                
                fstream << std::to_string(tcomp.order()) << ") integral components";
            }
            
            fstream << std::endl << std::endl;
            
            header = false;
        }
        
        fstream << space << "t_" << tcomp.label(!flg_hrr) << " = " ;
        
        if ((tint[0] + tint[1] + tint[2] + tint[3]) == 0)
        {
            fstream << _buffer_name(tint, flg_hrr) << ";";
        }
        else
        {
            fstream << _buffer_name(tint, flg_hrr) << ".data(";
            
            fstream << _indexes_name(tint, flg_hrr) << "(";
            
            fstream << std::to_string(idx) << "));";
        }
        
        fstream << std::endl << std::endl;
        
        idx++;
    }
}

void
EriCPUGenerator::_write_fractions(      std::ofstream& fstream,
                                  const R4Group&       recgroup) const
{
    const auto space = std::string(4, ' ');
    
    bool header = true;
    
    for (const auto& tval : recgroup.prefactors())
    {
        if ((tval == Fraction(1)) || (tval == Fraction(-1))) continue;
        
        if (header)
        {
            fstream << space << "// set up scaling factors";
            
            fstream << std::endl << std::endl;
            
            header = false;
        }
        
        fstream << space << "const auto " << _fraction_name(tval);
        
        fstream << " = static_cast<T>(" << tval.label() << ");";
        
        fstream << std::endl << std::endl;
    }
}

void
EriCPUGenerator::_write_vrr_loop(      std::ofstream& fstream,
                                 const R4Group&       recgroup) const
{
    const auto ncomps = static_cast<int32_t>(recgroup.expansions());
    
    auto ngroups = ncomps / 36;
    
    if ((ncomps % 36) != 0) ngroups++;
    
    const auto space = std::string(4, ' ');
    
    fstream << space << "if (useSummation)" << std::endl;
        
    fstream << space << "{" << std::endl;
        
    for (int32_t i = 0; i < ngroups; i++)
    {
        const auto bstart = i * 36;
        
        const auto bend = ((bstart + 36) > ncomps) ? ncomps : bstart + 36;
        
        _write_simd_loop(fstream, recgroup, bstart, bend, false, true);
    }
    
    fstream << space << "}" << std::endl;
    
    fstream << space << "else" << std::endl;
    
    fstream << space << "{" << std::endl;
    
    for (int32_t i = 0; i < ngroups; i++)
    {
        const auto bstart = i * 36;
        
        const auto bend = ((bstart + 36) > ncomps) ? ncomps : bstart + 36;
        
        _write_simd_loop(fstream, recgroup, bstart, bend, false, false);
        
        if (bend != ncomps) fstream << std::endl;
    }
    
    fstream << space << "}" << std::endl;
}

void
EriCPUGenerator::_write_hrr_loop(      std::ofstream& fstream,
                                 const R4Group&       recgroup) const
{
    const auto ncomps = static_cast<int32_t>(recgroup.expansions());
    
    auto ngroups = ncomps / 36;
    
    if ((ncomps % 36) != 0) ngroups++;
    
    const auto space = std::string(4, ' ');
            
    for (int32_t i = 0; i < ngroups; i++)
    {
        const auto bstart = i * 36;
        
        const auto bend = ((bstart + 36) > ncomps) ? ncomps : bstart + 36;
        
        _write_simd_loop(fstream, recgroup, bstart, bend, true, false);
        
        if (bend != ncomps) fstream << std::endl;
    }
}

void
EriCPUGenerator::_write_simd_loop(      std::ofstream& fstream,
                                  const R4Group&       recgroup,
                                  const int32_t        lstart,
                                  const int32_t        lend,
                                  const bool           flg_hrr, 
                                  const bool           flg_sum) const
{
    // omp header for loop
    
    _write_omp_header(fstream, recgroup, lstart, lend, flg_hrr);

    // loop body
    
    const auto space2x = (flg_hrr) ? std::string(4, ' ') : std::string(8, ' ');
    
    const auto space3x = (flg_hrr) ? std::string(8, ' ') : std::string(12, ' ');
    
    fstream << space2x << "for (int32_t i = 0; i < nBatchPairs; i++)" << std::endl;
    
    fstream << space2x << "{" << std::endl;
    
    for (int32_t i = lstart; i < lend; i++)
    {
        // reference integral
        
        fstream << space3x << _rec_term_name(recgroup[i].root(), "[i]", true, flg_hrr);
        
        fstream << ((flg_sum) ? " += " : " = ");
        
        // recursion terms
        
        const auto nterms = recgroup[i].terms();
        
        for (int32_t j = 0; j < nterms; j++)
        {
            fstream << _rec_term_name((recgroup[i])[j], "[i]", j == 0, flg_hrr);
            
            if ((j + 1) != nterms) fstream << " ";
        }
        
        fstream << ";" << std::endl;
        
        if ((i + 1) != lend) fstream << std::endl;
    }
    
    fstream << space2x << "}" << std::endl;
}

void
EriCPUGenerator::_write_omp_header(      std::ofstream& fstream,
                                   const R4Group&       recgroup,
                                   const int32_t        lstart,
                                   const int32_t        lend,
                                   const bool           flg_hrr) const
{
    const auto labels = _get_align_vars(recgroup, lstart, lend, flg_hrr);
    
    std::string vstr = (flg_hrr) ? std::string(4, ' ') : std::string(8, ' ');
    
    vstr += "#pragma omp simd align(";
    
    for (const auto& tlabel : labels)
    {
        vstr += tlabel;
        
        if (vstr.size() > 81)
        {
            if (tlabel != *labels.rbegin())
            {
                fstream << vstr << ",\\" << std::endl;
                
                vstr = (flg_hrr) ? std::string(27, ' ') : std::string(31, ' ');
            }
            else
            {
                fstream << vstr << " : VLX_ALIGN)" << std::endl;
                
                vstr.clear();
            }
        }
        else
        {
            if (tlabel != *labels.rbegin()) vstr += ", ";
        }
    }

    if (!vstr.empty())
    {
        fstream << vstr << " : VLX_ALIGN)" << std::endl;
    }
}

std::set<std::string>
EriCPUGenerator::_get_align_vars(const R4Group& recgroup,
                                 const int32_t  lstart,
                                 const int32_t  lend,
                                 const bool     flg_hrr) const
{
    std::set<std::string> avars;
    
    for (int32_t i = lstart; i < lend; i++)
    {
        // reference integral
        
        avars.insert("t_" + (recgroup[i].root()).label(!flg_hrr));
    
        // recursion terms
        
        const auto nterms = recgroup[i].terms();
        
        for (int32_t j = 0; j < nterms; j++)
        {
            const auto rterm = (recgroup[i])[j];
            
            avars.insert("t_" + rterm.label(!flg_hrr));
            
            // factors of recursion term
            
            for (const auto& tfact : rterm.factors())
            {
                avars.insert(tfact.label());
            }
        }
    }
    
    return avars;
}

std::string
EriCPUGenerator::_hvrr_func_name(const std::map<Signature<T4CIntegral>, R4Group>& signatures,
                                const Signature<T4CIntegral>&                    signature) const
{
    std::string label = "compHost";
 
    const auto tint = signature.base<I4CIntegral>();
    
    if (_is_vrr_rec(*tint))
    {
        label += "VRR";
    }
    else
    {
        label += "HRR";
    }
    
    label += "For" + tint->label();
        
    int32_t idx = 0;
    
    for (const auto& tval : signatures)
    {
        if (tval.first == signature)
        {
            return label + "_V" + std::to_string(idx);
        }
        
        idx++;
    }
    
    return std::string();
}



