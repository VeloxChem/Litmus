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
    // generate VRR and HRR recursions
    
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
    
    // generate integrals computation codes
    
    for (auto tgraph : repo.graphs())
    {
        _write_comp_cpp_header(tgraph);
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
    
    ost::write_hvrr_includes(fstream);
    
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
    
    ost::write_hvrr_includes(fstream);
    
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

void
EriCPUGenerator::_write_comp_cpp_header(const Graph<R4Group>* graph) const
{
    const auto tint = graph->base<I4CIntegral>();
    
    if (_is_aux_rec(tint)) return;
    
    std::string fname = _file_name(tint, "") + ".hpp";
    
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    ost::write_copyright(fstream);
    
    _write_diag_includes(fstream, graph);
    
    ost::write_namespace(fstream, "derirec", true);
    
    _write_comp_func_decl(fstream, graph);
    
    _write_comp_func_body(fstream, graph);
    
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
    // distances PB, QD, WP, WQ, AB, CD
    
    if (_is_distance(label))
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

int
EriCPUGenerator::_number_os_factors(const Graph<R4Group>* graph) const
{
    // get list of unique labels
    
    std::set<std::string> labels;
    
    for (const auto& tfact : graph->factors())
    {
        labels.insert(tfact.name());
    }
    
    int nfacts = 1;
    
    for (const auto& tlabel : labels)
    {
        // Obara-Saika factors
        
        if (tlabel == "1/zeta") nfacts++;
        
        if (tlabel == "1/eta") nfacts++;
        
        if (tlabel == "1/(zeta+eta)") nfacts++;
        
        if (tlabel == "rho/zeta^2") nfacts++;
        
        if (tlabel == "rho/eta^2") nfacts++;
    }
    
    return nfacts;
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

bool
EriCPUGenerator::_is_distance(const std::string& name) const
{
    if (name == "PQ") return true;
    
    if (name == "PB") return true;
    
    if (name == "QD") return true;
    
    if (name == "WP") return true;
    
    if (name == "WQ") return true;
    
    if (name == "AB") return true;
    
    if (name == "CD") return true;
    
    return false;
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
EriCPUGenerator::_write_comp_func_decl(      std::ofstream&  fstream,
                                       const Graph<R4Group>* graph) const
{
    // unique integrals
    
    const auto tint = graph->base<I4CIntegral>();
    
    // write function declaration
    
    fstream << "template <typename T>" << std::endl << "auto" << std::endl;
    
    fstream << "compHost" << tint.label() << "(" << std::string(6, ' ');
    
    fstream << "T*                                 intsBuffer,"  << std::endl;
    
    const auto space = std::string(13, ' ');
    
    fstream << space << "const CBinnedGtoPairBlock<T, mem::Host>* gtoPairBlock," << std::endl;
    
    fstream << space << "const int32_t                            bPosition," << std::endl;
    
    fstream << space << "const int32_t                            ePosition) -> void" << std::endl;
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
EriCPUGenerator::_write_comp_func_body(      std::ofstream&  fstream,
                                       const Graph<R4Group>* graph) const
{
    fstream << "{" << std::endl;
    
    ost::write_dimensions(fstream);
    
    _write_comp_factors(fstream, graph);
    
    _write_comp_buffers(fstream, graph);
    
    _write_comp_loop(fstream, graph); 
    
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
EriCPUGenerator::_write_comp_factors(      std::ofstream&  fstream,
                                     const Graph<R4Group>* graph) const
{
    const auto space = std::string(4, ' ');
    
    // write Obara-Saika factors
    
    if (const auto nfacts = _number_os_factors(graph); nfacts > 0)
    {
        fstream << space << "// allocate Obara-Saika factors" << std::endl << std::endl;
        
        fstream << space << "BufferHostMY<T, " << std::to_string(nfacts);
        
        fstream << "> osfacts(ncpairs);" << std::endl << std::endl;
    }
    
    // write distances
 
    fstream << space << "// allocate distances" << std::endl << std::endl;
    
    fstream << space << "BufferHostMY<T, 3> rpq(ncpairs); " << std::endl << std::endl;
    
    std::set<std::string> labels;
    
    for (const auto& tfact : graph->factors())
    {
        if (const auto name = tfact.name(); _is_distance(name))
        {
            if (labels.find(name) == labels.end())
            {
                fstream << space << "BufferHostMY<T, 3> " << tfact.label(true);
                
                fstream << "(ncpairs);" << std::endl << std::endl;
                
                labels.insert(name);
            }
        }
    }
    
    // write coordinates
    
    if ((labels.find("WP") != labels.end()) ||
        (labels.find("WQ") != labels.end()))
    {
        fstream << space << "// allocate coordinates" << std::endl << std::endl;
        
        fstream << space << "BufferHostMY<T, 3> rw(ncpairs); " << std::endl << std::endl;
    }
    
    // writes Boys function
    
    fstream << space << "// allocate Boys function data" << std::endl << std::endl;
    
    const auto tint = graph->base<I4CIntegral>();
    
    const auto border = tint[0] + tint[1] + tint[2] + tint[3];
    
    fstream << space << "BufferHostX<T> bargs(ncpairs);" << std::endl << std::endl;
    
    fstream << space << "BufferHostXY<T> bvals(" << border + 1;
    
    fstream << ", ncpairs);" << std::endl << std::endl;
    
    fstream << space << "CBoysFunc<T, " << border <<  "> bftable;";
    
    fstream << std::endl << std::endl;
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
EriCPUGenerator::_write_comp_buffers(      std::ofstream&  fstream,
                                     const Graph<R4Group>* graph) const
{
    const auto space = std::string(4, ' ');
    
    // auxilary integrals
    
    const auto rint = graph->base<I4CIntegral>();
    
    const auto border = rint[0] + rint[1] + rint[2] + rint[3] + 1;
    
    fstream << space << "// Primitive integral buffers" << std::endl << std::endl;
    
    fstream << space << "BufferHostXY<T> pbufSSSS(";
    
    fstream << border << ", ncpairs);" << std::endl << std::endl;
   
    // VRR integrals
    
    for (const auto& tint : _get_integrals(graph))
    {
        if (_is_vrr_rec(tint))
        {
            const auto tcomps = _get_components(tint, graph);
            
            if (rint == tint)
            {
                fstream << space << "// Contracted integral buffers" << std::endl << std::endl;
                
                fstream << space << "auto cbuf" << tint.label();
                
                fstream << " = BufferHostXY<T>::Zero(" << tcomps.size() << ", ncpairs);" << std::endl << std::endl;
            }
            else
            {
                fstream << space << "BufferHostXY<T> pbuf" << tint.label() << tint.order();
                
                fstream << "(" << tcomps.size() << ", ncpairs);" << std::endl << std::endl;
            }
        }
    }
    
    // HRR integrals
    
    if (_is_hrr_rec(rint))
    {
        fstream << space << "// Contracted integral buffers" << std::endl << std::endl;
    }
    
    for (const auto& tint : _get_hrr_integrals(graph))
    {
        const auto tcomps = _get_components(tint, graph);
        
        fstream << space << "BufferHostXY<T> cbuf" << tint.label();
            
        fstream << "(" << tcomps.size() << ", ncpairs);" << std::endl << std::endl;
    }
    
    for (const auto& tint : _get_integrals(graph))
    {
        if (_is_hrr_rec(tint))
        {
            const auto tcomps = _get_components(tint, graph);
            
            fstream << space << "BufferHostXY<T> cbuf" << tint.label();
                
            fstream << "(" << tcomps.size() << ", ncpairs);" << std::endl << std::endl;
        }
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
EriCPUGenerator::_write_comp_loop(      std::ofstream&  fstream,
                                  const Graph<R4Group>* graph) const
{
    
    const auto space = std::string(4, ' ');
    
    const auto space2x = std::string(8, ' ');
    
    const auto space3x = std::string(12, ' ');
    
    const auto space4x = std::string(16, ' ');
        
    // first loop pver primitive integrals
    
    fstream << space << "for (int32_t i = 0; i < nppairs; i++)" << std::endl;
    
    fstream << space << "{" << std::endl;
    
    // R(PB) distances
    
    if (_need_factor("PB", graph))
    {
        fstream << space2x;
        
        fstream << "derirec::compHostDistancesPT(rpb, gtoPairBlock, bPosition, ePosition, i);";
        
        fstream << std::endl << std::endl;
    }
    
    // R(QD) distances
    
    if (_need_factor("1/zeta", graph))
    {
        fstream << space2x;
        
        fstream << "derirec::compHostFactorsPartialZeta(fz, gtoPairBlock, bPosition, ePosition, i);";
        
        fstream << std::endl << std::endl;
    }
    
    // second loop over primitive integrals
    
    fstream << space2x << "for (int j = i; j < nppairs; j++)" << std::endl;
    
    fstream << space2x << "{" << std::endl;
    
    // R(PQ) distances
    
    fstream << space3x;
    
    fstream << "derirec::compHostDistancesPQ(rpq, gtoPairBlock, bPosition, ePosition, i, j);";
    
    fstream << std::endl << std::endl;
    
    // Obara-Saika factors: zeta*eta/(zeta+eta)
    
    fstream << space3x;
    
    fstream << "derirec::compHostFactorsRho(frho, gtoPairBlock, bPosition, ePosition, i, j);";
    
    fstream << std::endl << std::endl;
    
    // Obara-Saika factors: normalized overlaps
    
    fstream << space3x;
    
    fstream << "derirec::compHostFactorsNorm(fnorm, gtoPairBlock, bPosition, ePosition, i, j);";
    
    fstream << std::endl << std::endl;
    
    // Obara-Saika factors: 1/(zeta+eta)
    
    if (_need_factor("1/(zeta+eta)", graph))
    {
        fstream << space3x;
        
        fstream << "derirec::compHostFactorsZeta(fze, gtoPairBlock, bPosition, ePosition, i, j);";
        
        fstream << std::endl << std::endl;
    }
    
    // Obara-Saika factors: 1/eta
    
    if (_need_factor("1/eta", graph))
    {
        fstream << space3x;
        
        fstream << "derirec::compHostFactorsPartialZeta(fe, gtoPairBlock, bPosition, ePosition, j);";
        
        fstream << std::endl << std::endl;
    }
    
    // R(QD) distances
    
    if (_need_factor("QD", graph))
    {
        fstream << space3x;
        
        fstream << "derirec::compHostDistancesPT(rqd, gtoPairBlock, bPosition, ePosition, j);";
        
        fstream << std::endl << std::endl;
    }
    
    // W coordinates
    
    if (_need_factor("WP", graph) || _need_factor("WQ", graph))
    {
        fstream << space3x;
        
        fstream << "derirec::compHostCoordinatesW(rw, gtoPairBlock, bPosition, ePosition, i, j);";
        
        fstream << std::endl << std::endl;
    }
    
    // WP and WQ distances
    
    fstream << space3x << "if (i == j)" << std::endl;
    
    fstream << space3x << "{" << std::endl;
    
    if (_need_factor("WP", graph))
    {
        fstream << space4x << "rwp.setZero();" << std::endl;
        
        if (_need_factor("WQ", graph)) fstream << std::endl;
    }
    
    if (_need_factor("WQ", graph))
    {
        fstream << space4x << "rwq.setZero();" << std::endl;
    }
    
    fstream << space3x << "}" << std::endl;
    
    fstream << space3x << "else" << std::endl;
    
    fstream << space3x << "{" << std::endl;
    
    if (_need_factor("WP", graph))
    {
        fstream << space4x << "derirec::compHostDistancesWT(rwp, rw, gtoPairBlock, bPosition, ePosition, i);" << std::endl;
        
        if (_need_factor("WQ", graph)) fstream << std::endl;
    }
     
    if (_need_factor("WQ", graph))
    {
        fstream << space4x << "derirec::compHostDistancesWT(rwq, rw, gtoPairBlock, bPosition, ePosition, j);" << std::endl;
    }
        
    fstream << space3x << "}" << std::endl;
    
    fstream << space2x << "}" << std::endl;
    
    fstream << space << "}" << std::endl;
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

ST4CIntegrals
EriCPUGenerator::_get_components(const Graph<R4Group>* graph) const
{
    ST4CIntegrals tcomps;
    
    for (size_t i = 0; i < graph->vertices(); i++)
    {
        for (const auto& tval : (*graph)[i].components())
        {
            tcomps.insert(tval);
        }
    }
    
    return tcomps;
}

ST4CIntegrals
EriCPUGenerator::_get_components(const I4CIntegral&    integral,
                                 const Graph<R4Group>* graph) const
{
    ST4CIntegrals tcomps;
    
    for (size_t i = 0; i < graph->vertices(); i++)
    {
        for (const auto& tval : (*graph)[i].components())
        {
            if (integral == I4CIntegral(tval))
            {
                tcomps.insert(tval);
            }
        }
    }
    
    return tcomps;
}

std::set<I4CIntegral>
EriCPUGenerator::_get_integrals(const Graph<R4Group>* graph) const
{
    std::set<I4CIntegral> tints;
    
    for (const auto& tcomp : _get_components(graph))
    {
        tints.insert(I4CIntegral(tcomp));
    }
    
    return tints;
}

std::set<I4CIntegral>
EriCPUGenerator::_get_hrr_integrals(const Graph<R4Group>* graph) const
{
    std::set<I4CIntegral> tints;
    
    for (size_t i = 0; i < graph->vertices(); i++)
    {
        if (const auto tgraph = (*graph)[i]; _is_hrr_rec(*tgraph.base<I4CIntegral>()))
        {
            for (const auto& tcomp : tgraph.components())
            {
                if (const auto tint = I4CIntegral(tcomp); _is_vrr_rec(tint))
                {
                    tints.insert(tint);
                }
            }
        }
    }
    
    return tints;
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

void
EriCPUGenerator::_write_diag_includes(      std::ofstream&  fstream,
                                      const Graph<R4Group>* graph) const
{
    fstream << "#include <cstdint>" << std::endl << std::endl;
    
    fstream << "#include \"Buffer.hpp\"" << std::endl;
    
    fstream << "#include \"BinnedGtoPairBlock.hpp\"" << std::endl;
    
    fstream << "#include \"DiagEriRecFacts.hpp\"" << std::endl;
    
    // unique integrals
    
    const auto tints = graph->roots<I4CIntegral>();
    
    // VRR integrals
    
    for (const auto& tint : tints)
    {
        if (_is_vrr_rec(tint) || _is_aux_rec(tint))
        {
            fstream << "#include \"" + _file_name(tint, "VRR") + ".hpp\"" << std::endl;
        }
    }
    
    // HRR integrals
    
    for (const auto& tint : tints)
    {
        if (_is_hrr_rec(tint))
        {
            fstream << "#include \"" + _file_name(tint, "HRR") + ".hpp\"" << std::endl;
        }
    }
    
    fstream << std::endl;
}

bool
EriCPUGenerator::_need_factor(const std::string&    name,
                              const Graph<R4Group>* graph) const
{
    for (const auto& fact : graph->factors())
    {
        if (fact.name() == name)
        {
            return true;
        }
    }
    
    return false;
}
