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

#ifndef eri_cpu_generator_hpp
#define eri_cpu_generator_hpp

#include <set>
#include <fstream>

#include "eri_driver.hpp"
#include "repository.hpp"

/// Electron repulsion integrals code generator for CPU class.
class EriCPUGenerator
{
    // The diagonal form flag.
    bool _diag_form;
    
    /// Writes header file for VRR recursion.
    /// @param integral The base four center integral.
    /// @param rectype The recursion type.
    /// @return The file name.
     std::string _file_name(const I4CIntegral& integral,
                            const std::string& rectype) const;
    
    /// Writes header file for VRR recursion.
    /// @param integral The base four center integral.
    /// @param repo The repository of recursion graphs.
    void _write_vrr_cpp_header(const I4CIntegral&                      integral,
                               const Repository<R4Group, T4CIntegral>& repo) const;
    
    /// Writes header file for HRR recursion.
    /// @param integral The base four center integral.
    /// @param repo The repository of recursion graphs.
    void _write_hrr_cpp_header(const I4CIntegral&                      integral,
                               const Repository<R4Group, T4CIntegral>& repo) const;
    
    /// Writes header file for computation scheme of recursion.
    /// @param graph The recursion graph.
    void _write_comp_cpp_header(const Graph<R4Group>* graph) const;
    
    /// Gets buffer name for given integral
    /// @param integral The base four center integral.
    /// @param flg_hrr The HRR recursion flag.
    /// @return The buffer name.
    std::string _buffer_name(const I4CIntegral& integral,
                             const bool         flg_hrr) const;
    
    /// Gets buffer's indexes name for given integral
    /// @param integral The base four center integral.
    /// @param flg_hrr The HRR recursion flag.
    /// @return The buffer's indexes name.
    std::string _indexes_name(const I4CIntegral& integral,
                              const bool         flg_hrr) const;
    
    /// Gets factor name for given factor.
    /// @param label The label of factor.
    /// @return The factor name.
    std::string _factor_name(const std::string& label) const;
    
    /// Gets fraction name for given fraction.
    /// @param fraction The fraction.
    /// @return The fraction name.
    std::string _fraction_name(const Fraction& fraction) const;
    
    /// Gets fraction name for given fraction.
    /// @param graph The recursion graph. 
    /// @return The number of Obara-Saika factors.
    int _number_os_factors(const Graph<R4Group>* graph) const;
    
    /// Gets recursion term name for given fraction.
    /// @param recterm The recursion term.
    /// @param index The string with index.
    /// @param first The flag for first recursion term in expansion.
    /// @param flg_hrr The HRR recursion flag.
    /// @return The  name of recursion term.
    std::string _rec_term_name(const R4CTerm&     recterm,
                               const std::string& index,
                               const bool         first,
                               const bool         flg_hrr) const;
    
    /// Checks if integral is generated by horizontal recursion.
    /// @param integral The base four center integral.
    /// @return True  if integral is generated by horizontal recursion, false otherwise.
    bool _is_hrr_rec(const I4CIntegral& integral) const;
    
    /// Checks if integral is generated by vertical recursion.
    /// @param integral The base four center integral.
    /// @return True if integral is generated by vertical recursion, false otherwise.
    bool _is_vrr_rec(const I4CIntegral& integral) const;
    
    /// Checks if integral is auxilary integral.
    /// @param integral The base four center integral.
    /// @return True if auxilary false otherwise.
    bool _is_aux_rec(const I4CIntegral& integral) const;
    
    /// Checks if name of factor is distance factor name.
    /// @param name The name of factor.
    /// @return True if name is name of distance factor false otherwise.
    bool _is_distance(const std::string& name) const;
    
    /// Writes HRR/VRR  function declaration to file stream.
    /// @param fstream the file stream.
    /// @param signatures the list of signatures.
    /// @param signature the signature of integral.
    void _write_hvrr_func_decl(      std::ofstream&                             fstream,
                              const std::map<Signature<T4CIntegral>, R4Group>& signatures,
                              const Signature<T4CIntegral>&                    signature) const;
    
    /// Writes compute function declaration to file stream.
    /// @param fstream the file stream.
    /// @param graph The recursion graph. 
    void _write_comp_func_decl(      std::ofstream&  fstream,
                               const Graph<R4Group>* graph) const;
    
    /// Writes VRR  function body to file stream.
    /// @param fstream the file stream.
    /// @param signature the signature of integral.
    /// @param recgroup the recursion group.
    void _write_vrr_func_body(      std::ofstream&          fstream,
                              const Signature<T4CIntegral>& signature,
                              const R4Group&                recgroup) const;
    
    /// Writes HRR  function body to file stream.
    /// @param fstream the file stream.
    /// @param signature the signature of integral.
    /// @param recgroup the recursion group.
    void _write_hrr_func_body(      std::ofstream&          fstream,
                              const Signature<T4CIntegral>& signature,
                              const R4Group&                recgroup) const;
    
    /// Writes compute function body to file stream.
    /// @param fstream the file stream.
    /// @param graph The recursion graph. 
    void _write_comp_func_body(      std::ofstream&  fstream,
                               const Graph<R4Group>* graph) const;
    
    /// Writes Obara-Saika factors to file stream.
    /// @param fstream the file stream.
    /// @param signature the signature of integral.
    void _write_os_factors(      std::ofstream&          fstream,
                           const Signature<T4CIntegral>& signature) const;
    
    /// Writes distances used in recursion to file stream.
    /// @param fstream the file stream.
    /// @param signature the signature of integral.
    void _write_distances(      std::ofstream&          fstream,
                          const Signature<T4CIntegral>& signature) const;
    
    /// Writes compute factors to file stream.
    /// @param fstream the file stream.
    /// @param graph The recursion graph. 
    void _write_comp_factors(      std::ofstream&  fstream,
                             const Graph<R4Group>* graph) const;
    
    /// Writes integral buffers used in recursion to file stream.
    /// @param fstream the file stream.
    /// @param signature the signature of integral.
    /// @param flg_hrr  the HRR recursion flag.
    void _write_buffers(      std::ofstream&          fstream,
                        const Signature<T4CIntegral>& signature,
                        const bool                    flg_hrr) const;
    
    /// Writes integral components to file stream.
    /// @param fstream the file stream.
    /// @param integrals the set of integral components.
    /// @param flg_hrr  the HRR recursion flag.
    void _write_intetgrals(      std::ofstream&         fstream,
                           const std::set<T4CIntegral>& integrals,
                           const bool                   flg_hrr) const;
    
    /// Writes fraction factors to file stream.
    /// @param fstream the file stream.
    /// @param recgroup the recursion group to generate integral components.
    void _write_fractions(      std::ofstream& fstream,
                          const R4Group&       recgroup) const;
    
    /// Writes VRR recursion loop to file stream.
    /// @param fstream the file stream.
    /// @param recgroup the recursion group to generate integral components.
    void _write_vrr_loop(      std::ofstream& fstream,
                         const R4Group&       recgroup) const;
    
    /// Writes HRR recursion loop to file stream.
    /// @param fstream the file stream.
    /// @param recgroup the recursion group to generate integral components.
    void _write_hrr_loop(      std::ofstream& fstream,
                         const R4Group&       recgroup) const;
    
    /// Writes recursion loop to file stream.
    /// @param fstream the file stream.
    /// @param recgroup the recursion group to generate integral components.
    /// @param lstart the start of partial recursion loop.
    /// @param lend the end of partial recursion loop.
    /// @param flg_hrr  the HRR recursion flag.
    /// @param flg_sum  the summation form of recursion.
    void _write_simd_loop(      std::ofstream& fstream,
                          const R4Group&       recgroup,
                          const int32_t        lstart,
                          const int32_t        lend,
                          const bool           flg_hrr,
                          const bool           flg_sum) const;
    
    /// Writes omp header for simd loop to file stream.
    /// @param fstream the file stream.
    /// @param recgroup the recursion group to generate integral components.
    /// @param lstart the start of partial recursion loop.
    /// @param flg_hrr  the HRR recursion flag.
    void _write_omp_header(      std::ofstream& fstream,
                           const R4Group&       recgroup,
                           const int32_t        lstart,
                           const int32_t        lend,
                           const bool           flg_hrr) const;
    
    /// Gets set of unique variables in given range of recursion group.
    /// @param recgroup the recursion group to generate integral components.
    /// @param lstart the start of partial recursion loop.
    /// @param lend the end of partial recursion loop.
    /// @param flg_hrr  the HRR recursion flag.
    /// @return The set of unique variables.
    std::set<std::string> _get_align_vars(const R4Group& recgroup,
                                          const int32_t  lstart,
                                          const int32_t  lend,
                                          const bool     flg_hrr) const;
    
    /// Gets HRR/VRR  function name.
    /// @param signatures the list of signatures.
    /// @param signature the signature of integral.
    std::string _hvrr_func_name(const std::map<Signature<T4CIntegral>, R4Group>& signatures,
                                const Signature<T4CIntegral>&                    signature) const;
    
    
    /// Writes diagonal VRR  includes to file stream.
    /// @param fstream the file stream.
    /// @param graph The recursion graph.
    void _write_diag_includes(      std::ofstream&  fstream,
                              const Graph<R4Group>* graph) const;
    
public:
    /// Creates an electron repulsion integrals CPU code generator.
    EriCPUGenerator();
 
    /// Sets diagonal form of generated integrals.
    void set_diag_form(); 
    
    /// Generates electron repulsion integrals code for the given repository.
    /// @param repo The repository of two-electron integrals.
    void generate(const Repository<R4Group, T4CIntegral>& repo) const;
};

#endif /* eri_cpu_generator_hpp */
