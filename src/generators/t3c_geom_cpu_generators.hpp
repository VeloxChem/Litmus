#ifndef t3c_geom_cpu_generators_hpp
#define t3c_geom_cpu_generators_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <array>
#include <utility>

#include "t3c_defs.hpp"

// Geometrical derivatives of three-center integrals code generator for CPU.
class T3CGeomCPUGenerator
{
    /// Checks if recursion is available for four-center inetgral with given label.
    /// @param label The label of requested four-center integral.
    bool _is_available(const std::string& label) const;
    
    /// Gets four-center inetgral with requested label.
    /// @param label The label of requested four-center integral.
    /// @param ang_moms The angular momentum of  A, C, and D centers.
    /// @param geom_drvs The geometrical derivative of bra side, integrand, and  ket side.
    /// @return The four-center integral.
    I3CIntegral _get_integral(const std::string&        label,
                              const std::array<int, 3>& ang_moms,
                              const std::array<int, 3>& geom_drvs) const;
    
    /// Generates set of integrals required for geometrical derivatives.
    /// @param integral The base three center integral.
    /// @return The set of integrals.
    SI3CIntegrals _generate_geom_integral_group(const I3CIntegral& integral) const;
    
    /// Generates set of geometrical terms required for geometrical derivatives.
    /// @param integrals The set of four center integrals.
    /// @return The set of geometrical terms.
    SG3Terms _generate_geom_terms_group(const SI3CIntegrals& integrals,
                                        const I3CIntegral&   integral) const;
    
    /// Adds ket horizontal recursion to geometrical terms.
    /// @param terms The set of geometrical terms.
    void _add_ket_hrr_terms_group(SG3Terms& terms) const;
    
    /// Filters cbuffer terms from set of geometrical terms.
    /// @param terms The set of filtered geometrical terms.
    SG3Terms _filter_cbuffer_terms(const SG3Terms& terms) const;
    
    /// Filters ckbuffer terms from set of geometrical terms.
    /// @param terms The set of filtered geometrical terms.
    SG3Terms _filter_skbuffer_terms(const I3CIntegral& integral,
                                    const SG3Terms&    terms) const;
    
    /// Generates set of integrals required for geometrical derivatives.
    /// @param terms The set of geometrical terms.
    /// @return The set of integrals.
    SI3CIntegrals _generate_vrr_integral_group(const SG3Terms& terms) const;
    
    /// Writes header file for recursion.
    /// @param cterms The set of filtered geometrical terms.
    /// @param skterms The set of filtered geometrical terms.
    /// @param vrr_integrals The set of unique integrals for vertical recursion.
    /// @param integral The base two center integral.
    void _write_cpp_header(const SG3Terms&      cterms,
                           const SG3Terms&      skterms,
                           const SI3CIntegrals& vrr_integrals,
                           const I3CIntegral&   integral) const;
    
    /// Gets file name of file with recursion functions for two center integral.
    /// @param integral The base two center integral.
    /// @return The file name.
    std::string _file_name(const I3CIntegral& integral) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_defines(      std::ofstream& fstream,
                            const I3CIntegral&   integral,
                            const bool           start) const;
    
    /// Writes definitions of includes for header file.
    /// @param fstream the file stream.
    /// @param skterms The set of filtered geometrical terms.
    /// @param vrr_integrals The set of unique integrals for vertical recursion.
    /// @param integral The base two center integral.
    void _write_hpp_includes(      std::ofstream& fstream,
                             const SG3Terms&      skterms,
                             const SI3CIntegrals& vrr_integrals,
                             const I3CIntegral&   integral) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const I3CIntegral&   integral,
                          const bool           start) const;
    
    /// Prunes set of geometrical terms removing terms which matches one to one with ordinary integrals.
    /// @param terms The set of geometrical terms.
    void _prune_terms_group(SG3Terms& terms) const;
    
public:
    /// Creates a geometrical derivatives of three-center integrals CPU code generator.
    T3CGeomCPUGenerator() = default;
     
    /// Generates selected four-center integrals up to given angular momentum (inclusive)  on A, B, C, and D centers.
    /// @param label The label of requested two-center integral.
    /// @param max_ang_mom The maximum angular momentum of A, B, C and D centers.
    /// @param max_aux_ang_mom The maximum angular momentum of auxilary A center.
    /// @param geom_drvs The geometrical derivative of bra side, integrand, and  ket side.
    void generate(const std::string&        label,
                  const int                 max_ang_mom,
                  const int                 max_aux_ang_mom,
                  const std::array<int, 3>& geom_drvs) const;
};


#endif /* t3c_geom_cpu_generators_hpp */
