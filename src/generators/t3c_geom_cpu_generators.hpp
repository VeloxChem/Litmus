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
