#ifndef t1c_body_hpp
#define t1c_body_hpp

#include <string>
#include <vector>
#include <fstream>

#include "t2c_defs.hpp"
#include "file_stream.hpp"

// GTOs compute function body generators for CPU.
class T1CFuncBodyDriver
{
    /// Generates vector of strings with angular momentum transformation factors.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @return The vector of strings with angular momentum transformation factors.
    std::vector<std::string> _get_angmom_def(const int angmom) const; 
    
    /// Generates vector of strings with GTOs definitions in compute function.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    /// @return The vector of strings with GTOS definitions in compute function.
    std::vector<std::string> _get_gtos_def(const int angmom,
                                           const int gdrv) const;
    
    /// Adds bra loop body definitions to code lines container.
    /// @param lines The code lines container to which bra loop body definition are added.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    void _add_loop_body(      VCodeLines&  lines,
                        const int          angmom,
                        const int          gdrv) const;
    
    /// Adds bra loop body definitions to code lines container.
    /// @param lines The code lines container to which bra loop body definition are added.
    /// @param integral The base two center integral.
    /// @param angcomp The angular momentum componentn of Cartesian GTO.
    /// @param gdrv The geometrical derivative of GTOs.
    void _add_simd_code(      VCodeLines&      lines,
                        const I2CIntegral&     integral,
                        const TensorComponent& angcomp,
                        const int              gdrv) const;
    
    /// Generates  the recursion  group for given vector of integral components.
    /// @param components The vector of integral components.
    /// @return The recursion group.
    R2Group _generate_integral_group(const VT2CIntegrals& components) const;
    
    /// Creates integral representation of geometrical derivatives of GTOs.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    I2CIntegral _get_integral(const int angmom,
                              const int gdrv) const;
    
    /// Select integrals components with predefined components.
    /// @param component the tensor component.
    /// @param integral The base two center integral.
    VT2CIntegrals _select_integral_components(const TensorComponent& component,
                                              const I2CIntegral&     integral) const;
    
    /// Adds simd line for given recursion distribution.
    /// @param lines The code lines container to which bra loop body definition are added.
    /// @param rdist The recursion distribution.
    void _add_simd_line(      VCodeLines&  lines,
                        const R2CDist&     rdist) const;
    
   /// Gets polynomial representation of Cartesian GTO component.
   /// @param component the tensor component.
   std::string _polynomial_string(const TensorComponent& component) const;
    
    /// Adds GTOs values distribution code.
    /// @param lines The code lines container to which bra loop body definition are added.
    /// @param component The angular momentum componentn of Cartesian GTO.
    /// @param gdrv The geometrical derivative of GTOs.
    void _add_distribution_code(      VCodeLines&      lines,
                                const TensorComponent& component,
                                const int              gdrv) const;
  
public:
    /// Creates a GTOs compute function body generator.
    T1CFuncBodyDriver() = default;
    
    /// Writes body of compute function.
    /// @param fstream the file stream.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    void write_func_body(      std::ofstream& fstream,
                         const int            angmom,
                         const int            gdrv) const;
};


#endif /* t1c_body_hpp */
