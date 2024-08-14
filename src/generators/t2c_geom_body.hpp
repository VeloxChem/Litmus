#ifndef t2c_geom_body_hpp
#define t2c_geom_body_hpp

#include <string>
#include <array>
#include <vector>
#include <utility>
#include <fstream>

#include "t2c_defs.hpp"
#include "file_stream.hpp"

// Two-center compute function body generators for CPU.
class T2CGeomFuncBodyDriver
{
    /// Generates vector of buffer strings.
    /// @param geom_integrals The set of unique integrals for geometrical recursion.
    /// @param integral The base four center integral.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_buffers_str(const SI2CIntegrals& geom_integrals,
                                              const I2CIntegral&   integral) const;
    
    /// Generates vector of buffer strings.
    /// @param integral The base four center integral.
    /// @param components The vector of integral components.
    /// @param rec_range The recursion range [first, last) in integral components space.
    /// @return The vector of buffer strings.
    std::vector<std::string> _get_buffers_str(const I2CIntegral&        integral,
                                              const VT2CIntegrals&      components,
                                              const std::array<int, 2>& rec_range) const;
    
    /// Gets tensor label for integral.
    /// @param integral The base four center integral.
    /// @return The tensorial label.
    std::string _get_tensor_label(const I2CIntegral& integral) const;
    
    /// Gets tensor label for integral.
    /// @param integral The base four center integral.
    /// @return The tensorial label.
    std::string _get_tensor_label(const T2CIntegral& integral) const;
    
    /// Adds single loop computation of primitive integrals.
    /// @param lines The code lines container to which loop start definition are added.
    /// @param rgroup The recursion group.
    /// @param integral The base four center integral.
    /// @param rec_range The recursion range [first, last) in integral components space.
    void _add_recursion_loop(      VCodeLines&         lines,
                             const R2Group&            rgroup,
                             const I2CIntegral&        integral,
                             const std::array<int, 2>& rec_range) const;

    /// Gets pragma string for vector of recursion distributions.
    /// @param rgroup The recursion group.
    /// @param integral The base four center integral.
    /// @param rec_range The recursion range [first, last) in integral components space.
    std::string _get_pragma_str(const R2Group&            rgroup,
                                const I2CIntegral&        integral,
                                const std::array<int, 2>& rec_range) const;
    
    /// Creates code line for recursion expansion.
    /// @param rec_distribution The recursion distribution
    /// @return The string with code line.
    std::string _get_code_line(const R2CDist& rec_distribution) const;
    
    /// Creates code string for recursion term.
    /// @param rec_term The recursion distribution.
    /// @param is_first The flag to indicate first term in recursion expnasion.
    /// @return The string with code term.
    std::string _get_rterm_code(const R2CTerm& rec_term,
                                const bool     is_first) const;
    
    /// Gets integral component label.
    /// @param integral The base four center integral component.
    /// @return The string with integral component label.
    std::string _get_component_label(const T2CIntegral& integral) const;
    
    /// Creates recursion group for all components of given integral.
    /// @param components The vector of integral components.
    /// @param integral The diff. integral.
    /// @return The recursion group.
    R2Group _generate_integral_group(const VT2CIntegrals& components,
                                     const I2CIntegral&   integral) const;
    
    std::vector<std::string> _get_factors_str(const I2CIntegral& integral) const;
    
    R2CDist _get_geom_recursion(const T2CIntegral& integral) const;
    
    std::vector<std::string> _get_buffers_str(const std::vector<R2CDist>& rec_dists,
                                              const I2CIntegral&          integral) const;
    
    bool _find_integral(const std::vector<R2CDist>& rec_dists,
                        const T2CIntegral&          integral) const; 
    
    void _add_recursion_loop(      VCodeLines&         lines,
                             const I2CIntegral&        integral,
                             const VT2CIntegrals&      components,
                             const std::array<int, 2>& rec_range) const;
    
    std::string _get_pragma_str(const I2CIntegral&          integral,
                                const std::vector<R2CDist>& rec_distributions) const;
    
    void  _get_factor_lines(                VCodeLines& lines,
                            const I2CIntegral&          integral,
                            const std::vector<R2CDist>& rec_distributions) const; 

public:
    /// Creates a two-center compute function body generator.
    T2CGeomFuncBodyDriver() = default;
    
    /// @param geom_integrals The set of unique integrals for geometrical recursion.
    /// @param integral The base four center integral.
    void write_func_body(      std::ofstream& fstream,
                         const SI2CIntegrals& geom_integrals,
                         const I2CIntegral&   integral) const;
};

#endif /* t2c_geom_body_hpp */
