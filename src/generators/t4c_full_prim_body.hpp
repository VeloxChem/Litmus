#ifndef t4c_full_prim_body_hpp
#define t4c_full_prim_body_hpp

#include <string>
#include <vector>
#include <fstream>

#include "t4c_defs.hpp"
#include "file_stream.hpp"

// Four-center primitive compute function body generators for CPU.
class T4CFullPrimFuncBodyDriver
{
    /// Generates the common data definitions in primitive compute function.
    /// @return The vector of common primitives  data.
    std::vector<std::string> _get_common_data_str() const;
    
    /// Generates the common data definitions in primitive compute function.
    /// @return The vector of common primitives  data.
    std::vector<std::string> _get_vrr_common_data_str() const;
    
    /// Generates the common data definitions in primitive compute function.
    /// @return The vector of common primitives  data.
    std::vector<std::string> _get_hrr_common_data_str() const;
    
    /// Adds block of code for computation of P  and Q coordinates.
    /// @param lines The code lines container to which simd code are added.
    void _add_coords_compute(VCodeLines&  lines) const;
    
    /// Generates vector of Boys function variables strings.
    /// @param integral The base two center integral.
    /// @return The vector of special variable strings.
    std::vector<std::string> _get_boys_vars_str(const I4CIntegral& integral) const;
    
    /// Adds Boys function computation code lines.
    /// @param lines The code lines container to which simd code are added.
    /// @param integral The base two center integral.
    void _add_boys_compute_lines(      VCodeLines&  lines,
                                 const I4CIntegral& integral) const;
    
    /// Adds split simd code generated for selected set of integral components.
    /// @param lines The code lines container to which simd code are added.
    /// @param component The integral component.
    /// @param integral The base integral.
    void _add_split_simd_code(      VCodeLines&  lines,
                              const T4CIntegral& component,
                              const I4CIntegral& integral) const;
    
    /// Adds split simd code generated for selected set of integral components.
    /// @param lines The code lines container to which simd code are added.
    /// @param component The integral component.
    /// @param integral The base integral.
    void _add_split_hrr_simd_code(      VCodeLines&  lines,
                                  const T4CIntegral& component,
                                  const I4CIntegral& integral) const;
    
    /// Adds block of  simd code lines for given reccursion distribution.
    /// @param lines The code lines container to which simd code are added.
    /// @param integral The integral  component.
    /// @param rdist  The recursion distribution.
    void _add_split_simd_block(      VCodeLines&  lines,
                               const T4CIntegral& integral,
                               const R4CDist&     rdist) const;
    
    /// Adds split loop pragma definition.
    /// @param lines The code lines container to which simd code are added.
    /// @param integral The integral  component.
    /// @param rdist  The recursion distribution.
    void _add_split_pragma(      VCodeLines&  lines,
                           const T4CIntegral& integral,
                           const R4CDist&     rdist) const;
    
    /// Adds split loop start definition.
    /// @param lines The code lines container to which simd code are added.
    /// @param integral The integral  component.
    /// @param rdist  The recursion distribution.
    void _add_split_loop_start(      VCodeLines&  lines,
                               const T4CIntegral& integral,
                               const R4CDist&     rdist) const;
    
    /// Adds split loop end definition.
    /// @param lines The code lines container to which simd code are added.
    /// @param integral The integral  component.
    void _add_split_loop_end(      VCodeLines&  lines,
                             const T4CIntegral& integral) const;
    
    /// Adds block of  simd code lines for given reccursion distribution.
    /// @param lines The code lines container to which simd code are added.
    /// @param integral The integral  component.
    /// @param rdist  The recursion distribution.
    /// @param index The unique integral index.
    void _add_simd_lines_block(      VCodeLines&  lines,
                               const T4CIntegral& integral,
                               const R4CDist&     rdist,
                               const size_t       index) const;
    
public:
    /// Creates a diagonal four-center compute function body generator.
    T4CFullPrimFuncBodyDriver() = default;
    
    /// Writes body of primitive compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    void write_prim_func_body(      std::ofstream& fstream,
                              const T4CIntegral&   component,
                              const I4CIntegral&   integral) const;
    
    /// Writes body of primitive VRR compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    void write_vrr_func_body(      std::ofstream& fstream,
                             const T4CIntegral&   component,
                             const I4CIntegral&   integral) const;
    
    /// Writes body of primitive HRR compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    void write_hrr_func_body(      std::ofstream& fstream,
                             const T4CIntegral&   component,
                             const I4CIntegral&   integral) const;
};


#endif /* t4c_full_prim_body_hpp */
