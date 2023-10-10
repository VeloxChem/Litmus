#ifndef t4c_full_body_hpp
#define t4c_full_body_hpp

#include <string>
#include <vector>
#include <fstream>

#include "t4c_defs.hpp"
#include "file_stream.hpp"

// Four-center compute function body generators for CPU.
class T4CFullFuncBodyDriver
{
    /// Generates vector of strings with GTOs pairs definitions in compute function.
    /// @return The vector of strings with GTOS definitions in compute function.
    std::vector<std::string> _get_gtos_def() const;
    
    /// Generates vector of strings with variables definitions in compute function.
    /// @param integral The base four center integral.
    /// @return The vector of strings with variables definitions in compute function.
    std::vector<std::string> _get_vars_def(const I4CIntegral& integral) const;
    
    /// Generates vector of strings with main loop definition in compute function.
    /// @return The vector of strings with main loop definition in compute function.
    std::vector<std::string> _get_batches_def() const;
    
    /// Adds loop start definitions to code lines container.
    /// @param lines The code lines container to which loop start definition are added.
    void _add_batches_loop_start(VCodeLines& lines) const;
    
    /// Adds loop body definitions to code lines container.
    /// @param lines The code lines container to which loop body definition are added.
    /// @param integral The base four center integral.
    void _add_batches_loop_body(      VCodeLines&  lines,
                                const I4CIntegral& integral) const;
    
    /// Adds loop end definitions to code lines container.
    /// @param lines The code lines container to which loop end definition are added.
    void _add_batches_loop_end(VCodeLines& lines) const;
    
    /// Adds loop body definitions to code lines container.
    /// @param lines The code lines container to which loop body definition are added.
    /// @param integral The base four center integral.
    /// @param component The base four center integral component.
    void _add_component_body(      VCodeLines&  lines,
                             const I4CIntegral& integral,
                             const T4CIntegral& component) const;
    
    /// Adds definition of block distribution call tree for compute function.
    /// @param lines The code lines container to which primitives loop start definition are added.
    /// @param integral The base two center integral.
    /// @param component The base four center integral component.
    void _write_block_distributor(      VCodeLines&  lines,
                                  const I4CIntegral& integral,
                                  const T4CIntegral& component) const;
    
public:
    /// Creates a four-center compute function body generator.
    T4CFullFuncBodyDriver() = default;
    
    /// Writes body of compute function.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void write_func_body(      std::ofstream& fstream,
                         const I4CIntegral&   integral) const;
};

#endif /* t4c_full_body_hpp */
