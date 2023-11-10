#ifndef t4c_full_decl_hpp
#define t4c_full_decl_hpp

#include <string>
#include <vector>
#include <fstream>
#include <utility>

#include "t4c_defs.hpp"

// Four-center functions declaration generator for CPU.
class T4CFullDeclDriver
{
    /// Generates vector of variables strings.
    /// @param integral The base two center integral.
    /// @param terminus The flag to add termination symbol.
    /// @return The vector of variables strings.
    std::vector<std::string> _get_vars_str(const I4CIntegral& integral,
                                           const bool         terminus) const;
    
    /// Generates vector of variables strings.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param terminus The flag to add termination symbol.
    /// @return The vector of variables strings.
    std::vector<std::string> _get_prim_vars_str(const T4CIntegral& component,
                                                const I4CIntegral& integral,
                                                const bool         terminus) const;
    
    /// Generates vector of variables strings.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param terminus The flag to add termination symbol.
    /// @return The vector of variables strings.
    std::vector<std::string> _get_vrr_vars_str(const T4CIntegral& component,
                                               const I4CIntegral& integral,
                                               const bool         terminus) const;
    
public:
    /// Creates a four-center functions declaration generator.
    T4CFullDeclDriver() = default;
    
    /// Writes declaration for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param terminus The flag to add termination symbol.
    void write_func_decl(      std::ofstream& fstream,
                         const I4CIntegral&   integral,
                         const bool           terminus) const;
    
    /// Writes declaration of primitive compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param terminus The flag to add termination symbol.
    void write_prim_func_decl(      std::ofstream& fstream,
                              const T4CIntegral&   component,
                              const I4CIntegral&   integral,
                              const bool           terminus) const;
    
    /// Writes declaration of primitive VRR compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    /// @param terminus The flag to add termination symbol.
    void write_vrr_func_decl(      std::ofstream& fstream,
                              const T4CIntegral&   component,
                              const I4CIntegral&   integral,
                              const bool           terminus) const;
};

#endif /* t4c_full_decl_hpp */
