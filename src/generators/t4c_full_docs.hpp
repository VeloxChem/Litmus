#ifndef t4c_full_docs_hpp
#define t4c_full_docs_hpp

#include <string>
#include <vector>
#include <fstream>

#include "t4c_defs.hpp"

// Four-center documentation generator for CPU.
class T4CFullDocuDriver
{
    /// Generates compute string.
    /// @param integral The base two center integral.
    /// @return The compute string.
    std::string _get_compute_str(const I4CIntegral& integral) const;
    
    /// Generates primitive compute string.
    /// @param component the integral component.
    /// @param integral The base four center integral.
    /// @return The primitive compute string.
    std::string _get_prim_compute_str(const T4CIntegral& component,
                                      const I4CIntegral& integral) const;
    
    /// Generates vector of variable strings.
    /// @return The vector of variable strings.
    std::vector<std::string> _get_vars_str() const;
    
    /// Generates vector of variable strings.
    /// @return The vector of variable strings.
    std::vector<std::string> _get_prim_vars_str() const;
    
    /// Generates vector of variable strings.
    /// @return The vector of variable strings.
    std::vector<std::string> _get_vrr_vars_str() const;

public:
    /// Creates a four-center documentation generator.
    T4CFullDocuDriver() = default;
    
    /// Writes documentation string for compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    void write_doc_str(      std::ofstream& fstream,
                       const I4CIntegral&   integral) const;
    
    /// Writes documentation string for primitive compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    void write_prim_doc_str(      std::ofstream& fstream,
                            const T4CIntegral&   component,
                            const I4CIntegral&   integral) const;
    
    /// Writes documentation string for primitive VRR compute function.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    void write_vrr_doc_str(      std::ofstream& fstream,
                           const T4CIntegral&   component,
                           const I4CIntegral&   integral) const;
};

#endif /* t4c_full_docs_hpp */
