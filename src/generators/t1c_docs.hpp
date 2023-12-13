#ifndef t1c_docs_hpp
#define t1c_docs_hpp

#include <string>
#include <vector>
#include <fstream>

#include "t2c_defs.hpp"

// GTOs documentation generator for CPU.
class T1CDocuDriver
{
    /// Generates vector of variable strings.
    /// @return The vector of special variable strings.
    std::vector<std::string> _get_vars_str() const;
    
    /// Generates compute string.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    /// @return The  compute string.
    std::string _get_compute_str(const int angmom,
                                 const int gdrv) const;
    
public:
    /// Creates a GTOs documentation generator.
    T1CDocuDriver() = default;
    
    /// Writes documentation string for compute function.
    /// @param fstream the file stream.
    /// @param angmom The maximum angular momentum of GTOs.
    /// @param gdrv The geometrical derivative of GTOs.
    void write_doc_str(      std::ofstream& fstream,
                       const int            angmom,
                       const int            gdrv) const;
};


#endif /* t1c_docs_hpp */
