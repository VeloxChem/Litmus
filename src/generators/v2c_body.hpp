#ifndef v2c_body_hpp
#define v2c_body_hpp

#include <utility>
#include <string>
#include <vector>
#include <array>
#include <set>
#include <fstream>

#include "t2c_defs.hpp"
#include "t2c_utils.hpp"
#include "file_stream.hpp"

// Two-center compute function body generators for CPU.
class V2CFuncBodyDriver
{
    /// Generates vector of strings with GTOs definitions in compute function.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    /// @return The vector of strings with GTOS definitions in compute function.
    std::vector<std::string> _get_gtos_def(const bool diagonal) const;
    
    /// Generates vector of ket factors in compute function.
    /// @return The vector of ket factors in compute function.
    std::vector<std::string> _get_ket_variables_def(const bool diagonal) const;
    
    /// Generates vector of distances in compute function.
    /// @param integral The base two center integral.
    /// @return The vector of distances in compute function.
    std::vector<std::string> _get_coordinates_def(const I2CIntegral& integral) const;

public:
    /// Creates a two-center compute function body generator.
    V2CFuncBodyDriver() = default;
    
    /// Writes body of compute function.
    /// @param fstream the file stream.
    /// @param integral The base two center integral.
    /// @param sum_form The flag to used sum form for nuclear potential, multipoles, etc integrals.
    /// @param diagonal The flag to indicate diagonal or full form of compute function.
    void write_func_body(      std::ofstream& fstream,
                         const I2CIntegral&   integral,
                         const bool           sum_form,
                         const bool           diagonal) const;
};

#endif /* v2c_body_hpp */
