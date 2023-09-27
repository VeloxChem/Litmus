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
};


#endif /* t4c_full_prim_body_hpp */
