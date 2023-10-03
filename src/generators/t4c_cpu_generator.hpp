#ifndef t4c_cpu_generator_hpp
#define t4c_cpu_generator_hpp

#include <string>
#include <fstream>
#include <vector>
#include <map>

#include "operator.hpp"
#include "tensor_component.hpp"
#include "t4c_defs.hpp"
#include "file_stream.hpp"

// Four-center integrals code generator for CPU.
class T4CCPUGenerator
{
    /// Checks if recursion is available for four-center diagonal inetgral with given label.
    /// @param label The label of requested two-center integral.
    bool _is_available(const std::string& label) const;
    
    /// Gets four-center inetgral with requested label.
    /// @param label The label of requested two-center integral.
    /// @param ang_a The angular momentum of center A.
    /// @param ang_b The angular momentum of center B.
    /// @param ang_c The angular momentum of center C.
    /// @param ang_d The angular momentum of center D.
    /// @return The four-center integral.
    I4CIntegral _get_integral(const std::string& label,
                              const int          ang_a,
                              const int          ang_b,
                              const int          ang_c,
                              const int          ang_d) const;
    
    /// Gets file name of file with recursion functions for four center integral.
    /// @param integral The base two center integral.
    /// @return The file name.
    std::string _file_name(const I4CIntegral& integral) const;
    
    /// Writes header file for recursion.
    /// @param integral The base four center integral.
    void _write_cpp_header(const I4CIntegral& integral) const;
    
    /// Writes C++ code file for recursion.
    /// @param integral The base two center integral.
    void _write_cpp_file(const I4CIntegral& integral) const;
    
    /// Writes header files for primitive recursion.
    /// @param integral The base four center integral.
    void _write_cpp_prim_headers(const I4CIntegral& integral) const;
    
    /// Writes C++ code files for primitive recursion.
    /// @param integral The base four center integral.
    void _write_cpp_prim_files(const I4CIntegral& integral) const;
    
    /// Writes definitions of define for header file.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_defines(      std::ofstream& fstream,
                            const I4CIntegral&   integral,
                            const bool           start) const;
    
    /// Writes definitions of includes for header file.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void _write_hpp_includes(      std::ofstream& fstream,
                             const I4CIntegral&   integral) const;
    
    /// Writes definitions of includes for C++ code file.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void _write_cpp_includes(      std::ofstream& fstream,
                             const I4CIntegral&   integral) const;
    
    /// Writes namespace definition to file stream.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    /// @param start The flag to indicate position of namespace definition (start or end).
    void _write_namespace(      std::ofstream& fstream,
                          const I4CIntegral&   integral,
                          const bool           start) const;
    
    /// Writes definitions of define for primitive header file.
    /// @param fstream the file stream.
    /// @param fname The base file name.
    /// @param start The flag to indicate position of define (start or end).
    void _write_hpp_prim_defines(      std::ofstream& fstream,
                                 const std::string&   fname,
                                 const bool           start) const;
    
    /// Writes definitions of includes for primitives header file.
    /// @param fstream the file stream.
    /// @param integral The base four center integral.
    void _write_hpp_prim_includes(      std::ofstream& fstream,
                                  const I4CIntegral&   integral) const;
    
    /// Writes definitions of includes for C++ code file.
    /// @param fstream the file stream.
    /// @param component the integral component.
    /// @param integral The base two center integral.
    void _write_cpp_prim_includes(      std::ofstream& fstream,
                                  const T4CIntegral&   component,
                                  const I4CIntegral&   integral) const;
    
    /// Adds primitive functions headers to code lines container.
    /// @param lines The code lines container to which primitives loop start definition are added.
    /// @param integral The base four center integral.
    void _add_prim_call_includes(      VCodeLines&  lines,
                                 const I4CIntegral& integral) const;
    
public:
    /// Creates an electron repulsion integrals CPU code generator.
    T4CCPUGenerator() = default;
     
    /// Generates selected one-electron integrals up to given angular momentum (inclusive) )on A and B centers.
    /// @param label The label of requested two-center integral.
    /// @param angmom The maximum angular momentum of A and B centers.
    void generate(const std::string& label,
                  const int          angmom) const;
};


#endif /* t4c_cpu_generator_hpp */
