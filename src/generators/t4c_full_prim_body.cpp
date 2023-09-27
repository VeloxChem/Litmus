#include "t4c_full_prim_body.hpp"

void
T4CFullPrimFuncBodyDriver::write_prim_func_body(      std::ofstream& fstream,
                                                const T4CIntegral&   component,
                                                const I4CIntegral&   integral) const
{
    auto lines = VCodeLines();
    
    lines.push_back({0, 0, 1, "{"});
    
    
    lines.push_back({0, 0, 2, "}"});
    
    ost::write_code_lines(fstream, lines);
}
