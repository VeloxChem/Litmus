// LITMUS: An Automated Molecular Integrals Generator
// Copyright 2022 Z. Rinkevicius, KTH, Sweden.
// E-mail: rinkevic@kth.se
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "t4c_geom_deriv_cpu_generators.hpp"

#include "t4c_utils.hpp"
#include "t4c_center_driver.hpp"

void
T4CGeomDerivCPUGenerator::generate(const int max_ang_mom,
                                   const int max_geom_order) const
{
    for (int i = 0; i <= max_ang_mom; i++)
    {
        for (int j = i; j <= max_ang_mom; j++)
        {
            for (int k = 0; k <= max_ang_mom; k++)
            {
                for (int l = k; l <= max_ang_mom; l++)
                {
                    for (int m = 0; m <= max_geom_order; m++)
                    {
                        for (int n = 0; n <= max_geom_order; n++)
                        {
                            for (int p = 0; p <= max_geom_order; p++)
                            {
                                for (int q = 0; q <= max_geom_order; q++)
                                {
                                    const auto gorder = m + n + p + q;
                                    
                                    if ((gorder > 0) && (gorder <= max_geom_order))
                                    {
                                        if (n > m) continue;
                                        
                                        if (q > p) continue;
                                        
                                        if (p > m) continue;
                                        
                                        if (q > n) continue;
                                        
                                        const auto integral = _get_integral({i, j, k, l}, {m, n, p, q});
                                        
                                        const auto components = integral.components<T2CPair, T2CPair>();
                                        
                                        const auto rgroup = _generate_integral_group(components, integral);
                                        
                                        _write_cpp_header(integral);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

I4CIntegral
T4CGeomDerivCPUGenerator::_get_integral(const std::array<int, 4>& ang_moms,
                                        const std::array<int, 4>& geom_drvs) const
{
    // bra and ket sides

    const auto bpair = I2CPair("GA", ang_moms[0], "GB", ang_moms[1]);

    const auto kpair = I2CPair("GC", ang_moms[2], "GD", ang_moms[3]);

    VOperators prefixes;

    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[0])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[1])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[2])));
    
    prefixes.push_back(Operator("d/dR", Tensor(geom_drvs[3])));

    return I4CIntegral(bpair, kpair, Operator("1"), 0, prefixes);
}

R4Group
T4CGeomDerivCPUGenerator::_generate_integral_group(const VT4CIntegrals& components,
                                                   const I4CIntegral&   integral) const
{
    R4Group rgroup;
        
    T4CCenterDriver t4c_geom_drv;
        
    return t4c_geom_drv.create_recursion(components);
}

void
T4CGeomDerivCPUGenerator::_write_cpp_header(const I4CIntegral& integral) const
{
    auto fname = t4c::geom_file_name(integral) + ".hpp";
        
    std::ofstream fstream;
               
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    //_write_hpp_defines(fstream, integral, true);
    
    //_write_hpp_includes(fstream, geom_integrals, vrr_integrals, integral);
    
    //_write_namespace(fstream, integral, true);
    
//    T4CDocuDriver docs_drv;
//
//    T4CDeclDriver decl_drv;
//
//    T4CFuncBodyDriver func_drv;
//
//    docs_drv.write_doc_str(fstream, integral, false);
//
//    decl_drv.write_func_decl(fstream, integral, false, false);
//
//    func_drv.write_geom_func_body(fstream, geom_integrals, vrr_integrals, integral);
//
//    fstream << std::endl;
//
//    _write_namespace(fstream, integral, false);
//
//    _write_hpp_defines(fstream, integral, false);
    
    fstream.close();
}
