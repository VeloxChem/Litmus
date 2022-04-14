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

#include "eri_cpu_generator.hpp"

#include <fstream>

#include "file_stream.hpp"

EriCPUGenerator::EriCPUGenerator()

    : _diag_form(false)
{
    
}

void
EriCPUGenerator::set_diag_form()
{
    _diag_form = true;
}

void
EriCPUGenerator::generate(const Repository<R4Group, T4CIntegral>& repo) const
{
    for (const auto& tint : repo.base<I4CIntegral>())
    {
        _write_vrr_cpp_header(tint);
    }
}

std::string
EriCPUGenerator::_file_name(const I4CIntegral& integral,
                            const std::string& rectype) const
{
    std::string fname;
    
    if (integral.integrand() == Operator("1/|r-r'|")) fname += "Eri";
    
    if (_diag_form) fname += "Diag";
    
    fname += rectype + "For" + integral.label(); 
    
    return fname;
}

void
EriCPUGenerator::_write_vrr_cpp_header(const I4CIntegral& integral) const
{
    std::string fname = _file_name(integral, "VRR") + ".hpp";
        
    std::ofstream fstream;
           
    fstream.open(fname.c_str(), std::ios_base::trunc);
    
    ost::write_copyright(fstream);
    
//    _write_header_includes(fstream);
//
//    ost::write_namespace_start(_ref_namespace, fstream);
//
//    _write_header_data(fstream);
//
//    ost::write_namespace_end(_ref_namespace, fstream);
    
    fstream.close();
}

bool
EriCPUGenerator::is_hrr_rec_group(const R4Group& rgroup) const
{
    if (const auto facts = rgroup.factors(); facts.empty())
    {
        return false;
    }
    else
    {
        const std::set<Factor> hfacts({Factor("AB", "rab", TensorComponent(1, 0, 0)),
                                       Factor("AB", "rab", TensorComponent(0, 1, 0)),
                                       Factor("AB", "rab", TensorComponent(0, 0, 1)),
                                       Factor("CD", "rcd", TensorComponent(1, 0, 0)),
                                       Factor("CD", "rcd", TensorComponent(0, 1, 0)),
                                       Factor("CD", "rcd", TensorComponent(0, 0, 1))});
        
        for (const auto& fact : facts)
        {
            if (const auto idx = hfacts.find(fact) == hfacts.cend())
            {
                return false;
            }
        }
    }
        
    return true;
}


bool
EriCPUGenerator::is_vrr_rec_group(const R4Group& rgroup) const
{
    if (const auto facts = rgroup.factors(); facts.empty())
    {
        return false;
    }
    else
    {
        const std::set<Factor> vfacts({Factor("PB", "rpb", TensorComponent(1, 0, 0)),
                                       Factor("PB", "rpb", TensorComponent(0, 1, 0)),
                                       Factor("PB", "rpb", TensorComponent(0, 0, 1)),
                                       Factor("WP", "rwp", TensorComponent(1, 0, 0)),
                                       Factor("WP", "rwp", TensorComponent(0, 1, 0)),
                                       Factor("WP", "rwp", TensorComponent(0, 0, 1)),
                                       Factor("QD", "rqd", TensorComponent(1, 0, 0)),
                                       Factor("QD", "rqd", TensorComponent(0, 1, 0)),
                                       Factor("QD", "rqd", TensorComponent(0, 0, 1)),
                                       Factor("WQ", "rwq", TensorComponent(1, 0, 0)),
                                       Factor("WQ", "rwq", TensorComponent(0, 1, 0)),
                                       Factor("WQ", "rwq", TensorComponent(0, 0, 1)),
                                       Factor("1/zeta", "fz"),
                                       Factor("rho/zeta^2", "frz2"),
                                       Factor("1/(zeta+eta)", "fze"),
                                       Factor("1/eta", "fe"),
                                       Factor("rho/eta^2", "fre2")});
        
        for (const auto& fact : facts)
        {
            if (const auto idx = vfacts.find(fact) == vfacts.cend())
            {
                return false;
            }
        }
    }
    
    return true;
}

bool
EriCPUGenerator::is_aux_rec_group(const R4Group& rgroup) const
{
    if (const auto facts = rgroup.factors(); facts.empty())
    {
        return true;
    }
    else
    {
        return false;
    }
}
