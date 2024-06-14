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

#include "t2c_utils.hpp"

#include <iostream>

#include "string_formater.hpp"

#include "v2i_ovl_driver.hpp"
#include "v2i_dip_driver.hpp"
#include "v2i_kin_driver.hpp"
#include "v2i_npot_driver.hpp"
#include "v2i_linmom_driver.hpp"
#include "v2i_el_field_driver.hpp"
#include "v2i_center_driver.hpp"

namespace t2c { // t2c namespace

// MR: Need to amend this for any new integral labels
std::string
integral_label(const I2CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    const auto prefixes = integral.prefixes();
    
    auto border = std::string("0");
    
    auto korder = std::string("0");
    
    const auto iorder = std::to_string(integrand.shape().order());
    
    if (const auto nterms = prefixes.size(); nterms > 0)
    {
        if (nterms >= 1) border = std::to_string(prefixes[0].shape().order());
        
        if (nterms >= 2) korder = std::to_string(prefixes[1].shape().order());
    }
    
    std::string suffix = "Geom" + border + iorder  + korder;
    
    if (integrand.name() == "AG")
    {
        return (prefixes.empty()) ? "NuclearPotentialGeom0" + iorder + "0" : "NuclearPotential" + suffix;
    }
    
    if (integrand.name() == "A")
    {
        return (prefixes.empty()) ? "NuclearPotential" : "NuclearPotential" + suffix;
    }
    
    if (integrand.name() == "T")
    {
        return (prefixes.empty()) ? "KineticEnergy" : "KineticEnergy" + suffix;
    }
    
    if (integrand.name() == "1")
    {
        return (prefixes.empty()) ? "Overlap" : "Overlap" + suffix;
    }

    if (integrand.name() == "p")
    {
        return (prefixes.empty()) ? "LinearMomentum" : "LinearMomentum" + suffix;
    }


    if (integrand.name() == "A1")
    {
        return (prefixes.empty()) ? "ElectricField" : "ElectricField" + suffix;
    }


    if (integrand.name() == "r")
    {
        const auto iorder = integrand.shape().order();
        
        if (iorder == 1)
        {
            return "ElectricDipoleMomentum";
        }
        
        if (iorder == 2)
        {
            return "ElectricQuadrupoleMomentum";
        }
        
        if (iorder == 3)
        {
            return "ElectricOctupoleMomentum";
        }
    }
    
    if (integrand.name() == "G(r)")
    {
        const auto iorder = integrand.shape().order();
        
        if (iorder == 0)
        {
            return "ThreeCenterOverlap";
        }
        
        if (iorder == 1)
        {
            return "ThreeCenterOverlapGradient";
        }
    }
    
    return std::string();
}

// MR: Need to amend this for any new integral labels
// Split label means that the string representation of the integral has got more than one parts and the parts are separated by underscores
std::string
integral_split_label(const I2CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    if ((integrand.name() == "A") || (integrand.name() == "AG"))
    {
        return "Nuclear_Potential";
    }
    
    if (integrand.name() == "T")
    {
        return "Kinetic_Energy";
    }
    
    if (integrand.name() == "1")
    {
        return "Overlap";
    }

    if (integrand.name() == "r")
    {
        return "Electric_Dipole_Momentum";
    }

    if (integrand.name() == "p")
    {
        return "Linear_Momentum";
    }

    if (integrand.name() == "A1")
    {
        return "Electric_Field";
    }

    return std::string();
}

// MR: Need to amend this for any new integral labels (this one is the name of the namespace to be used)
std::string
namespace_label(const I2CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    const auto iorder = std::to_string(integrand.shape().order());
    
    if ((integrand.name() == "A") || (integrand.name() == "AG"))
    {
        return "npotrec";
    }
    
    if (integrand.name() == "T")
    {
        return "kinrec";
    }
    
    if (integrand.name() == "1")
    {
        return "ovlrec";
    }

    if (integrand.name() == "p")
    {
        return "linmomrec";
    }


    if (integrand.name() == "A1")
    {
        return "elfield" + std::to_string(integrand.shape().order()) + "rec";
    }
    
    if (integrand.name() == "r")
    {
        if (iorder == "1") return "diprec";
        
        if (iorder == "2") return "quadrec";
        
        if (iorder == "3") return "octurec";
    }
    
    if (integrand.name() == "G(r)")
    {
        if (iorder == "0") return "t3ovlrec";
        
        if (iorder == "1") return "g3ovlrec";
    }
    
    return std::string();
}

// MR: Need to amend this for any new integral labels
// This one is the name of the operator when used in caption/documentation
std::string
integrand_label(const Operator& integrand)
{
    const auto iname = integrand.name();
    
    const auto iorder = std::to_string(integrand.shape().order());
    
    if (iname == "AG")
    {
        return iname + "(" + iorder + ")";
    }
    
    if ((iname == "r") && (iorder != "1"))
    {
        return iname + "^" + iorder;
    }

    if ((iname == "p") && (iorder != "1"))
    {
        return iname + "^" + iorder;
    }

    return iname;
}

// Only relevant for geometric derivatives
std::pair<std::string, std::string>
prefixes_label(const I2CIntegral& integral)
{
    const auto prefixes = integral.prefixes();
    
    auto bra_geom = std::string("");
    
    auto ket_geom = std::string("");
    
    if (const auto nterms = prefixes.size(); nterms > 0)
    {
        if (nterms >= 1)
        {
            const auto border = std::to_string(prefixes[0].shape().order());
            
            bra_geom = "d^(" + border + ")/dA^(" + border + ")";
        }
        
        if (nterms >= 2)
        {
            const auto korder = std::to_string(prefixes[1].shape().order());
            
            ket_geom = "d^(" + korder + ")/dB^(" + korder + ")";
        }
    }
    
    return std::make_pair(bra_geom, ket_geom);
}

// No changes ever needed here
std::vector<std::string>
integrand_labels(const I2CIntegral& integral,
                 const std::string& prefix)
{
    const auto op = integral.integrand();
    
    if (const auto op_comps = op.components(); op_comps.size() == 1)
    {
        return {prefix, };
    }
    else
    {
        std::vector<std::string> labels;
        
        for (const auto& op_comp : op_comps)
        {
            labels.push_back(prefix + "_" + op_comp.label());
        }
        
        return labels;
    }
}

std::string
compute_func_name(const I2CIntegral&           integral,
                  const std::pair<bool, bool>& rec_form)
{
    std::string prefix = (rec_form.first) ? "comp_sum_" : "comp_";
    
    std::string geom_label;
    
    auto tint_prefixes = integral.prefixes();
    
    if ((!tint_prefixes.empty()) || (integral.integrand().name() == "AG"))
    {
        geom_label += "_geom";
        
        if (integral.integrand().name() == "AG")
        {
            geom_label += "0" + std::to_string(integral.integrand().shape().order()) + "0";
        }
        else
        {
            for (const auto& tint_prefix : tint_prefixes)
            {
                geom_label += std::to_string(tint_prefix.shape().order());
            }
        }
    }
        
    auto label = prefix  + t2c::integral_split_label(integral) + geom_label + "_" + integral.label();
        
    return fstr::lowercase(label);
}

std::string
prim_file_name(const I2CIntegral& integral)
{
 if (integral.integrand().name() == "A1")
    {
        return t2c::integral_label(integral) + "_A" + std::to_string(integral.integrand().shape().order()) + "_" + "PrimRec" + integral.label();
    }

    return t2c::integral_label(integral) + "PrimRec" + integral.label();
}

// May need to amend this for new integral cases
std::string
get_buffer_label(const I2CIntegral& integral,
                 const std::string& prefix)
{
    std::string label = prefix + "_buffer_";
    
    if (integral.integrand().name() == "1") label += "ovl_";
    
    if (integral.integrand().name() == "T") label += "kin_";

    if (integral.integrand().name() == "r") label += "dip_";
    
    if (integral.integrand().name() == "p") label += "linmom_";

    if (integral.integrand().name() == "A1")
    {
        label += "el_field_A" + std::to_string(integral.integrand().shape().order()) + "_" + std::to_string(integral.order()) + "_";
    }

    if (integral.integrand().name() == "A")
    {
        label += "npot_" + std::to_string(integral.order()) + "_";
    }
    
    auto tint_prefixes = integral.prefixes();
    
    if (!tint_prefixes.empty())
    {
        label += "geom";
        
        for (const auto& tint_prefix : tint_prefixes)
        {
            label += std::to_string(tint_prefix.shape().order());
        }
        
       label += "_";
    }
    
    label += fstr::lowercase(integral.label());

    return label;
}

std::string
get_index_label(const I2CIntegral& integral)
{
    std::string label = "id_";
    
    if (integral.integrand().name() == "1") label += "ovl_";
    
    if (integral.integrand().name() == "T") label += "kin_";

    if (integral.integrand().name() == "r") label += "dip_";

    if (integral.integrand().name() == "p") label += "linmom_";

    if (integral.integrand().name() == "AG")
    {
        label += "npot_geom0" + std::to_string(integral.integrand().shape().order()) + "0_" + std::to_string(integral.order()) + "_";
    }

    if (integral.integrand().name() == "A")
    {
        label += "npot_" + std::to_string(integral.order()) + "_";
    }
    
    label += fstr::lowercase(integral.label());

    return label;
}

std::string
prim_compute_func_name(const I2CIntegral& integral)
{
    std::string geom_label;
    
    auto tint_prefixes = integral.prefixes();
    
    if ((!tint_prefixes.empty()) || (integral.integrand().name() == "AG"))
    {
        geom_label += "_geom";
        
        if (integral.integrand().name() == "AG")
        {
            geom_label += "0" + std::to_string(integral.integrand().shape().order()) + "0";
        }
        else
        {
            for (const auto& tint_prefix : tint_prefixes)
            {
                geom_label += std::to_string(tint_prefix.shape().order());
            }
        }
    }
    
    auto label =  "comp_prim_" + t2c::integral_split_label(integral) +  geom_label + "_" + integral.label();
    
    return fstr::lowercase(label);
}

// May need to amend this for new integral cases
SI2CIntegrals
get_integrals(const I2CIntegral& integral)
{
    SI2CIntegrals tints;
    
    if (!integral.is_simple())
    {
        V2ICenterDriver geom_drv;
        
        const auto prefixes = integral.prefixes();
        
        if ((prefixes.size() == 2) || (prefixes.size() == 1))
        {
            if ((integral.prefixes()[0].shape().order() == 0) &&
                (integral.prefixes()[1].shape().order() >  0))
            {
                return geom_drv.bra_ket_vrr(integral, 1);
            }
            else
            {
                return geom_drv.bra_ket_vrr(integral, 0);
            }
        }
    }
    
    if (integral.integrand().name() == "1")
    {
        V2IOverlapDriver ovl_drv;
        
        if (integral[0] > 0)
        {
            tints = ovl_drv.bra_vrr(integral);
        }
        else
        {
            tints = ovl_drv.ket_vrr(integral);
        }
    }
    
    if (integral.integrand().name() == "T")
    {
        V2IKineticEnergyDriver kin_drv;
        
        if (integral[0] > 0)
        {
            tints = kin_drv.bra_vrr(integral);
        }
        else
        {
            tints = kin_drv.ket_vrr(integral);
        }
        
        if ((integral[0] + integral[1]) == 0)
        {
            tints.insert(integral.replace(Operator("1"))); 
        }
    }
    
    if (integral.integrand().name() == "A")
    {
        V2INuclearPotentialDriver npot_drv;
        
        if (integral[0] > 0)
        {
            tints = npot_drv.bra_vrr(integral);
        }
        else
        {
            tints = npot_drv.ket_vrr(integral);
        }
        
        if ((integral[0] + integral[1]) == 0)
        {
            auto xint = integral.replace(Operator("1"));
            
            xint.set_order(0);
            
            tints.insert(xint);
        }
    }

    if (integral.integrand().name() == "r")
    {
        V2IDipoleDriver dip_drv;

        if (integral[0] > 0)
        {
            tints = dip_drv.bra_vrr(integral);
        }
        else
        {
            tints = dip_drv.ket_vrr(integral);
        }

        if ((integral[0] + integral[1]) == 0)
        {
            auto xint = integral.replace(Operator("1"));

            xint.set_order(0);

            tints.insert(xint);
        }
    }

    if (integral.integrand().name() == "p")
    {
        V2ILinearMomentumDriver linmom_drv;

        tints = linmom_drv.op_vrr(integral);

        if ((integral[0] + integral[1]) == 0)
        {
            auto xint = integral.replace(Operator("1"));

            xint.set_order(0);

            tints.insert(xint);
        }
    }

    if (integral.integrand().name() == "AG")
    {
        V2IElectricFieldDriver el_field_drv;

        if (integral[0] > 0)
        {
            tints = el_field_drv.bra_vrr(integral);
        }
        else
        {
            tints = el_field_drv.ket_vrr(integral);
        }

        if ((integral[0] + integral[1]) == 0)
        {
            tints = el_field_drv.aux_vrr(integral);
        }
    }


    return tints;
}

} // t2c namespace
