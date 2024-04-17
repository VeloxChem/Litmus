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

namespace t2c { // t2c namespace

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
        return "NuclearPotentialGeom" + suffix;
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
    
    if (integrand.name() == "r")
    {
        const auto iorder = integrand.shape().order();
        
        if (iorder == 1)
        {
            return "Dipole";
        }
        
        if (iorder == 2)
        {
            return "Quadrupole";
        }
        
        if (iorder == 3)
        {
            return "Octupole";
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

std::string
namespace_label(const I2CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    const auto iorder = std::to_string(integrand.shape().order());
    
    if (integrand.name() == "AG")
    {
        return "npotg0" + iorder + "0rec";
    }
    
    if (integrand.name() == "A")
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

} // t2c namespace
