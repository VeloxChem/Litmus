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

#include "t4c_utils.hpp"

namespace t4c { // t2c namespace

std::string
integral_label(const I4CIntegral& integral)
{
    const auto integrand = integral.integrand();
    
    const auto prefixes = integral.prefixes();
    
    auto border = std::string("0");
    
    auto korder = std::string("0");
    
    const auto iorder = std::to_string(integrand.shape().order());
    
    if (const auto nterms = prefixes.size(); nterms > 0)
    {
        if (nterms >= 1) border =std::to_string(prefixes[0].shape().order());
        
        if (nterms >= 2) korder = std::to_string(prefixes[1].shape().order());
    }
    
    std::string suffix = "Geom" + border + iorder  + korder;
    
    if (integrand.name() == "1/|r-r'|")
    {
        return (prefixes.empty()) ? "ElectronRepulsion" : "ElectronRepulsion" + suffix;
    }
    
    return std::string();
}

} // t4c namespace
