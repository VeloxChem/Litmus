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

#include "file_stream.hpp"

namespace ost { // ost namespace
    
    void
    write_copyright(std::ofstream& fstream)
    {
        fstream << "//                                                                              " << std::endl;
        
        fstream << "//                           VELOXCHEM 1.0-RC2                                  " << std::endl;
        
        fstream << "//         ----------------------------------------------------                 " << std::endl;
        
        fstream << "//                     An Electronic Structure Code                             " << std::endl;
        
        fstream << "//                                                                              " << std::endl;
        
        fstream << "//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.         " << std::endl;
        
        fstream << "//  Contact: https://veloxchem.org/contact                                      " << std::endl;
        
        fstream << "//                                                                              " << std::endl;
        
        fstream << "//  SPDX-License-Identifier: LGPL-3.0-or-later                                  " << std::endl;
        
        fstream << "//                                                                              " << std::endl;
        
        fstream << "//  This file is part of VeloxChem.                                             " << std::endl;
        
        fstream << "//                                                                              " << std::endl;
        
        fstream << "//  VeloxChem is free software: you can redistribute it and/or modify it under  " << std::endl;
        
        fstream << "//  the terms of the GNU Lesser General Public License as published by the Free " << std::endl;
        
        fstream << "//  Software Foundation, either version 3 of the License, or (at your option)   " << std::endl;
        
        fstream << "//  any later version.                                                          " << std::endl;
        
        fstream << "//                                                                              " << std::endl;
        
        fstream << "//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT    " << std::endl;
        
        fstream << "//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       " << std::endl;
        
        fstream << "//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public         " << std::endl;
        
        fstream << "//  License for more details.                                                   " << std::endl;
        
        fstream << "//                                                                              " << std::endl;
        
        fstream << "//  You should have received a copy of the GNU Lesser General Public License    " << std::endl;
        
        fstream << "//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.          " << std::endl;
        
        fstream << std::endl; 
    }

    void
    write_namespace(      std::ofstream& fstream,
                    const std::string&   label,
                    const bool           start)
    {
        if (start)
        {
            fstream << "namespace "  << label << " { // " << label << " namespace" << std::endl;
            
            fstream << std::endl;
        }
        else
        {
            fstream << std::endl;
            
            fstream << "} // " << label << " namespace" << std::endl;
        }
    }

    void
    write_vrr_includes(std::ofstream& fstream)
    {
        fstream << "#include <cstdint>" << std::endl << std::endl;
        
        fstream << "#include \"Buffer.hpp\"" << std::endl << std::endl;
    }

} // ost namespace

