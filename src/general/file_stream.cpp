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
    
    void write_code_lines(      std::ofstream& fstream, 
                          const VCodeLines&    lines)
    {
        for (const auto& [nspacers, offset, nends, str] : lines)
        {
            fstream << std::string(4 * nspacers + offset, ' ') << str;
            
            for (int i = 0; i < nends; i++) fstream << std::endl;
        }
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

} // ost namespace

