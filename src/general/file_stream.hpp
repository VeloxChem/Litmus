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

#ifndef file_stream_hpp
#define file_stream_hpp

#include <fstream>
#include <string>
#include <tuple>
#include <vector>

using TCodeLine = std::tuple<int, int, int, std::string>;

using VCodeLines = std::vector<TCodeLine>;

namespace ost { // ost namespace
    
    /// Writes vector of code lines to file stream.
    /// @param lines the vector of code lines.
    /// @param fstream the file stream.
    void write_code_lines(      std::ofstream& fstream,
                          const VCodeLines&    lines);

    /// Writes VeloxChem copyright titler to file stream.
    /// @param fstream the file stream.
    void write_copyright(std::ofstream& fstream);

    /// Writes VeloxChem copyright titler to file stream.
    /// @param fstream the file stream.
    /// @param label the namespace label
    /// @param start the flag to write name space start or end.
    void write_namespace(      std::ofstream& fstream,
                         const std::string&   label,
                         const bool           start);

    /// Writes HRR/VRR  includes to file stream.
    /// @param fstream the file stream.
    void write_hvrr_includes(std::ofstream& fstream);

    /// Writes recursion dimensions to file stream.
    /// @param fstream the file stream.
    void write_dimensions(std::ofstream& fstream);

} // ost namespace

#endif /* file_stream_hpp */
