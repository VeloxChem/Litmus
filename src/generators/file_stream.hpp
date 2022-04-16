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

namespace ost { // ost namespace
    
    /// Writes VeloxChem copyright titler to file stream.
    /// @param fstream the file stream.
    void write_copyright(std::ofstream& fstream);
    
} // ost namespace

#endif /* file_stream_hpp */