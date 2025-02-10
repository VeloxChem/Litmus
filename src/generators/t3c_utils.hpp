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

#ifndef t3c_utils_hpp
#define t3c_utils_hpp

#include <string>
#include <array>

#include "t3c_defs.hpp"

namespace t3c { // t3c namespace

/// Gets standart capitalized label of integral.
/// @param integral The base two center integral.
/// @return The standart capitalized label of integral.
std::string integral_label(const I3CIntegral& integral);

/// Gets standart split label of integral.
/// @param integral The base two center integral.
/// @return The standart split label of integral.
std::string integral_split_label(const I3CIntegral& integral);

/// Generates primitive file name.
/// @param integral The base four center integral.
/// @return The primitive file name.
std::string prim_file_name(const I3CIntegral& integral);

/// Generates ket horizontal recursion file name.
/// @param integral The base four center integral.
/// @return The primitive file name.
std::string hrr_file_name(const I3CIntegral& integral);

/// Gets recursion namespace label of standart integral.
/// @param integral The base two center integral.
/// @return The recursion namespace label of standart integral.
std::string namespace_label(const I3CIntegral& integral);

/// Gets label of standart integrand.
/// @param integrand the integrand operator.
/// @return The label of standart integrand.
std::string integrand_label(const Operator& integrand);

/// Generates compute function  name.
/// @param integral The base four center integral.
/// @return The compute function name.
std::string compute_func_name(const I3CIntegral& integral);

/// Generates integral buffer label.
/// @param integral The base two center integral.
/// @return The string with integral label.
std::string get_buffer_label(const I3CIntegral& integral,
                             const std::string& prefix);

/// Generates primitive compute function name.
/// @param integral The base two center integral.
/// @return The primitive compute function name.
std::string prim_compute_func_name(const I3CIntegral& integral);

/// Gets arguments list for primitive vertical recursion function call.
/// @param integral The base two center integral.
SI3CIntegrals get_vrr_integrals(const I3CIntegral& integral);

/// Generates ket horizontal recursion compute function name.
/// @param integral The base two center integral.
/// @return The primitive compute function name.
std::string hrr_compute_func_name(const I3CIntegral& integral);

/// Gets arguments list for ket horizontal recursion function call.
/// @param integral The base two center integral.
SI3CIntegrals get_hrr_integrals(const I3CIntegral& integral);

/// Generates integral index label.
/// @param integral The base two center integral.
/// @return The string with index label.
std::string get_index_label(const I3CIntegral& integral);

/// Gets arguments list for primitive vertical recursion function call.
/// @param integral The base two center integral.
SI3CIntegrals get_vrr_integrals(const I3CIntegral& integral);

/// Generates horizontal recursion integral buffer index.
/// @param integral The base two center integral.
/// @return The string with integral buffer index.
std::string get_hrr_index(const I3CIntegral& integral);

/// Generates horizontal recursion integral buffer label.
/// @param integral The base two center integral.
/// @return The string with integral label.
std::string get_hrr_buffer_label(const I3CIntegral& integral);

/// Generates geometrical derrivatives labels.
/// @param integral The base four center integral.
/// @return The geometrical derivative labels.
std::string prefixes_label(const I3CIntegral& integral);

/// Prunes geometrical recursion term.
/// @param term The geometrical recursion term.
/// @return The pruned geometrical recursion term.
G3Term prune_term(const G3Term& term);

} // t3c namespace

#endif /* t3c_utils_hpp */
